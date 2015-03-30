/*
 * Copyright (c) 2015, Logan P. Evans <loganpevans@gmail.com>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "decaying_histogram.h"
#include <random>

#define NUM_BUCKETS 30
#define NUM_THREADS 1
#define ALPHA 0.00001
#define CYCLES (1ULL * 1024 * 1024 * 1024)
#define FRAMES_PER_SECOND 5
#define PERMITTED_BYTES (1ULL * 1024 * 1024 * 1024)
#define BYTES_PER_MALLOC (1024ULL * 1024)
#define MAX_NUM_THREADS 200
#define CONSTANT_SPACE_PER_TRIAL 0

static struct decaying_histogram *g_histogram;

static uint64_t rdtsc()
{
  uint32_t hi, lo;

  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return (((uint64_t)lo) | (((uint64_t)hi) << 32));
}

struct thread_func_args {
  uint64_t end_timestamp;
  uint64_t permitted_bytes;
  uint64_t bytes_per_malloc;
};

static void *
thread_func(void *args) {
  uint64_t start_timestamp, stop_timestamp, end_timestamp;
  uint64_t permitted_bytes;
  uint64_t bytes_per_malloc;
  char **buffer;
  int buffer_size;
  int idx;
  int malloc_idx;
  std::default_random_engine generator;
  std::uniform_int_distribution<int> *distribution;

  end_timestamp = ((struct thread_func_args *)args)->end_timestamp;
  permitted_bytes = ((struct thread_func_args *)args)->permitted_bytes;
  bytes_per_malloc = ((struct thread_func_args *)args)->bytes_per_malloc;

  buffer_size = permitted_bytes / bytes_per_malloc;

  distribution = new std::uniform_int_distribution<int>(0, buffer_size - 1);

  buffer = (char **)malloc(sizeof(char*) * buffer_size);
  for (idx = 0; idx < buffer_size; idx++) {
    buffer[idx] = NULL;
  }

  while (rdtsc() < end_timestamp) {
    malloc_idx = (*distribution)(generator);
    if (buffer[malloc_idx] == NULL) {
      start_timestamp = rdtsc();
      buffer[malloc_idx] = (char *)malloc(bytes_per_malloc);
      stop_timestamp = rdtsc();
      dh_insert(
          g_histogram, log2(stop_timestamp - start_timestamp),
          DHIST_MULTI_THREADED);
    } else {
      free(buffer[malloc_idx]);
      buffer[malloc_idx] = NULL;
    }
  }

  for (idx = 0; idx < buffer_size; idx++) {
    free(buffer[idx]);
  }

  delete distribution;
  free(buffer);
  return (void *)NULL;
}

int main() {
  pthread_t threads[MAX_NUM_THREADS];
  struct thread_func_args args;
  struct timespec tim, tim2;
  char title[1000];
  tim.tv_sec = 0;
  tim.tv_nsec = 1000000000 / FRAMES_PER_SECOND;

  args.permitted_bytes = PERMITTED_BYTES / MAX_NUM_THREADS;
  args.bytes_per_malloc = BYTES_PER_MALLOC;

  g_histogram = new struct decaying_histogram;
  init_decaying_histogram(g_histogram, NUM_BUCKETS, ALPHA);

  for (int num_threads = 1; num_threads < MAX_NUM_THREADS; num_threads++) {
    sprintf(title, "threads: %d", num_threads);
    args.end_timestamp = rdtsc() + CYCLES;
    if (CONSTANT_SPACE_PER_TRIAL)
      args.permitted_bytes = PERMITTED_BYTES / num_threads;

    for (int i = 0; i < num_threads; i++)
      pthread_create(
        &threads[i], NULL, (void *(*)(void *))thread_func, &args);

    while (rdtsc() < args.end_timestamp) {
      nanosleep(&tim , &tim2);
      print_histogram(g_histogram, true, title, "log_2(insertion time)",
          DHIST_MULTI_THREADED);
      fflush(stdout);
    }

    for (int i = 0; i < NUM_THREADS; i++)
      pthread_join(threads[i], NULL);
  }

  //clean_decaying_histogram(g_histogram);
  //delete g_histogram;

  return 0;
}

