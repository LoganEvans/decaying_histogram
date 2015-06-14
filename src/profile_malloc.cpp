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


#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <random>
#include "dhist.h"

#define SUCCESS 9001
#define FAILURE 9000
#define ANIMATE 1
#define PROFILE_FREE 0

uint64_t NUM_BUCKETS = 30;
double DECAY_RATE = 0.99999;
uint64_t FRAMES_PER_SECOND = 5;

uint64_t START_NUM_THREADS = 1;
uint64_t MAX_NUM_THREADS = 255;
// 2.397 * 2**30 is roughly one second.
uint64_t CYCLES_PER_TRIAL = 5ULL * 2.397 * 1024 * 1024 * 1024;
uint64_t PERMITTED_BYTES = 1ULL * 1024 * 1024 * 1024;

bool CONSTANT_SPACE_PER_TRIAL[] = {false};

uint64_t BYTES_PER_MALLOC[] = {
  8ULL,  // tiny
  512ULL,  // quantum-spaced
  1024ULL,  // sub-page
  4ULL * 1024,  // large
  1ULL * 1024 * 1024,  // large
  2ULL * 1024 * 1024,  // huge
};

static struct dhist *g_histogram;

static uint64_t rdtsc()
{
  uint32_t hi, lo;

  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return (((uint64_t)lo) | (((uint64_t)hi) << 32));
}

struct thread_func_args {
  uint64_t end_timestamp;
  uint64_t bytes_per_malloc;
  uint64_t bytes_per_thread;
  uint64_t num_mallocs_out;
  uint64_t total_cycles_out;
};

static void *
thread_func(void *args) {
  uint64_t start_timestamp, stop_timestamp, end_timestamp;
  uint64_t bytes_per_malloc;
  uint64_t bytes_per_thread;
  uint64_t num_mallocs_out;
  uint64_t total_cycles_out;
  char **buffer;
  int buffer_size;
  int idx;
  int malloc_idx;
  std::default_random_engine generator;
  std::uniform_int_distribution<int> *distribution;

  end_timestamp = ((struct thread_func_args *)args)->end_timestamp;
  bytes_per_thread = ((struct thread_func_args *)args)->bytes_per_thread;
  bytes_per_malloc = ((struct thread_func_args *)args)->bytes_per_malloc;
  num_mallocs_out = 0;
  total_cycles_out = 0;

  buffer_size = bytes_per_thread / bytes_per_malloc;

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
#if ANIMATE && !PROFILE_FREE
      dhist_insert(
          g_histogram, log2(stop_timestamp - start_timestamp),
          DHIST_MULTI_THREADED);
#endif  // ANIMATE
      num_mallocs_out++;
      total_cycles_out += stop_timestamp - start_timestamp;
    } else {
      start_timestamp = rdtsc();
      free(buffer[malloc_idx]);
      stop_timestamp = rdtsc();
      buffer[malloc_idx] = NULL;

#if ANIMATE && PROFILE_FREE
      dhist_insert(
          g_histogram, log2(stop_timestamp - start_timestamp),
          DHIST_MULTI_THREADED);
#endif
    }
  }

  for (idx = 0; idx < buffer_size; idx++) {
    free(buffer[idx]);
  }

  ((struct thread_func_args *)args)->num_mallocs_out = num_mallocs_out;
  ((struct thread_func_args *)args)->total_cycles_out = total_cycles_out;

  delete distribution;
  free(buffer);
  return (void *)SUCCESS;
}

int main() {
  pthread_t threads[MAX_NUM_THREADS];
  struct thread_func_args args[MAX_NUM_THREADS];
  struct timespec tim, tim2;
  char title[1000];
  char *histogram_json;
  uint64_t end_timestamp;
  uint64_t bytes_per_malloc;
  uint64_t bytes_per_thread;
  uint64_t num_mallocs;
  uint64_t total_cycles;
  void *status;
  bool constant_space_per_trial;

  g_histogram = dhist_init(NUM_BUCKETS, DECAY_RATE);
  tim.tv_sec = 0;
  tim.tv_nsec = 1000000000 / FRAMES_PER_SECOND;

#if !ANIMATE
  printf(
      "threads,cycles_per_malloc,num_mallocs,"
      "bytes_per_malloc,bytes_per_thread,constant_space_per_trial\n");
#endif

  for (int space_idx = 0;
      space_idx < sizeof(CONSTANT_SPACE_PER_TRIAL) /
                  sizeof(CONSTANT_SPACE_PER_TRIAL[0]);
      space_idx++) {
    constant_space_per_trial = CONSTANT_SPACE_PER_TRIAL[space_idx];
    for (int bytes_idx = 0;
         bytes_idx < sizeof(BYTES_PER_MALLOC) / sizeof(BYTES_PER_MALLOC[0]);
         bytes_idx++) {
      bytes_per_malloc = BYTES_PER_MALLOC[bytes_idx];

      for (int num_threads = START_NUM_THREADS;
           num_threads < MAX_NUM_THREADS; num_threads++) {
        sprintf(title,
            "threads: %d, bytes-per-malloc: %lu, bytes-per-thread: %lu",
            num_threads, bytes_per_malloc, bytes_per_thread);
        end_timestamp = rdtsc() + CYCLES_PER_TRIAL;
        if (constant_space_per_trial)
          bytes_per_thread = PERMITTED_BYTES / num_threads;
        else
          bytes_per_thread = PERMITTED_BYTES / MAX_NUM_THREADS;

        for (int i = 0; i < num_threads; i++) {
          args[i].end_timestamp = end_timestamp;
          args[i].bytes_per_thread = bytes_per_thread;
          args[i].bytes_per_malloc = bytes_per_malloc;
          pthread_create(
            &threads[i], NULL, (void *(*)(void *))thread_func, &args[i]);
        }

#if ANIMATE
        while (rdtsc() < end_timestamp) {
          nanosleep(&tim , &tim2);
          histogram_json = dhist_get_json(
              g_histogram, title, "log_2(insertion time)",
              DHIST_MULTI_THREADED);
          puts(histogram_json);
          free(histogram_json);
          fflush(stdout);
        }
#endif

        num_mallocs = 0;
        total_cycles = 0;
        for (int i = 0; i < num_threads; i++) {
          pthread_join(threads[i], &status);
          if ((uint64_t)status != SUCCESS) {
            printf("Thread %d returned status: %lu\n", i, (uint64_t)status);
            assert(false);
          }
          num_mallocs += args[i].num_mallocs_out;
          total_cycles += args[i].total_cycles_out;
        }
#if !ANIMATE
        printf("%d,%lu,%lu,%lu,%lu,%d\n",
            num_threads, total_cycles / num_mallocs, num_mallocs,
            bytes_per_malloc, bytes_per_thread, constant_space_per_trial);
#endif
      }
    }
  }

  dhist_destroy(g_histogram);

  return 0;
}

