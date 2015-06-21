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
#include <pthread.h>
#include "dhist.h"
#include <random>

#define NUM_BUCKETS 500
#define NUM_THREADS 4
#define DECAY_RATE 0.9999

// 2.397 * 1024 * 1024 * 1024 is roughly the number of cycles in a second.
#define CYCLES ((uint64_t)(10 * 2.397 * 1024 * 1024 * 1024))
#define DHIST_MP_FLAG \
    (NUM_THREADS > 1 ? DHIST_MULTI_THREADED : DHIST_SINGLE_THREADED)
#define FRAMES_PER_SECOND 5.01

#define ANIMATE 0
#define NORMAL_DISTRIBUTION 0

static struct dhist *g_histogram;

#if NORMAL_DISTRIBUTION
std::default_random_engine g_generator;
std::normal_distribution<double> g_distribution(0.0, 1.0);
#endif

static uint64_t rdtsc() {
  uint32_t hi, lo;

  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return (((uint64_t)lo) | (((uint64_t)hi) << 32));
}

struct thread_func_args {
  uint64_t stop_timestamp;
};

static void *
thread_func(void *args) {
  uint64_t last_timestamp, this_timestamp, stop_timestamp;

  stop_timestamp = ((struct thread_func_args *)args)->stop_timestamp;
  last_timestamp = rdtsc();
#if ANIMATE
  while (1) {
#else
  while (last_timestamp < stop_timestamp) {
#endif
    this_timestamp = rdtsc();
#if NORMAL_DISTRIBUTION
    dhist_insert(
        g_histogram, g_distribution(g_generator), DHIST_MP_FLAG);
#else
    dhist_insert(
        g_histogram, log2(this_timestamp - last_timestamp), DHIST_MP_FLAG);
#endif
    last_timestamp = this_timestamp;
  }

  return (void *)NULL;
}

int main() {
  pthread_t threads[NUM_THREADS];
  struct thread_func_args args;
  struct timespec tim, tim2;
  char *histogram_json = NULL;
  int json_len, needed_len;

  json_len = 0;
  tim.tv_sec = 0;
  tim.tv_nsec = 1000000000 / FRAMES_PER_SECOND;

  args.stop_timestamp = rdtsc() + CYCLES;

  g_histogram = dhist_init(NUM_BUCKETS, DECAY_RATE);

  for (int i = 0; i < NUM_THREADS; i++)
    pthread_create(
      &threads[i], NULL, (void *(*)(void *))thread_func, &args);

#if ANIMATE
  while (1) {
    nanosleep(&tim , &tim2);
    while (1) {
      needed_len = dhist_snprint_histogram(histogram_json, json_len,
          g_histogram, "Test", "log_2(insertion time)", DHIST_MP_FLAG);
      if (needed_len >= json_len) {
        json_len = needed_len * 2 + 1;
        histogram_json = (char *)realloc(histogram_json, json_len);
      } else {
        break;
      }
    }

    puts(histogram_json);
    fflush(stdout);
  }

  free(histogram_json);
#endif

  for (int i = 0; i < NUM_THREADS; i++)
    pthread_join(threads[i], NULL);

  histogram_json = dhist_get_json(
      g_histogram, "Test", "log_2(insertion time)", DHIST_MP_FLAG);
  puts(histogram_json);
  free(histogram_json);
  dhist_destroy(g_histogram);

  return 0;
}

