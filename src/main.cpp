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
#include "decaying_histogram.h"
#include <random>

#define NUM_BUCKETS 100
#define NUM_THREADS 50
#define ALPHA 0.000001
#define CYCLES 16ULL * 1024 * 1024 * 1024

struct decaying_histogram *g_histogram;

uint64_t rdtsc()
{
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
  while (last_timestamp < stop_timestamp) {
    this_timestamp = rdtsc();
    add_observation(g_histogram, log(this_timestamp - last_timestamp));
    last_timestamp = this_timestamp;
  }

  return (void *)NULL;
}

int main() {
  pthread_t threads[NUM_THREADS];
  struct thread_func_args args;

  args.stop_timestamp = rdtsc() + CYCLES;

  g_histogram = new struct decaying_histogram;
  init_decaying_histogram(g_histogram, NUM_BUCKETS, ALPHA);

  for (int i = 0; i < NUM_THREADS; i++)
    pthread_create(
      &threads[i], NULL, (void *(*)(void *))thread_func, &args);

  for (int i = 0; i < NUM_THREADS; i++)
    pthread_join(threads[i], NULL);

  print_histogram(g_histogram);
  clean_decaying_histogram(g_histogram);
  delete g_histogram;

  return 0;
}

