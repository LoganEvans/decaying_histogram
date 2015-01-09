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


#include <stdio.h>
#include <pthread.h>
#include "decaying_histogram.h"
#include <random>

#define NUM_BUCKETS 50
#define NUM_THREADS 2
#define ALPHA 0.0001
#define OBSERVATIONS 100000

struct decaying_histogram *g_histogram;

static void *
thread_func(void *args) {
  std::default_random_engine generator;
  std::normal_distribution<double> normal;

  for (int i = 0; i < OBSERVATIONS; i++) {
    add_observation(g_histogram, normal(generator));
  }

  return (void *)NULL;
}

int main() {
  pthread_t threads[NUM_BUCKETS];

  g_histogram = new struct decaying_histogram;
  init_decaying_histogram(g_histogram, NUM_BUCKETS, ALPHA);

  for (int i = 0; i < NUM_THREADS; i++)
    pthread_create(&threads[i], NULL, (void *(*)(void *))thread_func, NULL);

  for (int i = 0; i < NUM_THREADS; i++)
    pthread_join(threads[i], NULL);

  print_histogram(g_histogram);
  clean_decaying_histogram(g_histogram);
  delete g_histogram;

  return 0;
}

