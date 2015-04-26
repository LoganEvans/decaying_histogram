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
#include "decaying_histogram.h"
#include <random>

#define NUM_BUCKETS 50
#define ALPHA_SLOW 0.0001
#define ALPHA_FAST 0.000417128920021

#define COUNT 10000000
//#define COUNT 100

#define ANIMATE 0

int main() {
  int idx;
  char *histogram_json;
  struct decaying_histogram *dhist_slow, *dhist_fast;
  dhist_slow = dhist_fast = NULL;

  std::default_random_engine generator;
  std::normal_distribution<double> normal_0_1(0.0, 1.0);
  std::exponential_distribution<double> exponential_1(1.0);

  dhist_slow = new struct decaying_histogram;
  init_decaying_histogram(dhist_slow, NUM_BUCKETS, ALPHA_SLOW);

  dhist_fast = new struct decaying_histogram;
  init_decaying_histogram(dhist_fast, NUM_BUCKETS, ALPHA_FAST);

  uint64_t last_timestamp, this_timestamp;
  double observation;

  int iterations = 2;

  while (iterations--) {
    for (idx = 0; idx < COUNT; idx++) {
      if (idx % (COUNT / 100) == 0)
        fprintf(stderr, "%d / %d              \r", idx, COUNT);
      observation = normal_0_1(generator);
      //observation = exponential_1(generator);
      dhist_insert(dhist_slow, observation, DHIST_SINGLE_THREADED);
      dhist_insert(dhist_fast, observation, DHIST_SINGLE_THREADED);

#if 0
      if (idx % 10000 == 0) {
        histogram_json = get_new_histogram_json(
              dhist_slow, true, "dhist_slow", NULL, DHIST_SINGLE_THREADED);
        puts(histogram_json);
        free(histogram_json);

        histogram_json = get_new_histogram_json(
              dhist_fast, true, "dhist_fast", NULL, DHIST_SINGLE_THREADED);
        puts(histogram_json);
        free(histogram_json);
      }
#endif

#if 1
      if (iterations == 0) {
        printf("%lf\n",
            dhist_distance(
              //dhist_earth_movers_distance, dhist_slow, dhist_fast, true,
              dhist_Kolmogorov_Smirnov_statistic, dhist_slow, dhist_fast, true,
              //dhist_Jaccard_distance, dhist_slow, dhist_fast, true,
              DHIST_SINGLE_THREADED));
      }
#endif
    }
  }

  return 0;
}

