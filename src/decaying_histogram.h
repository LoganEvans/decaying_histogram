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


#ifndef DECAYING_HISTOGRAM_H_
#define DECAYING_HISTOGRAM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdbool.h>
#include <pthread.h>

struct bucket {
  double count;
  double mu;
  double lower_bound;
  double upper_bound;
  uint64_t last_decay_generation;
  struct bucket *below;
  struct bucket *above;
  pthread_mutex_t mutex;
  pthread_mutexattr_t mutexattr;
  char __padding[4];
};

struct decaying_histogram {
  double delete_bucket_threshold;
  double split_bucket_threshold;
  double alpha;
  uint64_t generation;
  struct bucket *bucket_list;
  int num_buckets;
  unsigned int max_num_buckets;
  double *pow_table;
  pthread_rwlock_t rwlock;
  pthread_mutex_t generation_mutex;
};

#define ABS_DIFF(x, y) (((x) - (y)) > 0 ? (x) - (y) : (y) - (x))
#define THRESH 0.00004

void init_bucket(
    struct bucket *to_init, struct bucket *below, struct bucket *above);
void init_decaying_histogram(
    struct decaying_histogram *histogram, int target_buckets,
    double alpha);
void clean_decaying_histogram(struct decaying_histogram *histogram);
void add_observation(struct decaying_histogram *histogram, double observation);
int find_bucket_idx(struct decaying_histogram *histogram, double observation);
double total_count(struct decaying_histogram *histogram);
// If lower_bound or upper_bound are non-NULL, this will store the lower_bound
// and upper bound at the provided address.
double density(
    struct decaying_histogram *histogram, struct bucket *bucket,
    double *lower_bound_output, double *upper_bound_output);
double Jaccard_distance(
    struct decaying_histogram *hist0, struct decaying_histogram *hist1);
double Kolomogorov_Smirnov_statistic(
    struct decaying_histogram *hist0, struct decaying_histogram *hist1);
void print_histogram(struct decaying_histogram *histogram);
void full_refresh(struct decaying_histogram *histogram);

#ifdef __cplusplus
}
#endif

#endif  /* DECAYING_HISTOGRAM_H_ */

