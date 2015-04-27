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

extern const int DHIST_SINGLE_THREADED;
extern const int DHIST_MULTI_THREADED;

enum dhist_distance_t {
  dhist_Jaccard_distance,
  dhist_Kolmogorov_Smirnov_statistic,
  dhist_earth_movers_distance
};

struct dhist {
  double delete_bucket_threshold;
  double split_bucket_threshold;
  double alpha;
  uint64_t generation;
  struct bucket *root;
  // The bucket_list is a pool of buckets that are allocated at initialization
  // time. The length will be max_num_buckets.
  struct bucket *bucket_list;
  uint32_t num_buckets;
  uint32_t max_num_buckets;
  double *pow_table;
  pthread_mutex_t *tree_mtx;
  pthread_mutex_t *generation_mtx;
};

struct dhist * dhist_init(int target_buckets, double alpha);
void dhist_destroy(struct dhist *histogram);
void dhist_insert(struct dhist *histogram, double observation, int mp_flag);
char * dhist_get_json(
    struct dhist *histogram, bool estimate_ok,
    const char *title, const char *xlabel, int mp_flag);
void dhist_set_target_buckets(struct dhist *histogram, int target_buckets);
void dhist_set_alpha(struct dhist *histogram, double alpha);

double dhist_distance(
    enum dhist_distance_t distance_name,
    struct dhist *hist1, struct dhist *hist2,
    bool estimate_ok, int mp_flag);

#ifdef __cplusplus
}
#endif

#endif  /* DECAYING_HISTOGRAM_H_ */

