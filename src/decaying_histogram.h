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

struct bucket {
  // The values for count, mu, update_generation and the below and above
  // pointers are protected by bucket->boundary_mtx and
  // bucket->above->boundary_mtx. If this thread holds either of those mutexes,
  // these values cannot be modified by any other thread.
  double count;
  double mu;
  uint64_t update_generation;
  struct bucket *below;
  struct bucket *above;
  // The parent and children pointers, height, and is_enabled fields are
  // protected by tree_mtx.
  struct bucket *parent;
  struct bucket *children[2];
  pthread_mutex_t *boundary_mtx;  /* lower boundary */
  int height;
  // The is_enabled field exists because the bucket lookup is not protected with
  // any atomics, so it is possible for a thread to obtain a reference to this
  // bucket while it is enabled, for another thread to recycle the bucket, and
  // then for the first thread to attempt to update the bucket.
  bool is_enabled;
  bool lock_held;
};

struct decaying_histogram {
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

struct bucket * init_bucket(struct decaying_histogram *histogram, int mp_flag);
void init_decaying_histogram(
    struct decaying_histogram *histogram, int target_buckets,
    double alpha);
void clean_decaying_histogram(struct decaying_histogram *histogram);
void dh_insert(
    struct decaying_histogram *histogram, double observation, int mp_flag);
/*
 * This returns the bucket with the greatest mu less than observation (which is
 * not necessarily the bucket where observation will be inserted) unless the
 * target bucket is the leftmost bucket, in which case observation may be less
 * than bucket->mu.
 * If mp_flag & DHIST_MULTI_THREADED, the boundary_mtx for this bucket and its
 * right hand neighbor will be locked.
 */
struct bucket * find_bucket(
    struct decaying_histogram *histogram, double observation, int mp_flag);
double total_count(struct decaying_histogram *histogram, uint64_t generation);
double weight(
    struct decaying_histogram *histogram, struct bucket *bucket,
    uint64_t generation,
    double *lower_bound_output, double *upper_bound_output);
void decay(
    struct decaying_histogram *histogram, struct bucket *bucket,
    uint64_t generation);
double Jaccard_distance(
    struct decaying_histogram *hist0, struct decaying_histogram *hist1);
double Kolomogorov_Smirnov_statistic(
    struct decaying_histogram *hist0, struct decaying_histogram *hist1);
void print_histogram(
    struct decaying_histogram *histogram, bool estimate_ok,
    const char *title, const char *xaxis, int mp_flag);

int assert_invariant(struct bucket *root);
int count_nodes(struct bucket *root);
void split_bucket(
    struct decaying_histogram *histogram, struct bucket *bucket, int mp_flag);
void delete_bucket(
    struct decaying_histogram *histogram, struct bucket *bucket, int mp_flag);
void print_tree(struct bucket *bucket);

#ifdef __cplusplus
}
#endif

#endif  /* DECAYING_HISTOGRAM_H_ */

