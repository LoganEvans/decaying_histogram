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


#ifndef DECAYING_HISTOGRAM_SRC_DHIST_H_
#define DECAYING_HISTOGRAM_SRC_DHIST_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdbool.h>
#include <pthread.h>

extern const int DHIST_SINGLE_THREADED;
extern const int DHIST_MULTI_THREADED;

struct dhist {
  double decay_rate;
  double total_count;
  uint64_t generation;
  struct bucket *root;
  uint32_t num_buckets;
  uint32_t target_num_buckets;
  double *pow_table;
  pthread_mutex_t *tree_mtx;
  pthread_mutex_t *generation_mtx;
  pthread_mutex_t *thread_info_mtx;
  struct bucket *fix_balance_stack;
  struct thread_info *thread_info_head;
  uint32_t num_precomputed_powers;
};

struct dhist * dhist_init(uint32_t target_buckets, double decay_rate);
void dhist_destroy(struct dhist *histogram);
void dhist_insert(struct dhist *histogram, double observation, int mp_flag);
int dhist_snprint_histogram(
    char *s_buffer, size_t n, struct dhist *histogram, const char *title,
    const char *xlabel, int mp_flag);
void dhist_set_num_buckets(
    struct dhist *histogram, uint32_t target_buckets, int mp_flag);
uint32_t dhist_get_num_buckets(
    struct dhist *histogram, bool get_actual_instead_of_target);
// This is only safe when called from a single-threaded context.
void dhist_set_decay_rate(struct dhist *histogram, double decay_rate);
double dhist_get_decay_rate(struct dhist *histogram);

#ifdef __cplusplus
}
#endif

#endif  // DECAYING_HISTOGRAM_SRC_DHIST_H_

