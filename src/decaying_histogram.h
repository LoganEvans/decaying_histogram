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

struct bucket {
    double alpha_mu;
    double alpha_count;
    double count;
    double mu;
    double min;
    double max;
    double lower_bound;
    double upper_bound;
    uint64_t last_decay_generation;
    struct bucket *below;
    struct bucket *above;
};

void init_bucket(
        struct bucket *to_init, struct bucket *below, struct bucket *above);
void add_observation_to_bucket(struct bucket *bucket, double observation);
bool is_in_bucket(struct bucket *bucket, double value);
void recompute_bound(
        struct bucket *lower, struct bucket *upper, uint64_t generation);
void decay(struct bucket *bucket, uint64_t generation);
double density(struct bucket *bucket, uint64_t generation, double total);
double height(struct bucket *bucket, uint64_t height, double total);


struct decaying_histogram {
    double total_count;
    double delete_bucket_threshold;
    double split_bucket_threshold;
    uint64_t generation;
    struct bucket *bucket_list;
    int num_buckets;
    int max_num_buckets;
};

void init_decaying_histogram(
        struct decaying_histogram *histogram, int target_buckets,
        double alpha_mu, double alpha_count);
void destroy_buckets(struct decaying_histogram *histogram);
struct bucket * find_bucket(
        struct decaying_histogram *histogram, double observation);
void full_refresh(struct decaying_histogram *histogram);
void add_observation(struct decaying_histogram *histogram, double observation);
void sprint_histogram_new(char **str);
void num_buckets(struct decaying_histogram *histogram);
void get_CDF(
        struct decaying_histogram *histogram, double *cdf, int *num_buckets);
void delete_bucket(
        struct decaying_histogram *histogram, struct bucket *bucket);
void split_bucket(struct decaying_histogram *histogram, struct bucket *bucket);
double Jaccard_distance(
        struct decaying_histogram *hist0, struct decaying_histogram *hist1);
double Kolomogorov_Smirnov_statistic(
        struct decaying_histogram *hist0, struct decaying_histogram *hist1);

#ifdef __cplusplus
}
#endif

#endif  /* DECAYING_HISTOGRAM_H_ */

