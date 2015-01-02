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
#include <stdlib.h>

#include "decaying_histogram.h"


void init_bucket(
        struct bucket *to_init, struct bucket *below, struct bucket *above) {
    to_init->alpha_mu = 0.0;
    to_init->alpha_count = 0.0;
    to_init->count = 0.0;
    to_init->mu = 0.0;
    to_init->min = 0.0;
    to_init->max = 0.0;
    to_init->lower_bound = 0.0;
    to_init->upper_bound = 0.0;
    to_init->last_decay_generation = 0;
    to_init->below = below;
    to_init->above = above;
}

void add_observation_to_bucket(struct bucket *bucket, double observation) {
    assert(false);
    return;
}

bool is_in_bucket(struct bucket *bucket, double observation) {
    if ((bucket->below == NULL || bucket->lower_bound <= observation)
        && (bucket->above == NULL || observation < bucket->upper_bound))
        return true;
    else
        return false;
}

void recompute_bound(
        struct bucket *lower, struct bucket *upper, uint64_t generation) {
    double bound;

    if (lower)
        decay(lower, generation);
    if (upper)
        decay(upper, generation);

    if (lower == NULL || upper == NULL)
        return;  // The bound is effectively infinate.

    bound = ((lower->mu * lower->count + upper->mu * upper->count) /
             (lower->count + upper->count));
    lower->upper_bound = bound;
    upper->lower_bound = bound;
    return;
}

/*
 * This can efficiently apply the decay for multiple generations.
 */
void decay(struct bucket *bucket, uint64_t generation) {
    if (bucket->last_decay_generation == generation)
        return;
    bucket->count *= pow(
            1.0 - bucket->alpha_count,
            generation - bucket->last_decay_generation);
    bucket->last_decay_generation = generation;
    return;
}

double density(struct bucket *bucket, uint64_t generation, double total) {
    decay(bucket, generation);
    return bucket->alpha_count * bucket->count;  // ???
}

double height(struct bucket *bucket, uint64_t generation, double total) {
    recompute_bound(bucket->below, bucket, generation);
    recompute_bound(bucket, bucket->above, generation);
    return (density(bucket, generation, total) /
            ((bucket->above ? bucket->upper_bound : bucket->max) -
             (bucket->below ? bucket->lower_bound : bucket->min)));
}

void init_decaying_histogram(
        struct decaying_histogram *histogram, int target_buckets,
        double alpha_mu, double alpha_count) {
    double max_count = 1.0 / alpha_count;
    double expected_count = 1.0 / (alpha_count * target_buckets);

    histogram->delete_bucket_threshold = expected_count / 2.0;
    histogram->split_bucket_threshold = expected_count * 2.0;
    histogram->max_num_buckets =
            ceil(max_count / histogram->delete_bucket_threshold);
    histogram->bucket_list = (struct bucket *)malloc(
            sizeof(struct bucket) * histogram->max_num_buckets);
    histogram->num_buckets = 0;
    histogram->generation = 0;

    return;
}

void destroy_buckets(struct decaying_histogram *histogram) {
    free(histogram->bucket_list);
}

struct bucket * find_bucket(
        struct decaying_histogram *histogram, double observation) {
    assert(false);
    return NULL;
}

void full_refresh(struct decaying_histogram *histogram) {
    assert(false);
    return;
}

void add_observation(
        struct decaying_histogram *histogram, double observation) {
    assert(false);
    return;
}

void sprint_histogram_new(char **str) {
    assert(false);
    return;
}

void num_buckets(struct decaying_histogram *histogram) {
    assert(false);
    return;
}

void get_CDF(
        struct decaying_histogram *histogram, double *cdf, int *num_buckets) {
    assert(false);
    return;
}

void delete_bucket(
        struct decaying_histogram *histogram, struct bucket *bucket) {
    assert(false);
    return;
}

void split_bucket(
        struct decaying_histogram *histogram, struct bucket *bucket) {
    assert(false);
    return;
}

double Jaccard_distance(
        struct decaying_histogram *hist0, struct decaying_histogram *hist1) {
    assert(false);
    return 0.0;
}

double Kolomogorov_Smirnov_statistic(
        struct decaying_histogram *hist0, struct decaying_histogram *hist1) {
    assert(false);
    return 0.0;
}

