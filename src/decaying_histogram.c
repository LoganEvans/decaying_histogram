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
#include <stdlib.h>
#include <stdio.h>
// XXX
#include <math.h>

#include "decaying_histogram.h"

// The same as ceil(x / y). Using this so that math.h is not a dependency.
#define CEIL(x, y) ((int)((x) + ((int)(y) - 1)) / (int)(y))


static bool is_in_bucket(struct bucket *bucket, double value);
static void recompute_bounds(struct bucket *bucket);
static void decay(
    struct decaying_histogram *histogram, struct bucket *bucket,
    uint64_t generation);
static void delete_bucket(
    struct decaying_histogram *histogram, int bucket_idx);
static void split_bucket(struct decaying_histogram *histogram, int bucket_idx);
static double ipow(double coefficient, uint64_t power);
static double get_decay(
    struct decaying_histogram *histogram, uint64_t missed_generations);
//static void assert_consistent(struct decaying_histogram *histogram);


void init_bucket(
    struct bucket *to_init, struct bucket *below, struct bucket *above) {
  to_init->count = 0.0;
  to_init->mu = 0.0;
  to_init->lower_bound = 0.0;
  to_init->upper_bound = 0.0;
  to_init->last_decay_generation = 0;
  to_init->below = below;
  to_init->above = above;
  pthread_mutexattr_init(&to_init->mutexattr);
  pthread_mutexattr_settype(&to_init->mutexattr, PTHREAD_MUTEX_RECURSIVE);
  pthread_mutex_init(&to_init->mutex, &to_init->mutexattr);
}

bool is_in_bucket(struct bucket *bucket, double observation) {
  if ((bucket->below == NULL || bucket->lower_bound <= observation) &&
      (bucket->above == NULL || observation < bucket->upper_bound)) {
    return true;
  } else {
    return false;
  }
}

/* Both buckets should be decayed before this is called. */
void recompute_bounds(struct bucket *bucket) {
  if (bucket == NULL)
    return;

  // TODO: This might need to check against DBL_MAX.

  // Note: the mutexes for bucket->below, bucket, and bucket->above
  // should be held (or else we have the wrlock);

  if (bucket->below != NULL) {
    bucket->lower_bound = (bucket->mu * bucket->count +
                           bucket->below->mu * bucket->below->count) /
                          (bucket->count + bucket->below->count);
    bucket->below->upper_bound = bucket->lower_bound;
  }

  if (bucket->above != NULL) {
    bucket->upper_bound = (bucket->mu * bucket->count +
                           bucket->above->mu * bucket->above->count) /
                          (bucket->count + bucket->above->count);
    bucket->above->lower_bound = bucket->upper_bound;
  }

  if (bucket->below == NULL && bucket->above == NULL) {
    // Only one bin exists. This is a bit of a fiction, but in order to
    // make the histogram's area sum to 1.0, we'll just center a square
    // at mu. This also depends on the density function returning 1.0 when
    // only one bucket exists.
    bucket->lower_bound = bucket->mu - 0.5;
    bucket->upper_bound = bucket->mu + 0.5;
  } else if (bucket->below == NULL) {
    // While the bound is effectively infinate, we need a sensible bound to
    // compute bucket heights. Here, we'll assume that the observations in
    // the bucket are uniform, so mu should be at the very center of the
    // bucket.
    bucket->lower_bound = bucket->mu - (bucket->upper_bound - bucket->mu);
  } else if (bucket->above == NULL) {
    bucket->upper_bound = bucket->mu + (bucket->mu - bucket->lower_bound);
  }
}

/*
 * This can efficiently apply the decay for multiple generations. This function
 * will not read histogram->generation since that may be updated by another
 * thread.
 */
void decay(
    struct decaying_histogram *histogram, struct bucket *bucket,
    uint64_t generation) {
  if (bucket == NULL || bucket->last_decay_generation == generation)
    return;

  bucket->count *= get_decay(
      histogram, generation - bucket->last_decay_generation);
  bucket->last_decay_generation = generation;
}

double density(
    struct decaying_histogram *histogram, struct bucket *bucket,
    double *lower_bound_output, double *upper_bound_output) {
  uint64_t generation;
  double retval;

  if (bucket->below)
    pthread_mutex_lock(&bucket->below->mutex);
  pthread_mutex_lock(&bucket->mutex);
  if (bucket->above)
    pthread_mutex_lock(&bucket->above->mutex);

  pthread_mutex_lock(&histogram->generation_mutex);
  generation = histogram->generation;
  pthread_mutex_unlock(&histogram->generation_mutex);

  decay(histogram, bucket->below, generation);
  decay(histogram, bucket, generation);
  decay(histogram, bucket->above, generation);
  recompute_bounds(bucket);
  if (lower_bound_output != NULL)
    *lower_bound_output = bucket->lower_bound;
  if (upper_bound_output != NULL)
    *upper_bound_output = bucket->upper_bound;

  if (histogram->num_buckets == 1 ||
      histogram->generation == 0) {
    // This is nonsensical, but we'll assume the bucket has width 1.0 (see
    // recompute_bounds()) so that the total area is 1.0.
    retval = 1.0;
  } else {
    // The total count should approach (1.0 / alpha),
    // but in the warmup phase, we won't have that many observations
    // recorded.
    //printf("??? %lf %lf %lf %lf %d\n",
    //  bucket->lower_bound, bucket->upper_bound,
    //  bucket->count, total_count(histogram),
    //  histogram->num_buckets);
    retval = (bucket->count / total_count(histogram)) /
             (bucket->upper_bound - bucket->lower_bound);
  }

  if (bucket->below)
    pthread_mutex_unlock(&bucket->below->mutex);
  pthread_mutex_unlock(&bucket->mutex);
  if (bucket->above)
    pthread_mutex_unlock(&bucket->above->mutex);

  return retval;
}

double ipow(double coefficient, uint64_t power) {
  double result;

  result = 1.0;
  while (power) {
    if (power & 1)
      result *= coefficient;
    power >>= 1;
    coefficient *= coefficient;
  }
  return result;
}

double get_decay(
    struct decaying_histogram *histogram, uint64_t missed_generations) {
  if (missed_generations < histogram->max_num_buckets)
    return histogram->pow_table[missed_generations];
  else
    return ipow(1.0 - histogram->alpha, missed_generations);
}

void init_decaying_histogram(
    struct decaying_histogram *histogram, int target_buckets, double alpha) {
  double max_total_count, expected_count;
  uint64_t idx;

  max_total_count = 1.0 / alpha;
  expected_count = 1.0 / (alpha * target_buckets);

  // TODO: What are good thresholds?
  histogram->delete_bucket_threshold = expected_count / 2.0;
  histogram->split_bucket_threshold = (3.0 * expected_count) / 2.0;
  histogram->alpha = alpha;
  histogram->max_num_buckets =
      (unsigned int)CEIL(max_total_count, histogram->delete_bucket_threshold);
  histogram->bucket_list = (struct bucket *)malloc(
      sizeof(struct bucket) * histogram->max_num_buckets);
  histogram->num_buckets = 1;
  for (idx = 0; idx < histogram->max_num_buckets; idx++)
    init_bucket(&histogram->bucket_list[idx], NULL, NULL);
  histogram->generation = 0;

  histogram->pow_table = (double *)malloc(
      histogram->max_num_buckets * sizeof(double));
  for (idx = 0; idx < histogram->max_num_buckets; idx++)
    histogram->pow_table[idx] = ipow(1.0 - alpha, idx);

  pthread_rwlock_init(&histogram->rwlock, NULL);
  pthread_mutex_init(&histogram->generation_mutex, NULL);

  return;
}

double total_count(struct decaying_histogram *histogram) {
  return (1 - get_decay(histogram, histogram->generation)) /
         histogram->alpha;
}

void clean_decaying_histogram(struct decaying_histogram *histogram) {
  free(histogram->bucket_list);
  free(histogram->pow_table);
}

int find_bucket_idx(
    struct decaying_histogram *histogram, double observation) {
  int low, mid, high;
  struct bucket *bucket;

  low = 0;
  high = histogram->num_buckets - 1;
  while (low < high) {
    mid = (low + high) / 2;
    bucket = &histogram->bucket_list[mid];
    if (is_in_bucket(bucket, observation))
      return mid;
    else if (bucket->mu < observation)
      low = mid + 1;
    else
      high = mid;
  }
  return low;
}

void full_refresh(struct decaying_histogram *histogram) {
  int idx;
  bool do_another_round, finished_deletes;
  struct bucket *bucket;

  pthread_rwlock_wrlock(&histogram->rwlock);
  for (idx = 0; idx < histogram->num_buckets; idx++)
    decay(histogram, &histogram->bucket_list[idx], histogram->generation);

  /*
   * In order to make sure we don't overrun the number of buckets available,
   * we need to condense before we expand.
   */
  finished_deletes = false;
  do {
    do_another_round = false;
    for (idx = 0; idx < histogram->num_buckets; idx++) {
      bucket = &histogram->bucket_list[idx];
      if (!finished_deletes &&
          bucket->count < histogram->delete_bucket_threshold) {
        delete_bucket(histogram, idx);
        do_another_round = true;
        break; /* delete_bucket reorders the buckets */
      } else if (finished_deletes &&
          bucket->count > histogram->split_bucket_threshold) {
        split_bucket(histogram, idx);
        do_another_round = true;
        break;
      }
    }
    if (!finished_deletes) {
      finished_deletes = true;
      do_another_round = true;
    }
  } while (do_another_round);

  for (idx = 0; idx < histogram->num_buckets; idx++)
    decay(histogram, &histogram->bucket_list[idx], histogram->generation);

  for (idx = 0; idx < histogram->num_buckets; idx++)
    recompute_bounds(&histogram->bucket_list[idx]);

  pthread_rwlock_unlock(&histogram->rwlock);
}

void add_observation(
    struct decaying_histogram *histogram, double observation) {
  int bucket_idx;
  struct bucket *bucket;
  uint64_t generation;

  pthread_rwlock_rdlock(&histogram->rwlock);

  for (;;) {
    bucket_idx = find_bucket_idx(histogram, observation);
    bucket = &histogram->bucket_list[bucket_idx];

    // We'll lock everything from left to right to avoid deadlocks.
    if (bucket->below)
      pthread_mutex_lock(&bucket->below->mutex);
    pthread_mutex_lock(&bucket->mutex);
    if (bucket->above)
      pthread_mutex_lock(&bucket->above->mutex);

    if (is_in_bucket(bucket, observation))
      break;

    // It is possible that the bucket boundaries could shift between when
    // the bucket is identified and when we lock it. If we fail, try again.
    if (bucket->above)
      pthread_mutex_unlock(&bucket->above->mutex);
    pthread_mutex_unlock(&bucket->mutex);
    if (bucket->below)
      pthread_mutex_unlock(&bucket->below->mutex);
  }

  pthread_mutex_lock(&histogram->generation_mutex);
  ++histogram->generation;
  generation = histogram->generation;
  pthread_mutex_unlock(&histogram->generation_mutex);

  decay(histogram, bucket->below, generation);
  decay(histogram, bucket, generation);
  decay(histogram, bucket->above, generation);
  bucket->mu =
      (bucket->count * bucket->mu + observation) / (bucket->count + 1);
  ++bucket->count;

  recompute_bounds(bucket);

  if (bucket->above)
    pthread_mutex_unlock(&bucket->above->mutex);
  pthread_mutex_unlock(&bucket->mutex);
  if (bucket->below)
    pthread_mutex_unlock(&bucket->below->mutex);

  pthread_rwlock_unlock(&histogram->rwlock);

  if (bucket->count < histogram->delete_bucket_threshold)
    full_refresh(histogram);
  if (bucket->count > histogram->split_bucket_threshold)
    full_refresh(histogram);
}

void print_histogram(struct decaying_histogram *histogram) {
  int idx;
  double bound;

  full_refresh(histogram);
  pthread_rwlock_rdlock(&histogram->rwlock);
  printf("{\"densities\": [");
  for (idx = 0; idx < histogram->num_buckets - 1; idx++) {
    printf(
        "%lf, ",
        density(histogram, &histogram->bucket_list[idx], NULL, NULL));
  }
  printf(
      "%lf], ",
      density(
        histogram, &histogram->bucket_list[histogram->num_buckets - 1],
        NULL, NULL));

  printf("\"boundaries\": [");
  for (idx = 0; idx < histogram->num_buckets; idx++) {
    density(histogram, &histogram->bucket_list[idx], &bound, NULL);
    printf("%lf, ", bound);
  }
  density(
      histogram, &histogram->bucket_list[histogram->num_buckets - 1],
      NULL, &bound);
  printf("%lf]}\n", bound);
  pthread_rwlock_unlock(&histogram->rwlock);
}

/*
 * All buckets must be decayed before calling this.
 */
void delete_bucket(
    struct decaying_histogram *histogram, int bucket_idx) {
  int lucky_idx, idx;
  struct bucket *lucky_bucket, *dying_bucket;

  if (histogram->num_buckets <= 2)
    return;  // Let's not.

  dying_bucket = &histogram->bucket_list[bucket_idx];
  if (bucket_idx == 0)
    lucky_idx = bucket_idx + 1;
  else if (bucket_idx == histogram->num_buckets - 1)
    lucky_idx = bucket_idx - 1;
  else if (dying_bucket->below->count < dying_bucket->above->count)
    lucky_idx = bucket_idx - 1;
  else
    lucky_idx = bucket_idx + 1;

  lucky_bucket = &histogram->bucket_list[lucky_idx];
  lucky_bucket->mu =
      (lucky_bucket->mu * lucky_bucket->count +
       dying_bucket->mu * dying_bucket->count) /
      (lucky_bucket->count + dying_bucket->count);
  lucky_bucket->count += dying_bucket->count;

  // Shift everything left.
  for (idx = bucket_idx; idx < histogram->num_buckets - 1; idx++) {
    histogram->bucket_list[idx].count = histogram->bucket_list[idx + 1].count;
    histogram->bucket_list[idx].mu = histogram->bucket_list[idx + 1].mu;
    histogram->bucket_list[idx].lower_bound =
        histogram->bucket_list[idx + 1].lower_bound;
    histogram->bucket_list[idx].upper_bound =
        histogram->bucket_list[idx + 1].upper_bound;
    histogram->bucket_list[idx].last_decay_generation =
        histogram->bucket_list[idx + 1].last_decay_generation;
  }
  --histogram->num_buckets;
  histogram->bucket_list[histogram->num_buckets].below = NULL;
  histogram->bucket_list[histogram->num_buckets - 1].above = NULL;

  recompute_bounds(&histogram->bucket_list[bucket_idx]);

  return;
}

void split_bucket(
    struct decaying_histogram *histogram, int bucket_idx) {
  int idx;
  double diameter, median;
  struct bucket *left, *right, *far_right;

  if (histogram->num_buckets == 2 &&
      histogram->bucket_list[0].mu == histogram->bucket_list[1].mu) {
    // This will happen if we observe a stream of constants. Let's avoid
    // making more than two buckets until we have observed two unique values.
    return;
  }

  // Shift everything right.
  for (idx = histogram->num_buckets; idx > bucket_idx; idx--) {
    histogram->bucket_list[idx].count = histogram->bucket_list[idx - 1].count;
    histogram->bucket_list[idx].mu = histogram->bucket_list[idx - 1].mu;
    histogram->bucket_list[idx].lower_bound =
        histogram->bucket_list[idx - 1].lower_bound;
    histogram->bucket_list[idx].upper_bound =
        histogram->bucket_list[idx - 1].upper_bound;
    histogram->bucket_list[idx].last_decay_generation =
        histogram->bucket_list[idx - 1].last_decay_generation;
  }
  histogram->bucket_list[histogram->num_buckets].above = NULL;
  if (histogram->num_buckets > 0) {
    histogram->bucket_list[histogram->num_buckets - 1].above =
        &histogram->bucket_list[histogram->num_buckets];
    histogram->bucket_list[histogram->num_buckets].below =
        &histogram->bucket_list[histogram->num_buckets - 1];
  } else {
    histogram->bucket_list[histogram->num_buckets].below = NULL;
  }
  ++histogram->num_buckets;

  // left is the bucket we are splitting.
  left = &histogram->bucket_list[bucket_idx];
  right = &histogram->bucket_list[bucket_idx + 1];
  if (bucket_idx + 2 < histogram->num_buckets)
    far_right = &histogram->bucket_list[bucket_idx + 2];
  else
    far_right = NULL;

  init_bucket(right, left, far_right);
  left->above = right;
  left->count /= 2.0;
  right->count = left->count;
  if (histogram->num_buckets == 2) {
    // This is an awkward split because we don't have enough information to
    // compute a diameter. Instead, we can set mu to the same value in both
    // buckets, which will also make the boundary between them be mu. As soon
    // as we observe a non-mu entry, the split between the buckets
    // will be more sensical. However, we will need to take care to not split
    // a bucket again until this happens.
    left->count /= 2.0;
    right->count = left->count;
    right->mu = left->mu;
  } else {
    diameter = left->upper_bound - left->lower_bound;
    median = left->lower_bound + diameter / 2.0;
    left->mu = median - diameter / 6.0;
    right->mu = median + diameter / 6.0;
  }

  right->last_decay_generation = left->last_decay_generation;

  recompute_bounds(left);
  recompute_bounds(right);

  return;
}

double Jaccard_distance(
    struct decaying_histogram *hist0, struct decaying_histogram *hist1) {
  assert(hist0 && hist1);
  assert(false);
  return 0.0;
}

double Kolomogorov_Smirnov_statistic(
    struct decaying_histogram *hist0, struct decaying_histogram *hist1) {
  assert(hist0 && hist1);
  assert(false);
  return 0.0;
}

#if 0
void assert_consistent(struct decaying_histogram *histogram) {
  int idx;
  struct bucket *bucket;

  assert(histogram->num_buckets <= histogram->max_num_buckets);
  assert(histogram->num_buckets > 0);

  for (idx = 0; idx < histogram->num_buckets; idx++)
    decay(histogram, &histogram->bucket_list[idx], histogram->generation);
  for (idx = 0; idx < histogram->num_buckets; idx++) {
    bucket = &histogram->bucket_list[idx];
    recompute_bounds(bucket);
    assert(
        bucket->lower_bound <= bucket->mu &&
        bucket->mu <= bucket->upper_bound);
    assert(bucket->mu != NAN);
    assert(bucket->mu != INFINITY);
    assert(bucket->lower_bound != NAN);
    assert(bucket->lower_bound != INFINITY);
    assert(bucket->upper_bound != NAN);
    assert(bucket->upper_bound != INFINITY);
  }

  for (idx = 1; idx < histogram->num_buckets - 1; idx++) {
    assert(
        &histogram->bucket_list[idx - 1] == histogram->bucket_list[idx].below);
    if (&histogram->bucket_list[idx + 1] !=
        histogram->bucket_list[idx].above) {
      printf("(%d / %d) &bucket_list[%d] (==%lx) !="
             "bucket_list[%d].above (==%lx)\n",
          idx, histogram->num_buckets,
          idx + 1, (uint64_t)&histogram->bucket_list[idx + 1],
          idx, (uint64_t)histogram->bucket_list[idx].above);
      assert(false);
    }
    if (histogram->bucket_list[idx].upper_bound !=
        histogram->bucket_list[idx + 1].lower_bound) {
      printf("(%d / %d) bucket_list[%d].upper_bound (==%lf) != "
             "bucket_list[%d].lower_bound (==%lf)\n",
          idx, histogram->num_buckets,
          idx, histogram->bucket_list[idx].upper_bound,
          idx + 1, histogram->bucket_list[idx + 1].lower_bound);
      printf("!!!!!\n");
      for (idx = 0; idx < histogram->num_buckets; idx++) {
        printf("??? %lf | %lf\n",
            histogram->bucket_list[idx].lower_bound,
            histogram->bucket_list[idx].upper_bound);
      }
      assert(false);
    }
  }
  assert(NULL == histogram->bucket_list[0].below);
  assert(NULL == histogram->bucket_list[histogram->num_buckets].above);
  if (histogram->num_buckets > 1) {
    assert(
        histogram->bucket_list[0].upper_bound ==
        histogram->bucket_list[1].lower_bound);
    assert(
        histogram->bucket_list[histogram->num_buckets - 2].upper_bound ==
        histogram->bucket_list[histogram->num_buckets - 1].lower_bound);
  }
}
#endif

