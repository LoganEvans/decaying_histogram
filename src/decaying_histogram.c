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

#include "decaying_histogram.h"

// The same as ceil(x / y). Using this so that math.h is not a dependency.
#define CEIL(x, y)                                                            \
    fixedpt_toint(fixedpt_xdiv((x) + ((y) - fixedpt_rconst(1)), (y)))


static bool is_in_bucket(struct bucket *bucket, fixedpt value);
static void recompute_bounds(struct bucket *bucket);
static void decay(
    struct decaying_histogram *histogram, struct bucket *bucket,
    uint64_t generation);
static void delete_bucket(
    struct decaying_histogram *histogram, int bucket_idx);
static void split_bucket(struct decaying_histogram *histogram, int bucket_idx);
static fixedpt get_decay(
    struct decaying_histogram *histogram, int missed_generations);
static void assert_consistent(struct decaying_histogram *histogram);

void
fixedpt_print(fixedpt A)
{
  char num[100];

  fixedpt_str(A, num, -2);
  printf("%s", num);
}

void init_bucket(
    struct bucket *to_init, struct bucket *below, struct bucket *above,
    fixedpt alpha) {
  to_init->count = fixedpt_rconst(0);
  to_init->mu = fixedpt_rconst(0);
  to_init->lower_bound = fixedpt_rconst(0);
  to_init->upper_bound = fixedpt_rconst(0);
  to_init->last_decay_generation = 0;
  to_init->below = below;
  to_init->above = above;
  pthread_mutexattr_init(&to_init->mutexattr);
  pthread_mutexattr_settype(&to_init->mutexattr, PTHREAD_MUTEX_RECURSIVE);
  pthread_mutex_init(&to_init->mutex, &to_init->mutexattr);
}

bool is_in_bucket(struct bucket *bucket, fixedpt observation) {
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
    if (bucket->mu == bucket->below->mu) {
      bucket->lower_bound = bucket->mu;
    } else {
      bucket->lower_bound = fixedpt_xdiv(
          (fixedpt_xmul(bucket->mu, bucket->count) +
           fixedpt_xmul(bucket->below->mu, bucket->below->count)),
          (bucket->count + bucket->below->count));
    }
    bucket->below->upper_bound = bucket->lower_bound;
  }

  if (bucket->above != NULL) {
    if (bucket->mu == bucket->above->mu) {
      bucket->upper_bound = bucket->mu;
    } else {
      bucket->upper_bound = fixedpt_xdiv(
          (fixedpt_xmul(bucket->mu, bucket->count) +
           fixedpt_xmul(bucket->above->mu, bucket->above->count)),
          (bucket->count + bucket->above->count));
    }
    bucket->above->lower_bound = bucket->upper_bound;
  }

  if (bucket->below == NULL && bucket->above == NULL) {
    // Only one bin exists. This is a bit of a fiction, but in order to
    // make the histogram's area sum to 1.0, we'll just center a square
    // at mu. This also depends on the density function returning 1.0 when
    // only one bucket exists.
    bucket->lower_bound = bucket->mu - fixedpt_rconst(0.5);
    bucket->upper_bound = bucket->mu + fixedpt_rconst(0.5);
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

  bucket->count = fixedpt_xmul(
      bucket->count,
      get_decay(histogram, generation - bucket->last_decay_generation));
  bucket->last_decay_generation = generation;
}

fixedpt density(
    struct decaying_histogram *histogram, struct bucket *bucket,
    fixedpt *lower_bound_output, fixedpt *upper_bound_output) {
  uint64_t generation;
  fixedpt retval;

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

  if (bucket->lower_bound == bucket->upper_bound) {
    // This is nonsensical, but we'll assume the bucket has width 1.0 (see
    // recompute_bounds()) so that the total area is 1.0.
    retval = fixedpt_rconst(1.0);
  } else {
    // The total count should approach (1.0 / alpha),
    // but in the warmup phase, we won't have that many observations
    // recorded.
    retval = fixedpt_xdiv(
        fixedpt_xdiv(bucket->count, total_count(histogram)),
        bucket->upper_bound - bucket->lower_bound);
  }

  if (bucket->below)
    pthread_mutex_unlock(&bucket->below->mutex);
  pthread_mutex_unlock(&bucket->mutex);
  if (bucket->above)
    pthread_mutex_unlock(&bucket->above->mutex);

  return retval;
}

fixedpt get_decay(
    struct decaying_histogram *histogram, int missed_generations) {
  return fixedpt_pow(
      fixedpt_rconst(1.0) - histogram->alpha,
      fixedpt_rconst(missed_generations));
}

void init_decaying_histogram(
    struct decaying_histogram *histogram, int target_buckets, double alpha) {
  fixedpt max_total_count, expected_count;
  int idx;
  fixedpt fixedpt_alpha;

  fixedpt_alpha = fixedpt_rconst(alpha);
  max_total_count = fixedpt_xdiv(fixedpt_rconst(1.0), fixedpt_alpha);
  expected_count =
      fixedpt_xdiv(
        fixedpt_rconst(1.0),
        fixedpt_xmul(fixedpt_alpha, fixedpt_fromint(target_buckets)));

  // TODO: What are good thresholds?
  histogram->delete_bucket_threshold =
      fixedpt_xdiv(expected_count, fixedpt_rconst(2.0));
  histogram->split_bucket_threshold =
      fixedpt_xdiv(
        fixedpt_xmul(fixedpt_rconst(3.0), expected_count),
        fixedpt_rconst(2.0));
  histogram->alpha = fixedpt_alpha;
  histogram->max_num_buckets =
      CEIL(max_total_count, histogram->delete_bucket_threshold);
  histogram->bucket_list = (struct bucket *)malloc(
      sizeof(struct bucket) * histogram->max_num_buckets);
  histogram->num_buckets = 1;
  for (idx = 0; idx < histogram->max_num_buckets; idx++)
    init_bucket(&histogram->bucket_list[idx], NULL, NULL, fixedpt_alpha);
  histogram->generation = 0;

  pthread_rwlock_init(&histogram->rwlock, NULL);
  pthread_mutex_init(&histogram->generation_mutex, NULL);

  return;
}

fixedpt total_count(struct decaying_histogram *histogram) {
  return fixedpt_xdiv(
      (fixedpt_rconst(1.0) - get_decay(histogram, histogram->generation)),
      histogram->alpha);
}

void clean_decaying_histogram(struct decaying_histogram *histogram) {
  free(histogram->bucket_list);
}

int find_bucket_idx(
    struct decaying_histogram *histogram, fixedpt observation) {
  int low, mid, high;
  struct bucket *bucket;

  low = 0;
  high = histogram->num_buckets - 1;
  while (low <= high) {
    mid = (low + high) / 2;
    bucket = &histogram->bucket_list[mid];
    if (is_in_bucket(bucket, observation))
      return mid;
    else if (bucket->mu < observation)
      low = mid + 1;
    else
      high = mid;
  }
  assert(false);
  return low;
}

void full_refresh(struct decaying_histogram *histogram) {
  int idx;
  bool do_another_round, finished_deletes;
  struct bucket *bucket;

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
}

void add_observation(
    struct decaying_histogram *histogram, double observation) {
  int bucket_idx;
  struct bucket *bucket;
  uint64_t generation;
  fixedpt fixedpt_observation;

  fixedpt_observation = fixedpt_rconst(observation);

  pthread_rwlock_rdlock(&histogram->rwlock);

  bucket_idx = find_bucket_idx(histogram, fixedpt_observation);
  bucket = &histogram->bucket_list[bucket_idx];

  // We'll lock everything from left to right to avoid deadlocks.
  if (bucket->below)
    pthread_mutex_lock(&bucket->below->mutex);
  pthread_mutex_lock(&bucket->mutex);
  if (bucket->above)
    pthread_mutex_lock(&bucket->above->mutex);

  while (!is_in_bucket(bucket, fixedpt_observation)) {
    // It is possible that the bucket boundaries could shift between when
    // the bucket is identified and when we lock it. If we fail, try again.
    if (bucket->above)
      pthread_mutex_unlock(&bucket->above->mutex);
    pthread_mutex_unlock(&bucket->mutex);
    if (bucket->below)
      pthread_mutex_unlock(&bucket->below->mutex);

    bucket_idx = find_bucket_idx(histogram, fixedpt_observation);
    bucket = &histogram->bucket_list[bucket_idx];

    if (bucket->below)
      pthread_mutex_lock(&bucket->below->mutex);
    pthread_mutex_lock(&bucket->mutex);
    if (bucket->above)
      pthread_mutex_lock(&bucket->above->mutex);
  }

  pthread_mutex_lock(&histogram->generation_mutex);
  ++histogram->generation;
  generation = histogram->generation;
  pthread_mutex_unlock(&histogram->generation_mutex);

  decay(histogram, bucket->below, generation);
  decay(histogram, bucket, generation);
  decay(histogram, bucket->above, generation);

  bucket->mu = fixedpt_xdiv(
      fixedpt_xmul(bucket->count, bucket->mu) + fixedpt_observation,
      (bucket->count + fixedpt_rconst(1.0)));
  bucket->count = bucket->count + fixedpt_rconst(1.0);
  recompute_bounds(bucket);

  if (bucket->above)
    pthread_mutex_unlock(&bucket->above->mutex);
  pthread_mutex_unlock(&bucket->mutex);
  if (bucket->below)
    pthread_mutex_unlock(&bucket->below->mutex);

  pthread_rwlock_unlock(&histogram->rwlock);

  if (bucket->count < histogram->delete_bucket_threshold ||
      bucket->count > histogram->split_bucket_threshold) {
    pthread_rwlock_wrlock(&histogram->rwlock);
    full_refresh(histogram);
    pthread_rwlock_unlock(&histogram->rwlock);
  }
}

void print_histogram(struct decaying_histogram *histogram) {
  int idx;
  fixedpt bound;

  pthread_rwlock_wrlock(&histogram->rwlock);
  full_refresh(histogram);

  printf("{\"densities\": [");
  for (idx = 0; idx < histogram->num_buckets - 1; idx++) {
    fixedpt_print(
        density(histogram, &histogram->bucket_list[idx], NULL, NULL));
    printf(", ");
  }
  fixedpt_print(
      density(histogram, &histogram->bucket_list[histogram->num_buckets - 1],
      NULL, NULL));
  printf("],\n\"boundaries\": [");
  for (idx = 0; idx < histogram->num_buckets; idx++) {
    density(histogram, &histogram->bucket_list[idx], &bound, NULL);
    fixedpt_print(bound);
    printf(", ");
  }
  density(
      histogram, &histogram->bucket_list[histogram->num_buckets - 1],
      NULL, &bound);
  fixedpt_print(bound);
  printf("]}\n");
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
      fixedpt_xdiv(
        (fixedpt_xmul(lucky_bucket->mu, lucky_bucket->count) +
         fixedpt_xmul(dying_bucket->mu, dying_bucket->count)),
        (lucky_bucket->count + dying_bucket->count));
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
  fixedpt diameter, median;
  struct bucket *left, *right, *far_right;

  if (histogram->num_buckets == 2 &&
      histogram->bucket_list[0].mu == histogram->bucket_list[1].mu) {
    // This will happen if we observe a stream of constants. Let's avoid
    // making more than two buckets until we have observed two unique values.
    return;
  }

  // Shift everything right.
  // XXX Does this shift one more than it needs to?
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
  histogram->bucket_list[histogram->num_buckets - 1].above =
      &histogram->bucket_list[histogram->num_buckets];
  histogram->bucket_list[histogram->num_buckets].below =
      &histogram->bucket_list[histogram->num_buckets - 1];
  ++histogram->num_buckets;

  // left is the bucket we are splitting.
  left = &histogram->bucket_list[bucket_idx];
  right = &histogram->bucket_list[bucket_idx + 1];
  if (bucket_idx + 2 < histogram->num_buckets)
    far_right = &histogram->bucket_list[bucket_idx + 2];
  else
    far_right = NULL;

  init_bucket(right, left, far_right, histogram->alpha);
  left->above = right;
  left->count = fixedpt_xdiv(left->count, fixedpt_rconst(2.0));
  right->count = left->count;

  if (histogram->num_buckets == 2) {
    // This is an awkward split because we don't have enough information to
    // compute a diameter. Instead, we can set mu to the same value in both
    // buckets, which will also make the boundary between them be mu. As soon
    // as we observe a non-mu entry, the split between the buckets
    // will be more sensical. However, we will need to take care to not split
    // a bucket again until this happens.
    right->mu = left->mu;
  } else {
    diameter = left->upper_bound - left->lower_bound;
    median = left->lower_bound + fixedpt_xdiv(diameter, fixedpt_rconst(2.0));
    left->mu = median - fixedpt_xdiv(diameter, fixedpt_rconst(6.0));
    right->mu = median + fixedpt_xdiv(diameter, fixedpt_rconst(6.0));
  }

  right->last_decay_generation = left->last_decay_generation;

  recompute_bounds(left);
  recompute_bounds(right);

  return;
}

fixedpt Jaccard_distance(
    struct decaying_histogram *hist0, struct decaying_histogram *hist1) {
  assert(false);
  return 0.0;
}

fixedpt Kolmogorov_Smirnov_statistic(
    struct decaying_histogram *hist0, struct decaying_histogram *hist1) {
  assert(false);
  return 0.0;
}

void assert_consistent(struct decaying_histogram *histogram) {
  int idx;

  assert(histogram->num_buckets <= histogram->max_num_buckets);
  assert(histogram->num_buckets > 0);

  for (idx = 0; idx < histogram->num_buckets; idx++)
    decay(histogram, &histogram->bucket_list[idx], histogram->generation);

  for (idx = 1; idx < histogram->num_buckets; idx++) {
    if (histogram->bucket_list[idx - 1].mu > histogram->bucket_list[idx].mu) {
      assert(false);
    }
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
      printf("(%d / %d) bucket_list[%d].upper_bound (==%s) != "
             "bucket_list[%d].lower_bound (==%s)\n",
          idx, histogram->num_buckets,
          idx,
          fixedpt_cstr(histogram->bucket_list[idx].upper_bound, -2),
          idx + 1,
          fixedpt_cstr(histogram->bucket_list[idx + 1].lower_bound, -2));
      printf("!!!!!\n");
      for (idx = 0; idx < histogram->num_buckets; idx++) {
        printf("??? %s | %s\n",
            fixedpt_cstr(histogram->bucket_list[idx].lower_bound, -2),
            fixedpt_cstr(histogram->bucket_list[idx].upper_bound, -2));
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

