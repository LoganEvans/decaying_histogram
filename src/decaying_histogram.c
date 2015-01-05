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
#include <stdio.h>

#include "decaying_histogram.h"


void init_bucket(
    struct bucket *to_init, struct bucket *below, struct bucket *above,
    double alpha) {
  to_init->alpha = alpha;
  to_init->count = 0.0;
  to_init->mu = 0.0;
  to_init->lower_bound = 0.0;
  to_init->upper_bound = 0.0;
  to_init->last_decay_generation = 0;
  to_init->below = below;
  to_init->above = above;
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
void recompute_bound(struct bucket *lower, struct bucket *upper) {
  double bound;
  // TODO: This might need to check against DBL_MAX.

  if (lower == NULL && upper->above == NULL) {
    // Only one bin exists. Since we don't have any other way to compute
    // the bound, we'll insert the only number that we have.
    upper->lower_bound = upper->mu;
  } else if (lower == NULL) {
    // While the bound is effectively infinate, we need a sensible bound to
    // compute bucket heights. Here, we'll assume that the observations in
    // the bucket are uniform, so mu should be at the very center of the
    // bucket.
    // We can't assume that the bounds are accurate here, so we'll
    // recompute them.
    bound = (upper->mu * upper->count +
         upper->above->mu * upper->above->count) /
        (upper->count + upper->above->count);
    upper->lower_bound = upper->mu - (bound - upper->mu);
  } else if (upper == NULL && lower->below == NULL) {
    lower->upper_bound = lower->mu;
  } else if (upper == NULL) {
    bound = (lower->below->mu * lower->below->count +
         lower->mu * lower->count) /
        (lower->below->count + lower->count);
    lower->upper_bound = lower->mu + (lower->mu - bound);
  } else {
    // Both buckets exist.
    bound = ((lower->mu * lower->count + upper->mu * upper->count) /
         (lower->count + upper->count));
    lower->upper_bound = bound;
    upper->lower_bound = bound;
  }
}

/*
 * This can efficiently apply the decay for multiple generations.
 */
void decay(struct decaying_histogram *histogram, struct bucket *bucket) {
  if (bucket == NULL || bucket->last_decay_generation == histogram->generation)
    return;
  bucket->count *= get_decay(
      histogram, histogram->generation - bucket->last_decay_generation);
  //bucket->count *= pow(
  //    1.0 - bucket->alpha,
  //    histogram->generation - bucket->last_decay_generation);
  bucket->last_decay_generation = histogram->generation;
  return;
}

double density(struct decaying_histogram *histogram, struct bucket *bucket) {
  double left, right;

  decay(histogram, bucket->below);
  decay(histogram, bucket);
  decay(histogram, bucket->above);
  recompute_bound(bucket->below, bucket);
  recompute_bound(bucket, bucket->above);
  if (bucket->upper_bound == bucket->lower_bound || histogram->generation == 0)
    return 1.0;  // This is nonsensical, but it graphs better than MAX_DOUBLE.
  else
    // The total count should approach (1.0 / bucket->alpha),
    // but in the warmup phase, we won't have that many observations
    // recorded.
    return (bucket->count / total_count(histogram)) /
         (bucket->upper_bound - bucket->lower_bound);
}

double get_decay(
    struct decaying_histogram *histogram, int missed_generations) {
  double decay, base;
  int exp = missed_generations;

  if (missed_generations < histogram->max_num_buckets) {
    return histogram->pow_table[missed_generations];
  } else {
    decay = 1.0;
    base = 1.0 - histogram->alpha;
    while (missed_generations) {
      if (missed_generations & 1)
        decay *= base;
      missed_generations >>= 1;
      base *= base;
    }
    // XXX pow(1.0 - histogram->alpha, 10000000) == 1.0, but
    // this function returns 0.0.
    printf("??? %lf\n", decay);
    return decay;
  }
}

void init_decaying_histogram(
    struct decaying_histogram *histogram, int target_buckets, double alpha) {
  double max_total_count, expected_count;
  int idx;

  max_total_count = 1.0 / alpha;
  expected_count = 1.0 / (alpha * target_buckets);

  // TODO: What are good thresholds?
  histogram->delete_bucket_threshold = expected_count / 2.0;
  histogram->split_bucket_threshold = (3.0 * expected_count) / 2.0;
  histogram->alpha = alpha;
  histogram->max_num_buckets =
      ceil(max_total_count / histogram->delete_bucket_threshold);
  histogram->bucket_list = (struct bucket *)malloc(
      sizeof(struct bucket) * histogram->max_num_buckets);
  histogram->num_buckets = 1;
  for (idx = 0; idx < histogram->max_num_buckets; idx++)
    init_bucket(&histogram->bucket_list[idx], NULL, NULL, alpha);
  histogram->generation = 0;

  histogram->pow_table = (double *)malloc(
      histogram->max_num_buckets * sizeof(double));
  for (idx = 0; idx < histogram->max_num_buckets; idx++)
    histogram->pow_table[idx] = pow(alpha, idx);

  return;
}

double total_count(struct decaying_histogram *histogram) {
  return (1 - get_decay(histogram, histogram->generation)) / histogram->alpha;
}

void clean_decaying_histogram(struct decaying_histogram *histogram) {
  free(histogram->bucket_list);
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

  for (idx = 0; idx < histogram->num_buckets; idx++)
    decay(histogram, &histogram->bucket_list[idx]);

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

  for (idx = 0; idx < histogram->num_buckets - 1; idx++) {
    recompute_bound(
        &histogram->bucket_list[idx], &histogram->bucket_list[idx + 1]);
  }

  recompute_bound(NULL, &histogram->bucket_list[0]);
  recompute_bound(&histogram->bucket_list[histogram->num_buckets - 1], NULL);
}


void add_observation(
    struct decaying_histogram *histogram, double observation) {
  int bucket_idx;
  struct bucket *bucket;

  ++histogram->generation;
  bucket_idx = find_bucket_idx(histogram, observation);
  bucket = &histogram->bucket_list[bucket_idx];
  decay(histogram, bucket);
  bucket->mu =
      (bucket->count * bucket->mu + observation) / (bucket->count + 1);
  ++bucket->count;
  if (bucket_idx > 0) {
    decay(histogram, &histogram->bucket_list[bucket_idx - 1]);
    recompute_bound(
        &histogram->bucket_list[bucket_idx - 1],
        &histogram->bucket_list[bucket_idx]);
  }
  if (bucket_idx < histogram->num_buckets - 1) {
    decay(histogram, &histogram->bucket_list[bucket_idx + 1]);
    recompute_bound(
        &histogram->bucket_list[bucket_idx],
        &histogram->bucket_list[bucket_idx + 1]);
  }
  if (bucket->count < histogram->delete_bucket_threshold)
    full_refresh(histogram);
  if (bucket->count > histogram->split_bucket_threshold)
    full_refresh(histogram);
}

void print_histogram(struct decaying_histogram *histogram) {
  int idx;

  full_refresh(histogram);
  printf("[");
  for (idx = 0; idx < histogram->num_buckets - 1; idx++) {
    printf("%lf, ", density(histogram, &histogram->bucket_list[idx]));
  }
  printf(
      "%lf]\n",
      density(histogram, &histogram->bucket_list[histogram->num_buckets - 1]));

  printf("[");
  for (idx = 0; idx < histogram->num_buckets; idx++) {
    printf("%lf, ", histogram->bucket_list[idx].lower_bound);
  }
  printf(
      "%lf]\n",
      histogram->bucket_list[histogram->num_buckets - 1].upper_bound);
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
    //XXX
    //histogram->bucket_list[idx].below = histogram->bucket_list[idx + 1].below;
    //histogram->bucket_list[idx].above = histogram->bucket_list[idx + 1].above;
  }
  --histogram->num_buckets;
  histogram->bucket_list[histogram->num_buckets].below = NULL;
  histogram->bucket_list[histogram->num_buckets - 1].above = NULL;

  if (bucket_idx > 0) {
    recompute_bound(
        &histogram->bucket_list[bucket_idx - 1],
        &histogram->bucket_list[bucket_idx]);
  }
  if (bucket_idx + 1 < histogram->num_buckets) {
    recompute_bound(
        &histogram->bucket_list[bucket_idx],
        &histogram->bucket_list[bucket_idx + 1]);
  }

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
  assert(histogram->num_buckets < histogram->max_num_buckets);

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
    //XXX
    //histogram->bucket_list[idx].below = histogram->bucket_list[idx - 1].below;
    //histogram->bucket_list[idx].above = histogram->bucket_list[idx - 1].above;
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

  init_bucket(right, left, far_right, histogram->alpha);
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

  recompute_bound(left->below, left);
  recompute_bound(left, right);
  recompute_bound(right, right->above);

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

void assert_consistent(struct decaying_histogram *histogram) {
  int idx;
  struct bucket *bucket;

  assert(histogram->num_buckets <= histogram->max_num_buckets);
  assert(histogram->num_buckets > 0);

  for (idx = 0; idx < histogram->num_buckets; idx++)
    decay(histogram, &histogram->bucket_list[idx]);
  for (idx = 0; idx < histogram->num_buckets - 1; idx++) {
    bucket = &histogram->bucket_list[idx];
    recompute_bound(bucket, &histogram->bucket_list[idx + 1]);
    assert(
        bucket->lower_bound <= bucket->mu &&
        bucket->mu <= bucket->upper_bound);
  }

  recompute_bound(NULL, &histogram->bucket_list[0]);
  recompute_bound(&histogram->bucket_list[histogram->num_buckets - 1], NULL);

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

