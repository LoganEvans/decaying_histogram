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
#ifndef CEIL
#define CEIL(x, y) ((int)((x) + ((int)(y) - 1)) / (int)(y))
#endif


static void delete_bucket(
    struct decaying_histogram *histogram, uint32_t bucket_idx);
static void split_bucket(
    struct decaying_histogram *histogram, uint32_t bucket_idx);
static double ipow(double coefficient, uint64_t power);
static double get_decay(
    struct decaying_histogram *histogram, uint64_t missed_generations);
static double compute_bound(
    struct decaying_histogram *histogram,
    struct bucket *left, struct bucket *right);
static double compute_count(
    struct decaying_histogram *histogram, struct bucket *bucket,
    uint64_t generation);
static bool perform_add(
    struct decaying_histogram *histogram, struct bucket *bucket,
    double observation, bool recheck_left_boundary,
    bool recheck_right_boundary, int mp_flag);
static bool find_bucket_helper(
    struct bucket *bucket, double observation, int mp_flag);
//static void assert_consistent(struct decaying_histogram *histogram);


void init_bucket(struct bucket *to_init) {
  to_init->count = 0.0;
  to_init->mu = 0.0;
  to_init->update_generation = 0;
  to_init->below = NULL;
  to_init->above = NULL;

  to_init->boundary_mtx =
    (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));
  pthread_mutex_init(to_init->boundary_mtx, NULL);
}

/*
 * Returns true if the observation is between the lower bucket mu and the
 * current bucket mu.
 */
bool is_target_boundary(struct bucket *bucket, double observation) {
  double left_mu, mu;
  if (bucket->below == NULL) {
    left_mu = observation - 1.0;
  } else {
    left_mu = bucket->below->mu;
  }

  mu = bucket->mu;
  if (bucket->above == NULL) {
    mu = observation + 1.0;
  } else if (left_mu == mu && observation == mu) {
    mu = observation + 1.0;
  }

  if (left_mu <= observation && observation < mu) {
    return true;
  } else {
    return false;
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
  if (bucket == NULL || bucket->update_generation == generation)
    return;

  bucket->count *= get_decay(
      histogram, generation - bucket->update_generation);
  bucket->update_generation = generation;
}

/*
 * If generation == 0, this will provide the density estimate using the last
 * update generation as to compute the count.
 */
double density(
    struct decaying_histogram *histogram, struct bucket *bucket,
    uint64_t generation,
    double *lower_bound_output, double *upper_bound_output) {
  double lower_bound, upper_bound, retval;

  lower_bound = compute_bound(histogram, bucket->below, bucket);
  upper_bound = compute_bound(histogram, bucket, bucket->above);

  if (lower_bound_output != NULL)
    *lower_bound_output = lower_bound;
  if (upper_bound_output != NULL)
    *upper_bound_output = upper_bound;

  if (bucket->update_generation <= generation) {
    decay(histogram, bucket, generation);
    retval = bucket->count / total_count(histogram, generation);
  } else {
    retval =
        bucket->count / total_count(histogram, bucket->update_generation);
  }
  retval /= upper_bound - lower_bound;

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
  if (missed_generations < (uint32_t)histogram->max_num_buckets)
    return histogram->pow_table[missed_generations];
  else
    return ipow(1.0 - histogram->alpha, missed_generations);
}

double compute_count(
    struct decaying_histogram *histogram, struct bucket *bucket,
    uint64_t generation) {
  return (
      bucket->count *
      get_decay(histogram, generation - bucket->update_generation));
}

// We can't increment generation until we finally make our update, but
// since the boundaries don't move upon decay, we can use any generation
// value to compute the bounds. We'll use the largest generation so that
// we don't need to recompute the decay for that bucket.
double compute_bound(
    struct decaying_histogram *histogram,
    struct bucket *left, struct bucket *right) {
  uint64_t generation;
  double count_left, mu_left, count_right, mu_right;

  if (!left) {
    return right->mu - 0.5;
  } else if (!right) {
    return left->mu + 0.5;
  } else if (left->update_generation > right->update_generation) {
    generation = left->update_generation;
    count_left = left->count;
    mu_left = left->mu;
    count_right = compute_count(histogram, right, generation);
    mu_right = right->mu;
  } else {
    generation = right->update_generation;
    count_left = compute_count(histogram, left, generation);
    mu_left = left->mu;
    count_right = right->count;
    mu_right = right->mu;
  }

  return (
      (mu_left * count_left + mu_right * count_right) /
      (count_left + count_right));
}

bool perform_add(
    struct decaying_histogram *histogram, struct bucket *bucket,
    double observation, bool recheck_left_boundary,
    bool recheck_right_boundary, int mp_flag) {
  double boundary;
  uint64_t update_generation;

  if (recheck_left_boundary) {
    boundary = compute_bound(histogram, bucket->below, bucket);
    if (observation < boundary) {
      if (mp_flag & DHIST_MULTI_THREADED)
        pthread_mutex_unlock(&histogram->generation_mutex);
      return false;
    }
  } else if (recheck_right_boundary) {
    boundary = compute_bound(histogram, bucket, bucket->above);
    if (boundary < observation) {
      if (mp_flag & DHIST_MULTI_THREADED)
        pthread_mutex_unlock(&histogram->generation_mutex);
      return false;
    }
  }

  if (mp_flag & DHIST_MULTI_THREADED) {
    pthread_mutex_lock(&histogram->generation_mutex);
    update_generation = ++histogram->generation;
    pthread_mutex_unlock(&histogram->generation_mutex);
  } else {
    update_generation = ++histogram->generation;
  }

  bucket->count = 1.0 + compute_count(
      histogram, bucket, update_generation);
  bucket->mu =
      (bucket->count * bucket->mu + observation) / (bucket->count + 1);
  bucket->update_generation = update_generation;

  return true;
}

void init_decaying_histogram(
    struct decaying_histogram *histogram, int target_buckets, double alpha) {
  double max_total_count, expected_count;
  uint64_t idx;

  max_total_count = 1.0 / alpha;
  expected_count = 1.0 / (alpha * target_buckets);

  // TODO: What are good thresholds?
  histogram->delete_bucket_threshold = expected_count * (1.0 / 4.0);
  histogram->split_bucket_threshold = expected_count * (7.0 / 4.0);
  histogram->alpha = alpha;
  histogram->max_num_buckets = 2 *
      (uint32_t)CEIL(max_total_count, histogram->delete_bucket_threshold);
  histogram->bucket_list = (struct bucket *)malloc(
      sizeof(struct bucket) * histogram->max_num_buckets);
  histogram->num_buckets = 1;
  for (idx = 0; idx < histogram->max_num_buckets; idx++)
    init_bucket(&histogram->bucket_list[idx]);
  histogram->generation = 0;

  histogram->pow_table = (double *)malloc(
      histogram->max_num_buckets * sizeof(double));
  for (idx = 0; idx < histogram->max_num_buckets; idx++)
    histogram->pow_table[idx] = ipow(1.0 - alpha, idx);

  pthread_rwlock_init(&histogram->rwlock, NULL);
  pthread_mutex_init(&histogram->generation_mutex, NULL);

  return;
}

double total_count(struct decaying_histogram *histogram, uint64_t generation) {
  return (1 - get_decay(histogram, generation)) /
         histogram->alpha;
}

void clean_decaying_histogram(struct decaying_histogram *histogram) {
  free(histogram->bucket_list);
  free(histogram->pow_table);
}

bool find_bucket_helper(
    struct bucket *bucket, double observation, int mp_flag) {
  if (mp_flag & DHIST_SINGLE_THREADED) {
    return true;
  } else {
    pthread_mutex_lock(bucket->boundary_mtx);

    if (!is_target_boundary(bucket, observation)) {
      // If we raced, restart.
      pthread_mutex_unlock(bucket->boundary_mtx);
      return false;
    } else {
      return true;
    }
  }
}

struct bucket * find_bucket(
    struct decaying_histogram *histogram, double observation, int mp_flag) {
  uint32_t low, mid, high;
  struct bucket *bucket;
  int tries = 0, retries = 0, loops = 0;

  while (histogram->num_buckets == 1) {
    bucket = &histogram->bucket_list[0];
    if (find_bucket_helper(bucket, observation, mp_flag))
      return bucket;
  }

  do {
    low = 0;
    high = histogram->num_buckets;

    while (low < high) {
      mid = (low + high) / 2;
      bucket = &histogram->bucket_list[mid];
      if (is_target_boundary(bucket, observation)) {
        if (find_bucket_helper(bucket, observation, mp_flag)) {
          return bucket;
        } else {
          // We raced, so start over.
          retries++;
          low = 0;
          high = histogram->num_buckets;
        }
      } else if (bucket->mu <= observation) {
        low = mid + 1;
      } else {
        high = mid;
      }
      tries++;
    }
    loops++;
    // We can fail to find a bucket if the bucket boundaries moved around
    // during the search. Start over.
  } while (1);
  return NULL;
}

void full_refresh(struct decaying_histogram *histogram, int mp_flag) {
  uint32_t idx;
  bool do_another_round, finished_deletes;
  struct bucket *bucket;

  if (mp_flag & DHIST_MULTI_THREADED)
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

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_rwlock_unlock(&histogram->rwlock);
}

void add_observation(
    struct decaying_histogram *histogram, double observation, int mp_flag) {
  struct bucket *bucket;
  bool add_succeeded;
  double boundary;

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_rwlock_rdlock(&histogram->rwlock);

  do {
    bucket = find_bucket(histogram, observation, mp_flag);

    // We have a read lock on the bucket that we need to update, but we don't
    // have a write lock until we have both boundary mutexes. If we can't grab
    // the other lock, we might be in danger of a deadlock. To avoid that, if
    // we fail to take the other lock and that lock is below->boundary_mtx,
    // we'll drop our lock and try again.

    boundary = compute_bound(histogram, bucket->below, bucket);

    if (!bucket->below ) {
      boundary = observation;
    } else if (!bucket->above ) {
      boundary = observation;
    }

    if (mp_flag & DHIST_SINGLE_THREADED) {
      // Single threaded, handles both cases.
      add_succeeded = perform_add(
          histogram, (boundary <= observation) ? bucket : bucket->below,
          observation, false, false, mp_flag);
    } else if (observation < boundary) {
      // We need to insert into bucket->below. That bucket exists because we
      // would have created an artificial boundary if it didn't.
      if (pthread_mutex_trylock(bucket->below->boundary_mtx) == 0) {
        add_succeeded = perform_add(
            histogram, bucket->below, observation, true, false, mp_flag);
        pthread_mutex_unlock(bucket->below->boundary_mtx);
      } else {
        add_succeeded = false;
      }
    } else if (bucket->above) {
      // We need to insert into this bucket, and since we have a neighbor to
      // the right, we need to extend our lock.
      pthread_mutex_lock(bucket->above->boundary_mtx);
      add_succeeded = perform_add(
          histogram, bucket, observation, false, true, mp_flag);
      pthread_mutex_unlock(bucket->above->boundary_mtx);
    } else {
      // No bucket exists to the right, and since we already have the boundary
      // mutex for this bucket, nothing else can add a bucket to the right.
      add_succeeded = perform_add(
          histogram, bucket, observation, false, false, mp_flag);
    }

    if (mp_flag & DHIST_MULTI_THREADED)
      pthread_mutex_unlock(bucket->boundary_mtx);
  } while (!add_succeeded);

  if (mp_flag & DHIST_MULTI_THREADED)
      pthread_rwlock_unlock(&histogram->rwlock);

  if (bucket->count < histogram->delete_bucket_threshold)
    full_refresh(histogram, mp_flag);
  if (bucket->count > histogram->split_bucket_threshold)
    full_refresh(histogram, mp_flag);
}

void print_histogram(
    struct decaying_histogram *histogram, bool estimate_ok,
    const char *title, const char *xlabel, int mp_flag) {
  uint32_t idx, num_buckets;
  uint64_t generation;
  double *boundaries, *densities;
  struct bucket *left, *right;

  if (mp_flag & DHIST_SINGLE_THREADED) {
    /* pass */
  } else if (!estimate_ok) {
    pthread_rwlock_wrlock(&histogram->rwlock);
  } else {
    pthread_rwlock_rdlock(&histogram->rwlock);
  }

  num_buckets = histogram->num_buckets;
  densities = (double *)malloc(sizeof(double) * num_buckets);
  boundaries = (double *)malloc(sizeof(double) * (num_buckets + 1));

  if (mp_flag & DHIST_MULTI_THREADED) {
    pthread_mutex_lock(&histogram->generation_mutex);
    generation = histogram->generation;
    pthread_mutex_unlock(&histogram->generation_mutex);
  } else {
    generation = histogram->generation;
  }

  right = NULL;
  for (idx = 0; idx < num_buckets; idx++) {
    left = right;
    right = &histogram->bucket_list[idx];
    if (mp_flag & DHIST_MULTI_THREADED) {
      pthread_mutex_lock(right->boundary_mtx);
      densities[idx] = density(
          histogram, right, generation, &boundaries[idx], NULL);
      pthread_mutex_unlock(right->boundary_mtx);
    } else {
      densities[idx] = density(
          histogram, right, generation, &boundaries[idx], NULL);
    }
  }

  if (mp_flag & DHIST_MULTI_THREADED) {
    pthread_mutex_lock(right->boundary_mtx);
    boundaries[idx] = compute_bound(histogram, right, NULL);
    pthread_mutex_unlock(right->boundary_mtx);
    pthread_rwlock_unlock(&histogram->rwlock);
  } else {
    boundaries[idx] = compute_bound(histogram, right, NULL);
  }

  printf("{");
  if (title)
    printf("\"title\": \"%s\", ", title);
  if (xlabel)
    printf("\"xlabel\": \"%s\", ", xlabel);

  printf("\"generation\": %lu, ", generation);

  printf("\"densities\": [");
  for (idx = 0; idx < num_buckets - 1; idx++) {
    printf("%lf, ", densities[idx]);
  }
  printf("%lf], ", densities[idx]);

  printf("\"boundaries\": [");
  for (idx = 0; idx < num_buckets; idx++) {
    printf("%lf, ", boundaries[idx]);
  }
  printf("%lf]}\n", boundaries[idx]);

  free(densities);
  free(boundaries);
}

/*
 * All buckets must be decayed before calling this.
 */
void delete_bucket(
    struct decaying_histogram *histogram, uint32_t bucket_idx) {
  uint32_t lucky_idx, idx;
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
    histogram->bucket_list[idx].update_generation =
        histogram->bucket_list[idx + 1].update_generation;
  }
  --histogram->num_buckets;
  histogram->bucket_list[histogram->num_buckets].below = NULL;
  histogram->bucket_list[histogram->num_buckets - 1].above = NULL;

  return;
}

void split_bucket(
    struct decaying_histogram *histogram, uint32_t bucket_idx) {
  uint32_t idx;
  double lower_bound, upper_bound, diameter, median;
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
    histogram->bucket_list[idx].update_generation =
        histogram->bucket_list[idx - 1].update_generation;
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
    lower_bound = compute_bound(histogram, left->below, left);
    upper_bound = compute_bound(histogram, left, right);

    diameter = upper_bound - lower_bound;
    median = lower_bound + diameter / 2.0;
    left->mu = median - diameter / 6.0;
    right->mu = median + diameter / 6.0;
  }

  right->update_generation = left->update_generation;

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

