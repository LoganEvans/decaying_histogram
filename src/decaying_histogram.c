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

#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif

static void delete_bucket(
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


struct bucket * init_bucket(int name) {
  struct bucket * bucket;

  bucket = (struct bucket *)malloc(sizeof(struct bucket));
  memset(bucket, 0, sizeof(struct bucket));
  bucket->height = 1;
  bucket->boundary_mtx =
    (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));
  pthread_mutex_init(bucket->boundary_mtx, NULL);
  bucket->name = name;

  return bucket;
}

void destroy_bucket(struct bucket *bucket) {
  free(bucket->boundary_mtx);
  free(bucket);
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
 * This is a bit of a misnomer... this doesn't return the height of the bar,
 * but instead returns the area of the bar.
 *
 * If generation == 0, this will provide the density estimate using the last
 * update generation as to compute the count.
 */
double weight(
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
  histogram->bucket_list = (struct bucket **)malloc(
      sizeof(struct bucket*) * histogram->max_num_buckets);

  for (idx = 0; idx < histogram->max_num_buckets; idx++)
    histogram->bucket_list[idx] = NULL;

  histogram->bucket_list[0] = init_bucket(histogram->namer++);
  histogram->num_buckets = 1;

  histogram->generation = 0;

  histogram->pow_table = (double *)malloc(
      histogram->max_num_buckets * sizeof(double));
  for (idx = 0; idx < histogram->max_num_buckets; idx++)
    histogram->pow_table[idx] = ipow(1.0 - alpha, idx);

  pthread_rwlock_init(&histogram->rwlock, NULL);
  pthread_mutex_init(&histogram->generation_mutex, NULL);

  histogram->namer = 10;

  return;
}

double total_count(struct decaying_histogram *histogram, uint64_t generation) {
  return (1 - get_decay(histogram, generation)) /
         histogram->alpha;
}

void clean_decaying_histogram(struct decaying_histogram *histogram) {
  uint32_t idx;

  for (idx = 0; idx < histogram->max_num_buckets; idx++)
    free(histogram->bucket_list[idx]);

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
    bucket = histogram->bucket_list[0];
    if (find_bucket_helper(bucket, observation, mp_flag))
      return bucket;
  }

  do {
    low = 0;
    high = histogram->num_buckets;

    while (low < high) {
      mid = (low + high) / 2;
      bucket = histogram->bucket_list[mid];
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
    decay(histogram, histogram->bucket_list[idx], histogram->generation);

  /*
   * In order to make sure we don't overrun the number of buckets available,
   * we need to condense before we expand.
   */
  finished_deletes = false;
  do {
    do_another_round = false;
    for (idx = 0; idx < histogram->num_buckets; idx++) {
      bucket = histogram->bucket_list[idx];
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

  for (int i = 0; i < histogram->num_buckets - 1; i++) {
    assert(histogram->bucket_list[i]->above == histogram->bucket_list[i + 1]);
    assert(histogram->bucket_list[i]->mu <= histogram->bucket_list[i + 1]->mu);
  }
  for (int i = 1; i < histogram->num_buckets; i++) {
    assert(histogram->bucket_list[i]->below == histogram->bucket_list[i - 1]);
  }

  if (bucket->count < histogram->delete_bucket_threshold && histogram->num_buckets > 2) {
    full_refresh(histogram, mp_flag);
    //printf("delete\n");
    //print_histogram(histogram, true, NULL, NULL, mp_flag);
  }
  if (bucket->count > histogram->split_bucket_threshold) {
    //printf("split\n");
    full_refresh(histogram, mp_flag);
    //print_histogram(histogram, true, NULL, NULL, mp_flag);
  }

  for (int i = 0; i < histogram->num_buckets - 1; i++) {
    assert(histogram->bucket_list[i]->above == histogram->bucket_list[i + 1]);
    assert(histogram->bucket_list[i]->mu <= histogram->bucket_list[i + 1]->mu);
  }
  for (int i = 1; i < histogram->num_buckets; i++) {
    assert(histogram->bucket_list[i]->below == histogram->bucket_list[i - 1]);
  }
}

void print_histogram(
    struct decaying_histogram *histogram, bool estimate_ok,
    const char *title, const char *xlabel, int mp_flag) {
  uint32_t idx, num_buckets;
  uint64_t generation;
  double *boundaries, *weights;
  struct bucket *left, *right;

  if (mp_flag & DHIST_SINGLE_THREADED) {
    /* pass */
  } else if (!estimate_ok) {
    pthread_rwlock_wrlock(&histogram->rwlock);
  } else {
    pthread_rwlock_rdlock(&histogram->rwlock);
  }

  num_buckets = histogram->num_buckets;
  weights = (double *)malloc(sizeof(double) * num_buckets);
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
    right = histogram->bucket_list[idx];
    if (mp_flag & DHIST_MULTI_THREADED) {
      pthread_mutex_lock(right->boundary_mtx);
      weights[idx] = weight(
          histogram, right, generation, &boundaries[idx], NULL);
      pthread_mutex_unlock(right->boundary_mtx);
    } else {
      weights[idx] = weight(
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

  printf("\"weights\": [");
  for (idx = 0; idx < num_buckets - 1; idx++) {
    printf("%lf, ", weights[idx]);
  }
  printf("%lf], ", weights[idx]);

  printf("\"boundaries\": [");
  for (idx = 0; idx < num_buckets; idx++) {
    printf("%lf, ", boundaries[idx]);
  }
  printf("%lf]}\n", boundaries[idx]);

  free(weights);
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

  dying_bucket = histogram->bucket_list[bucket_idx];
  if (bucket_idx == 0)
    lucky_idx = bucket_idx + 1;
  else if (bucket_idx == histogram->num_buckets - 1)
    lucky_idx = bucket_idx - 1;
  else if (dying_bucket->below->count < dying_bucket->above->count)
    lucky_idx = bucket_idx - 1;
  else
    lucky_idx = bucket_idx + 1;

  lucky_bucket = histogram->bucket_list[lucky_idx];
  lucky_bucket->mu =
      (lucky_bucket->mu * lucky_bucket->count +
       dying_bucket->mu * dying_bucket->count) /
      (lucky_bucket->count + dying_bucket->count);
  lucky_bucket->count += dying_bucket->count;
  if (dying_bucket->below == lucky_bucket) {
    lucky_bucket->above = dying_bucket->above;
    if (lucky_bucket->above)
      lucky_bucket->above->below = lucky_bucket;
  } else {
    lucky_bucket->below = dying_bucket->below;
    if (lucky_bucket->below)
      lucky_bucket->below->above = lucky_bucket;
  }
  destroy_bucket(dying_bucket);

  // Shift everything left.
  for (idx = bucket_idx; idx < histogram->num_buckets - 1; idx++) {
    histogram->bucket_list[idx] = histogram->bucket_list[idx + 1];
  }
  histogram->bucket_list[histogram->num_buckets] = NULL;
  --histogram->num_buckets;

  return;
}

#if 0
void split_bucket(
    struct decaying_histogram *histogram, uint32_t bucket_idx) {
  uint32_t idx;
  double lower_bound, upper_bound, diameter, median;
  struct bucket *left, *right;

  if (histogram->num_buckets == 2 &&
      histogram->bucket_list[0]->mu == histogram->bucket_list[1]->mu) {
    // This will happen if we observe a stream of constants. Let's avoid
    // making more than two buckets until we have observed two unique values.
    return;
  }

  right = histogram->bucket_list[bucket_idx];
  lower_bound = compute_bound(histogram, right->below, right);
  upper_bound = compute_bound(histogram, right, right->above);

  // Shift everything right.
  for (idx = histogram->num_buckets; idx > bucket_idx; idx--) {
    histogram->bucket_list[idx] = histogram->bucket_list[idx - 1];
  }

  histogram->bucket_list[bucket_idx] = init_bucket();
  ++histogram->num_buckets;
  left = histogram->bucket_list[bucket_idx];

  left->below = right->below;
  if (left->below)
    left->below->above = left;
  left->above = right;
  right->below = left;
  right->count /= 2.0;
  left->count = right->count;
  if (histogram->num_buckets == 2) {
    // This is an awkward split because we don't have enough information to
    // compute a diameter. Instead, we can set mu to the same value in both
    // buckets, which will also make the boundary between them be mu. As soon
    // as we observe a non-mu entry, the split between the buckets
    // will be more sensical. However, we will need to take care to not split
    // a bucket again until this happens.
    left->mu = right->mu;
  } else {
    diameter = upper_bound - lower_bound;
    median = lower_bound + diameter / 2.0;
    left->mu = median - diameter / 6.0;
    right->mu = median + diameter / 6.0;
  }
  assert(left->mu <= right->mu);
  if (left->below)
    assert(left->below->mu <= left->mu);
  if (right->above)
    assert(right->mu <= right->above->mu);

  right->update_generation = left->update_generation;

  return;
}
#endif

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

int count_nodes(struct bucket *root) {
  if (root == NULL) {
    return 0;
  } else {
    return 1 + count_nodes(root->children[0]) + count_nodes(root->children[1]);
  }
}

int assert_invariant(struct bucket *root) {
  int left, right;

  if (root == NULL) {
    return 0;
  } else if (root->children[0] == NULL && root->children[1] == NULL) {
    assert(root->height == 1);
    return 1;
  } else {
    left = assert_invariant(root->children[0]);
    right = assert_invariant(root->children[1]);
    if (root->height != left + 1 && root->height != right + 1 ||
        (root->height <= left || root->height <= right)) {
      printf(
        "root height is not correct. heights -- root: %d left: %d right: %d\n",
        root->height, left, right);
      assert(false);
    }

    if (ABS(right - left) > 1) {
      printf("invariant not met. heights -- left: %d right %d\n",
        left, right);
      assert(false);
    }

    return 1 + (left > right ? left : right);
  }
}

void fix_height(struct bucket *bucket) {
  struct bucket *left, *right;

  left = bucket->children[0];
  right = bucket->children[1];

  if (left == NULL && right == NULL) {
    bucket->height = 1;
  } else if (left == NULL) {
    bucket->height = 1 + right->height;
  } else if (right == NULL) {
    bucket->height = 1 + left->height;
  } else {
    bucket->height =
        1 + (left->height > right->height ? left->height : right->height);
  }
}

struct bucket * rotate_single(struct bucket *root, int dir) {
  struct bucket *new_root;
  int dir_from_parent;
  printf(" > rotate_single(%d, %d)\n", root->name, dir);

  new_root = root->children[!dir];
  new_root->parent = root->parent;

  root->children[!dir] = new_root->children[dir];
  if (new_root->children[dir])
    new_root->children[dir]->parent = root->children[!dir];
  new_root->children[dir] = root;
  root->parent = new_root->children[dir];

  if (new_root->parent != NULL) {
    if (new_root->parent->children[0] == root)
      dir_from_parent = 0;
    else
      dir_from_parent = 1;
    new_root->parent->children[dir_from_parent] = new_root;
  }

  fix_height(root);
  fix_height(new_root);

  return new_root;
}

struct bucket * rotate_double(struct bucket *root, int dir) {
  struct bucket *new_root;
  int dir_from_parent;
  printf(" > rotate_double(%d, %d)\n", root->name, dir);

  new_root = root->children[!dir]->children[dir];
  new_root->parent = root->parent;

  root->children[!dir]->children[dir] = new_root->children[!dir];
  new_root->children[!dir]->parent = root->children[!dir]->children[dir];

  new_root->children[!dir] = root->children[!dir];
  root->children[!dir]->parent = new_root->children[!dir];

  root->children[!dir] = new_root->children[dir];
  if (new_root->children[dir])
    new_root->children[dir]->parent = root->children[!dir];

  new_root->children[dir] = root;
  root->parent = new_root->children[dir];

  if (new_root->parent != NULL) {
    if (new_root->parent->children[0] == root)
      dir_from_parent = 0;
    else
      dir_from_parent = 1;
    new_root->parent->children[dir_from_parent] = new_root;
  }

  fix_height(new_root->children[0]);
  fix_height(new_root->children[1]);
  fix_height(new_root);

  return new_root;
}

int compute_balance(struct bucket *bucket) {
  struct bucket *left, *right;
  int dir;
  int heights[2];

  for (dir = 0; dir <= 1; dir++) {
    if (bucket->children[dir] == NULL)
      heights[dir] = 0;
    else
      heights[dir] = bucket->children[dir]->height;
  }

  return heights[1] - heights[0];
}

void fix_balance(struct decaying_histogram *histogram, struct bucket *bucket) {
  int balance;
  int prior_height;
  int dir, dir_from_parent;
  printf(" > fix_balance(~, %d)\n", bucket->name);

  while (bucket != NULL) {
    fix_height(bucket);
    prior_height = bucket->height;
    balance = compute_balance(bucket);
    printf(" <<loop>> (%d | %d, %d)\n",
      bucket->name, bucket->children[0] ? bucket->children[0]->height : 0,
      bucket->children[1] ? bucket->children[1]->height : 0);
    printf("<<<<<<<<<<<\n");
    print_tree(histogram->root);
    printf(">>>>>>>>>>>\n");
    while (ABS(compute_balance(bucket)) >= 2) {
      printf("balance: %d\n", compute_balance(bucket));
      if (bucket->children[0] == NULL) {
        dir = 0;
      } else if (bucket->children[1] == NULL) {
        dir = 1;
      } else if (bucket->children[0]->height < bucket->children[1]->height) {
        dir = 0;
      } else {
        dir = 1;
      }
      printf("rebalance! %d\n", dir);

      if (bucket->children[!dir]->children[!dir] == NULL) {
        bucket = rotate_double(bucket, dir);
        fix_balance(histogram, bucket->children[0]);
        fix_balance(histogram, bucket->children[1]);
      } else if (bucket->children[!dir]->children[dir] == NULL) {
        bucket = rotate_single(bucket, dir);
      } else if (bucket->children[!dir]->children[dir]->height >
          bucket->children[!dir]->children[!dir]->height) {
        bucket = rotate_double(bucket, dir);
        fix_balance(histogram, bucket->children[0]);
        fix_balance(histogram, bucket->children[1]);
      } else {
        bucket = rotate_single(bucket, dir);
      }

      if (bucket->parent == NULL) {
        histogram->root = bucket;
      } else {
        if (bucket->parent->children[0] == bucket) {
          dir_from_parent = 0;
        } else {
          dir_from_parent = 1;
        }

        bucket->parent->children[dir_from_parent] = bucket;
      }
    }

    if (bucket->height != prior_height) {
      bucket = bucket->parent;
    } else {
      break;
    }
  }
}

struct bucket * split_bucket(
    struct decaying_histogram *histogram, struct bucket *bucket) {
  int dir;
  struct bucket *new_bucket;

  new_bucket = init_bucket(histogram->namer++);
  printf(" > split_bucket %d => %d\n", bucket->name, new_bucket->name);
  new_bucket->parent = bucket;

  if (compute_balance(bucket) <= 0) {
    dir = 0;
  } else {
    dir = 1;
  }

  new_bucket->children[dir] = bucket->children[dir];
  bucket->children[dir] = new_bucket;

  fix_balance(histogram, new_bucket);
  fix_balance(histogram, bucket);
}

void _print_tree(struct bucket *bucket, int depth) {
  printf("!!!!%d\n", depth);

  printf("%d(%d)", bucket->name, bucket->height);

  if (bucket->children[0]) {
    printf("\t >");
    _print_tree(bucket->children[0], depth + 1);
  }

  if (bucket->children[1]) {
    printf("\n");
    for (int i = 0; i < depth; i++)
      printf("\t");
    printf("\t\\>");
    _print_tree(bucket->children[1], depth + 1);
  }
}

void print_tree(struct bucket *bucket) {
  _print_tree(bucket, 0);
  printf("\n");
}

