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
#include <string.h>
// XXX
#include <math.h>

#include "decaying_histogram.h"

// The same as ceil(x / y). Using this so that math.h is not a dependency.
#ifndef CEIL
#define CEIL(x, y) ((int)((x) + ((int)(y) - 1)) / (int)(y))
#endif

#ifndef ABS
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#endif

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
static void set_child(struct bucket *root, struct bucket *child, int dir);
static void fix_height(struct bucket *bucket);
static void fix_balance(
    struct decaying_histogram *histogram, struct bucket *bucket, bool upward);
static struct bucket * rotate_single(
    struct decaying_histogram *histogram, struct bucket *root, int dir);
static struct bucket * rotate_double(
    struct decaying_histogram *histogram, struct bucket *root, int dir);
static int compute_balance(struct bucket *bucket);
static void _print_tree(struct bucket *bucket, int depth);
static void assert_consistent(struct decaying_histogram *histogram);
static bool is_target_boundary(struct bucket *bucket, double observation);
static void obtain_write_lock(struct bucket *bucket, int mp_flag, int line);
static void release_write_lock(struct bucket *bucket, int mp_flag);
static void lock_boundary(struct bucket *bucket, int mp_flag, int line);
static bool trylock_boundary_succeeded(
    struct bucket *bucket, int mp_flag, int lint);
static void unlock_boundary(struct bucket *bucket, int mp_flag);
static uint64_t get_generation(
    struct decaying_histogram *histogram, bool increment, int mp_flag);


struct bucket * init_bucket(
    struct decaying_histogram *histogram, int mp_flag) {
  struct bucket * bucket;
  int idx;

  for (idx = 0; idx < histogram->max_num_buckets; idx++) {
    bucket = &histogram->bucket_list[idx];
    if (bucket->is_enabled)
      continue;
    if (mp_flag & DHIST_SINGLE_THREADED)
      break;
    if (!trylock_boundary_succeeded(bucket, mp_flag, __LINE__))
      continue;
    if (bucket->is_enabled) {
      // We raced or had a stale reference to bucket->is_enabled.
      unlock_boundary(bucket, mp_flag);
      continue;
    }
    break;
  }

  bucket->is_enabled = true;
  bucket->count = 0.0;
  bucket->mu = 0.0;
  bucket->update_generation = 0;
  bucket->below = NULL;
  bucket->above = NULL;
  bucket->parent = NULL;
  bucket->children[0] = NULL;
  bucket->children[1] = NULL;
  bucket->height = NULL;
  unlock_boundary(bucket, mp_flag);

  return bucket;
}

void recycle_bucket(struct bucket *bucket, int mp_flag) {
  bucket->is_enabled = false;
  release_write_lock(bucket, mp_flag);
}

void lock_boundary(struct bucket *bucket, int mp_flag, int line) {
  //printf("lock_boundary\n");
  if (mp_flag & DHIST_MULTI_THREADED) {
    pthread_mutex_lock(bucket->boundary_mtx);
    assert(bucket->lock_held == false);
    bucket->lock_held = true;
    bucket->line = line;
  }
}

bool trylock_boundary_succeeded(struct bucket *bucket, int mp_flag, int line) {
  if (mp_flag & DHIST_SINGLE_THREADED) {
    //printf("trylock_boundary\n");
    return true;
  } else if (pthread_mutex_trylock(bucket->boundary_mtx) == 0) {
    assert(bucket->lock_held == false);
    bucket->lock_held = true;
    bucket->line = line;
    return true;
  } else {
    return false;
  }
}

void unlock_boundary(struct bucket *bucket, int mp_flag) {
  //printf("unlock_boundary\n");
  if (mp_flag & DHIST_MULTI_THREADED) {
    assert(bucket->lock_held == true);
    bucket->lock_held = false;
    pthread_mutex_unlock(bucket->boundary_mtx);
  }
}

static uint64_t get_generation(
    struct decaying_histogram *histogram, bool increment, int mp_flag) {
  uint64_t generation;

  if (mp_flag & DHIST_MULTI_THREADED) {
    pthread_mutex_lock(&histogram->generation_mtx);
    if (increment)
      ++histogram->generation;
    generation = histogram->generation;
    pthread_mutex_unlock(&histogram->generation_mtx);
  } else {
    if (increment)
      ++histogram->generation;
    generation = histogram->generation;
  }

  return generation;
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
    if (observation < boundary)
      return false;
  } else if (recheck_right_boundary) {
    boundary = compute_bound(histogram, bucket, bucket->above);
    if (boundary < observation)
      return false;
  }

  update_generation = get_generation(histogram, true, mp_flag);

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
  struct bucket *bucket;

  max_total_count = 1.0 / alpha;
  expected_count = 1.0 / (alpha * target_buckets);

  // TODO: What are good thresholds?
  histogram->delete_bucket_threshold = expected_count * (1.0 / 4.0);
  histogram->split_bucket_threshold = expected_count * (7.0 / 4.0);
  histogram->alpha = alpha;
  histogram->max_num_buckets =
      (uint32_t)CEIL(max_total_count, histogram->delete_bucket_threshold);
  histogram->bucket_list = (struct bucket *)malloc(
      histogram->max_num_buckets * sizeof(struct bucket));
  for (idx = 0; idx < histogram->max_num_buckets; idx++) {
    bucket = &histogram->bucket_list[idx];
    bucket->is_enabled = false;
    bucket->boundary_mtx =
      (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(bucket->boundary_mtx, NULL);
  }

  // We're single threaded in initialization.
  histogram->root = init_bucket(histogram, DHIST_SINGLE_THREADED);
  histogram->num_buckets = 1;

  histogram->generation = 0;

  histogram->pow_table = (double *)malloc(
      histogram->max_num_buckets * sizeof(double));
  for (idx = 0; idx < histogram->max_num_buckets; idx++)
    histogram->pow_table[idx] = ipow(1.0 - alpha, idx);

  pthread_mutex_init(&histogram->tree_mtx, NULL);
  pthread_mutex_init(&histogram->generation_mtx, NULL);

  return;
}

double total_count(struct decaying_histogram *histogram, uint64_t generation) {
  return (1 - get_decay(histogram, generation)) /
         histogram->alpha;
}

void clean_decaying_histogram(struct decaying_histogram *histogram) {

  free(histogram->pow_table);
}

bool find_bucket_helper(
    struct bucket *bucket, double observation, int mp_flag) {
  if (mp_flag & DHIST_SINGLE_THREADED) {
    return true;
  } else {
    lock_boundary(bucket, mp_flag, __LINE__);

    if (!is_target_boundary(bucket, observation) || !bucket->is_enabled) {
      // If we raced, restart.
      unlock_boundary(bucket, mp_flag);
      return false;
    } else {
      return true;
    }
  }
}

struct bucket * find_bucket(
    struct decaying_histogram *histogram, double observation, int mp_flag) {
  struct bucket *bucket;
  int tries = 0, retries = 0;

  do {
    bucket = histogram->root;
    do {
      if (is_target_boundary(bucket, observation)) {
        if (find_bucket_helper(bucket, observation, mp_flag)) {
          return bucket;
        } else {
          // We raced, so start over.
          retries++;
          bucket = histogram->root;
        }
      }
      bucket = bucket->children[bucket->mu <= observation];
      tries++;
    } while (bucket != NULL);
    // We can fail to find a bucket if the bucket boundaries moved around
    // during the search. Start over.
  } while (1);
  return NULL;
}

void add_observation(
    struct decaying_histogram *histogram, double observation, int mp_flag) {
  struct bucket *bucket;
  bool add_succeeded;
  double boundary;

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
      if (trylock_boundary_succeeded(bucket->below, mp_flag, __LINE__)) {
        add_succeeded = perform_add(
            histogram, bucket->below, observation, true, false, mp_flag);
        unlock_boundary(bucket->below, mp_flag);
      } else {
        add_succeeded = false;
      }
    } else if (bucket->above) {
      // We need to insert into this bucket, and since we have a neighbor to
      // the right, we need to extend our lock.
      lock_boundary(bucket->above, mp_flag, __LINE__);
      add_succeeded = perform_add(
          histogram, bucket, observation, false, true, mp_flag);
      unlock_boundary(bucket->above, mp_flag);
    } else {
      // No bucket exists to the right, and since we already have the boundary
      // mutex for this bucket, nothing else can add a bucket to the right.
      add_succeeded = perform_add(
          histogram, bucket, observation, false, false, mp_flag);
    }

    unlock_boundary(bucket, mp_flag);
  } while (!add_succeeded);

  if (bucket->count < histogram->delete_bucket_threshold &&
      histogram->num_buckets > 1) {
    delete_bucket(histogram, bucket, mp_flag);
  }
  if (bucket->count > histogram->split_bucket_threshold) {
    split_bucket(histogram, bucket, mp_flag);
  }
}

void print_histogram(
    struct decaying_histogram *histogram, bool estimate_ok,
    const char *title, const char *xlabel, int mp_flag) {
  uint32_t idx, num_buckets;
  uint64_t generation;
  double *boundaries, *weights;
  struct bucket *bucket;

  weights = (double *)malloc(sizeof(double) * histogram->max_num_buckets);
  boundaries = (double *)malloc(
      sizeof(double) * (histogram->max_num_buckets + 1));

  generation = get_generation(histogram, false, mp_flag);

  bucket = histogram->root;
  while (bucket->children[0])
    bucket = bucket->children[0];

  idx = 0;
  while (1) {
    lock_boundary(bucket, mp_flag, __LINE__);
    weights[idx] = weight(
        histogram, bucket, generation,
        &boundaries[idx], &boundaries[idx + 1]);
    unlock_boundary(bucket, mp_flag);
    if (bucket->above == NULL)
      break;
    bucket = bucket->above;
    idx++;
  }
  num_buckets = idx;

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

void set_child(struct bucket *root, struct bucket *child, int dir) {
  root->children[dir] = child;
  if (child != NULL)
    child->parent = root;
  fix_height(root);
}

int assert_invariant(struct bucket *root) {
  int left, right;

  if (root == NULL) {
    return 0;
  } else if (root->children[0] == NULL && root->children[1] == NULL) {
    assert(root->height == 1);
    return 1;
  } else {
    if (root->children[0] && root->mu < root->children[0]->mu) {
      printf("ORDER ERROR(0): root->mu: %lf < root->children[0]->mu: %lf ...\n",
          root->mu, root->children[0]->mu);
      assert(false);
    }

    if (root->children[1] && root->mu > root->children[1]->mu) {
      printf("ORDER ERROR(1): root->mu: %lf > root->children[1]->mu: %lf ...\n",
          root->mu, root->children[1]->mu);
      assert(false);
    }

    left = assert_invariant(root->children[0]);
    if (root->children[0] == NULL && left != 0) {
      printf("ERROR(1): root->children[0] == NULL && left(%d) != 0\n", left);
      assert(false);
    } else if (root->children[0] && (root->children[0]->height != left)) {
      printf("ERROR(2): root->children[0]->hieght(%d) != left(%d)\n",
          root->children[0]->height, left);
      assert(false);
    }

    right = assert_invariant(root->children[1]);
    if (root->children[1] == NULL && right != 0) {
      printf("ERROR(3): root->children[1] == NULL && right(%d) != 0\n", right);
      assert(false);
    } else if (root->children[1] && (root->children[1]->height != right)) {
      printf("ERROR(4): root->children[1]->hieght(%d) != right(%d)\n",
          root->children[1]->height, right);
      assert(false);
    }

    if (root == root->children[0] || root == root->children[1]) {
      printf("root == a child\n");
      assert(false);
    }

    for (int dir = 0; dir <= 1; dir++) {
      if (root->children[dir] && root->children[dir]->parent != root) {
        assert(false);
      }
    }

    if ((root->height != left + 1 && root->height != right + 1) ||
        (root->height <= left || root->height <= right)) {
      printf(
        "root height is not correct. heights -- root: %d left: %d right: %d\n",
        root->height, left, right);
      assert(false);
    }

    if (ABS(right - left) > 1) {
      assert(false);
    }

    return 1 + (left > right ? left : right);
  }
}

void fix_height(struct bucket *bucket) {
  struct bucket *left, *right;
  int prior_height;

  while (bucket != NULL) {
    left = bucket->children[0];
    right = bucket->children[1];

    prior_height = bucket->height;
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

    bucket = bucket->parent;
  }
}

struct bucket * rotate_single(
    struct decaying_histogram *histogram, struct bucket *root, int dir) {
  struct bucket *new_root;
  int dir_from_parent;

  new_root = root->children[!dir];
  new_root->parent = root->parent;

  set_child(root, new_root->children[dir], !dir);
  set_child(new_root, root, dir);

  if (new_root->parent != NULL) {
    dir_from_parent = (new_root->parent->children[1] == root);
    new_root->parent->children[dir_from_parent] = new_root;
  } else {
    histogram->root = new_root;
  }

  fix_balance(histogram, new_root->children[0], false);
  fix_balance(histogram, new_root->children[1], false);
  fix_balance(histogram, new_root, false);

  return new_root;
}

struct bucket * rotate_double(
    struct decaying_histogram *histogram, struct bucket *root, int dir) {
  struct bucket *new_root;
  int dir_from_parent;

  new_root = root->children[!dir]->children[dir];
  new_root->parent = root->parent;

  set_child(root->children[!dir], new_root->children[!dir], dir);
  set_child(new_root, root->children[!dir], !dir);
  set_child(root, new_root->children[dir], !dir);
  set_child(new_root, root, dir);

  if (new_root->parent != NULL) {
    if (new_root->parent->children[0] == root)
      dir_from_parent = 0;
    else
      dir_from_parent = 1;
    new_root->parent->children[dir_from_parent] = new_root;
  } else {
    histogram->root = new_root;
  }

  fix_balance(histogram, new_root->children[0], false);
  fix_balance(histogram, new_root->children[1], false);
  fix_balance(histogram, new_root, false);

  return new_root;
}

int compute_balance(struct bucket *bucket) {
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

void fix_balance(
    struct decaying_histogram *histogram, struct bucket *bucket, bool upward) {
  int balance;
  int prior_height;
  int dir;

  fix_height(bucket);
  while (bucket != NULL) {
    prior_height = bucket->height;
    balance = compute_balance(bucket);
    while (ABS(compute_balance(bucket)) >= 2) {
      balance = compute_balance(bucket);
      dir = (balance < 0);

      if (bucket->children[!dir]->children[!dir] == NULL) {
        bucket = rotate_double(histogram, bucket, dir);
      } else if (bucket->children[!dir]->children[dir] == NULL) {
        bucket = rotate_single(histogram, bucket, dir);
      } else if (bucket->children[!dir]->children[dir]->height >
          bucket->children[!dir]->children[!dir]->height) {
        bucket = rotate_double(histogram, bucket, dir);
      } else {
        bucket = rotate_single(histogram, bucket, dir);
      }
    }
    fix_height(bucket);

    if (prior_height == bucket->height)
      bucket = bucket->parent;
    if (upward == false)
      break;
  }
}

void split_bucket(
    struct decaying_histogram *histogram, struct bucket *bucket, int mp_flag) {
  double lower_bound, upper_bound, diameter, median;
  struct bucket *new_bucket, *memo;

  // TODO: This needs to check if it would exceed max_num_buckets. If so, it
  // needs ot first reclaim some old buckets.

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_lock(&histogram->tree_mtx);

  obtain_write_lock(bucket, mp_flag, __LINE__);
  if (bucket->count <= histogram->split_bucket_threshold) {
    // We raced and lost.
    release_write_lock(bucket, mp_flag);
    if (mp_flag & DHIST_MULTI_THREADED)
      pthread_mutex_unlock(&histogram->tree_mtx);
    return;
  }

  new_bucket = init_bucket(histogram, mp_flag);
  ++histogram->num_buckets;
  new_bucket->parent = bucket;

  bucket->count /= 2.0;
  new_bucket->count = bucket->count;
  new_bucket->update_generation = bucket->update_generation;
  new_bucket->below = bucket->below;
  bucket->below = new_bucket;

  if (histogram->num_buckets == 2) {
    // This is an awkward split because we don't have enough information to
    // compute a diameter. Instead, we can set mu to the same value in both
    // buckets, which will also make the boundary between them be mu. As soon
    // as we observe a non-mu entry, the split between the buckets
    // will be more sensical. However, we will need to take care to not split
    // a bucket again until this happens.
    new_bucket->mu = bucket->mu;
  } else {
    lower_bound = compute_bound(histogram, new_bucket->below, bucket);
    upper_bound = compute_bound(histogram, bucket, bucket->above);
    diameter = upper_bound - lower_bound;
    median = lower_bound + diameter / 2.0;
    new_bucket->mu = median - diameter / 6.0;
    bucket->mu = median + diameter / 6.0;
  }
  release_write_lock(bucket, mp_flag);

  new_bucket->above = bucket;
  if (new_bucket->below) {
    obtain_write_lock(new_bucket->below, mp_flag, __LINE__);
    memo = new_bucket->below->above;
    new_bucket->below->above = new_bucket;
    unlock_boundary(new_bucket->below, mp_flag);
    unlock_boundary(memo, mp_flag);
  }

  set_child(new_bucket, bucket->children[0], 0);
  set_child(bucket, new_bucket, 0);

  fix_balance(histogram, new_bucket, true);
  pthread_mutex_unlock(&histogram->tree_mtx);
}

void delete_bucket(
    struct decaying_histogram *histogram, struct bucket *bucket, int mp_flag) {
  struct bucket *cursor, *attach, *lucky_bucket;
  int dir_from_parent, dir, balance;
  uint64_t generation;
  double below_count, above_count;
  pthread_mutex_t *third_boundary_mtx;

 restart:
  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_lock(&histogram->tree_mtx);
  obtain_write_lock(bucket, mp_flag, __LINE__);
  if (bucket->count >= histogram->delete_bucket_threshold ||
      histogram->num_buckets == 1) {
    // We raced and lost.
    release_write_lock(bucket, mp_flag);
    pthread_mutex_unlock(&histogram->tree_mtx);
    return;
  }

  generation = get_generation(histogram, false, mp_flag);

  if (bucket->below == NULL) {
    lucky_bucket = bucket->above;
  } else if (bucket->above == NULL) {
    lucky_bucket = bucket->below;
  } else {
    below_count = compute_count(histogram, bucket->below, generation);
    above_count = compute_count(histogram, bucket->above, generation);
    if (below_count < above_count) {
      lucky_bucket = bucket->below;
    } else {
      lucky_bucket = bucket->above;
    }
  }

  // Now we need to extend our write lock to also cover the lucky bucket.
  if (lucky_bucket == bucket->below)
    third_boundary_mtx = lucky_bucket->boundary_mtx;
  else if (lucky_bucket->above != NULL)
    third_boundary_mtx = lucky_bucket->above->boundary_mtx;
  else
    third_boundary_mtx = NULL;

  if (mp_flag & DHIST_MULTI_THREADED &&
      third_boundary_mtx != NULL &&
      pthread_mutex_trylock(third_boundary_mtx) != 0) {
    // Nuts. It isn't available. Let's restart once it is available.
    release_write_lock(bucket, mp_flag);
    pthread_mutex_unlock(&histogram->tree_mtx);
    pthread_mutex_lock(third_boundary_mtx);
    pthread_mutex_unlock(third_boundary_mtx);
    goto restart;
  }

  lucky_bucket->mu =
      (lucky_bucket->mu * lucky_bucket->count +
       bucket->mu * bucket->count) /
      (lucky_bucket->count + bucket->count);
  lucky_bucket->count += bucket->count;
  if (bucket->below == lucky_bucket) {
    lucky_bucket->above = bucket->above;
    if (lucky_bucket->above)
      lucky_bucket->above->below = lucky_bucket;
  } else {
    lucky_bucket->below = bucket->below;
    if (lucky_bucket->below)
      lucky_bucket->below->above = lucky_bucket;
  }
  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_unlock(third_boundary_mtx);

  if (bucket->children[0] == NULL || bucket->children[1] == NULL) {
    dir_from_parent = (bucket->parent->children[1] == bucket);
    dir = (bucket->children[0] == NULL);
    set_child(bucket->parent, bucket->children[dir], dir_from_parent);
    attach = bucket->parent;
    recycle_bucket(bucket, mp_flag);
  } else {
    balance = compute_balance(bucket);
    dir = (balance <= 0);
    attach = bucket->children[dir];
    while (attach->children[!dir] != NULL)
      attach = attach->children[!dir];
    set_child(attach, bucket->children[!dir], !dir);

    if (bucket->parent == NULL) {
      histogram->root = bucket->children[dir];
      bucket->children[dir]->parent = NULL;
      recycle_bucket(bucket, mp_flag);
      fix_height(histogram->root);
    } else {
      dir_from_parent = (bucket->parent->children[1] == bucket);
      cursor = bucket->children[dir];
      cursor->parent = bucket->parent;
      recycle_bucket(bucket, mp_flag);
      cursor->parent->children[dir_from_parent] = cursor;
      fix_height(cursor->parent);
    }
  }
  fix_balance(histogram, attach, true);

  --histogram->num_buckets;

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_unlock(&histogram->tree_mtx);
}

// This is only safe if histogram->tree_mtx is held.
void obtain_write_lock(struct bucket *bucket, int mp_flag, int line) {
  struct bucket *first, *second, *swap;

  if (mp_flag & DHIST_SINGLE_THREADED)
    return;

  if (bucket->above == NULL) {
    lock_boundary(bucket, mp_flag, line);
    return;
  }

  first = bucket;
  second = bucket->above;
  while (1) {
    lock_boundary(first, mp_flag, line);
    if (trylock_boundary_succeeded(second, mp_flag, __LINE__)) {
      return;
    }
    unlock_boundary(first, mp_flag);
    swap = first;
    first = second;
    second = swap;
  }
}

void release_write_lock(struct bucket *bucket, int mp_flag) {
  if (mp_flag & DHIST_SINGLE_THREADED)
    return;
  unlock_boundary(bucket, mp_flag);
  if (bucket->above)
    unlock_boundary(bucket->above, mp_flag);
}

void _print_tree(struct bucket *bucket, int depth) {
  printf("%.2lf", bucket->mu);

  if (bucket->children[0]) {
    printf("\t  ");
    _print_tree(bucket->children[0], depth + 1);
  }

  if (bucket->children[1]) {
    printf("\n");
    for (int i = 0; i < depth; i++)
      printf("\t");
    printf("\t\\ ");
    _print_tree(bucket->children[1], depth + 1);
  }
}

void print_tree(struct bucket *bucket) {
  _print_tree(bucket, 0);
  printf("\n");
}

void assert_consistent(struct decaying_histogram *histogram) {
  double upper_bound, lower_bound;
  struct bucket *cursor = histogram->root;
  while (cursor->children[0])
    cursor = cursor->children[0];

  while (cursor != NULL) {
    upper_bound = compute_bound(histogram, cursor, cursor->above);
    lower_bound = compute_bound(histogram, cursor->below, cursor);
    if (cursor->above && cursor->mu > cursor->above->mu) {
      printf("ERROR: cursor->mu(%lf) > cursor->above->mu(%lf)\n",
          cursor->mu, cursor->above->mu);
      print_tree(histogram->root);
      assert(0);
    }
    if (upper_bound < lower_bound) {
      printf("ERROR: upper_bound(%lf) < lower_bound(%lf)\n",
          upper_bound, lower_bound);
      print_tree(histogram->root);
      assert(0);
    }
    cursor = cursor->above;
  }

  assert_invariant(histogram->root);
}

