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

#include "dhist.h"

// The same as ceil(x / y). Using this so that math.h is not a dependency.
#ifndef CEIL
#define CEIL(x, y) ((int)((x) + ((int)(y) - 1)) / (int)(y))
#endif

#ifndef ABS
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#endif

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif


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
  // The is_enabled field exists because the bucket lookup is not protected
  // with any atomics, so it is possible for a thread to obtain a reference to
  // this bucket while it is enabled, for another thread to recycle the bucket,
  // and then for the first thread to attempt to update the bucket.
  bool is_enabled;
  bool lock_held;
  char __padding[2];
};

struct dhist_info {
  struct dhist *histogram;
  double *weights;
  double *boundaries;
  uint64_t generation;
  int num_buckets;
  int num_boundaries;
};

static double ipow(double coefficient, uint64_t power);

static struct bucket * init_bucket(struct dhist *histogram, int mp_flag);
static bool is_target_boundary(struct bucket *bucket, double observation);
static int compute_balance(struct bucket *bucket);
static void obtain_write_lock(struct bucket *bucket, int mp_flag);
static void release_write_lock(struct bucket *bucket, int mp_flag);
static void lock_boundary(struct bucket *bucket, int mp_flag);
static bool trylock_boundary_succeeded(struct bucket *bucket, int mp_flag);
static void unlock_boundary(struct bucket *bucket, int mp_flag);
static void _print_tree(struct bucket *bucket, int depth);
static void print_tree(struct bucket *bucket);
static void fix_height(struct bucket *bucket);
static void set_child(struct bucket *root, struct bucket *child, int dir);
static int count_buckets_in_tree(struct bucket *root);
static int assert_invariant(struct bucket *root);

static double get_decay(
    struct dhist *histogram, uint64_t missed_generations);
static uint64_t get_generation(
    struct dhist *histogram, bool increment, int mp_flag);
static double compute_bound(
    struct dhist *histogram, struct bucket *left, struct bucket *right);
static double compute_count(
    struct dhist *histogram, struct bucket *bucket, uint64_t generation);
static double compute_total_count(
    struct dhist *histogram, uint64_t generation);
static bool perform_add(
    struct dhist *histogram, struct bucket *bucket,
    double observation, bool recheck_left_boundary,
    bool recheck_right_boundary, int mp_flag);
static struct bucket * find_bucket(
    struct dhist *histogram, double observation, int mp_flag);
static void fix_balance(
    struct dhist *histogram, struct bucket *bucket, bool upward);
static struct bucket * rotate_single(
    struct dhist *histogram, struct bucket *root, int dir);
static struct bucket * rotate_double(
    struct dhist *histogram, struct bucket *root, int dir);
static void split_bucket(
    struct dhist *histogram, struct bucket *bucket, int mp_flag);
static void _split_bucket(
    struct dhist *histogram, struct bucket *bucket, int mp_flag);
static void delete_bucket(
    struct dhist *histogram, struct bucket *bucket, int mp_flag);
static void _delete_bucket(
    struct dhist *histogram, struct bucket *bucket, int mp_flag);
static void decay(
    struct dhist *histogram, struct bucket *bucket, uint64_t generation);
static int snprint_histogram(
    char *s_buffer, size_t n, struct dhist_info *info, const char *title,
    const char *xlabel);
static void assert_consistent(struct dhist *histogram);

static void extract_info(
    struct dhist *histogram, int mp_flag, struct dhist_info *info);
static void union_of_boundaries(
    struct dhist_info *info1, struct dhist_info *info2,
    struct dhist_info *info_union);
static void redistribute(
    struct dhist_info *info_orig, struct dhist_info *info_union,
    struct dhist_info *info_redist);
static void clean_info(struct dhist_info *info);


const int DHIST_SINGLE_THREADED = (1 << 0);
const int DHIST_MULTI_THREADED = (1 << 1);


struct dhist *
dhist_init(int target_buckets, double alpha) {
  struct dhist *histogram;
  double max_count, expected_count;
  uint64_t idx;
  struct bucket *bucket;
  double radius;

  histogram = (struct dhist *)malloc(sizeof(struct dhist));

  max_count = 1.0 / alpha;
  expected_count = 1.0 / (alpha * target_buckets);

  // XXX(LPE): What are good thresholds?
  radius = 0.8;
  histogram->delete_bucket_threshold = expected_count * (1.0 - radius);
  histogram->split_bucket_threshold = expected_count * (1.0 + radius);
  histogram->alpha = alpha;
  histogram->max_num_buckets =
      (uint32_t)CEIL(max_count, histogram->delete_bucket_threshold);
  histogram->bucket_list = (struct bucket *)malloc(
      histogram->max_num_buckets * sizeof(struct bucket));
  for (idx = 0; idx < histogram->max_num_buckets; idx++) {
    bucket = &histogram->bucket_list[idx];
    bucket->is_enabled = false;
    bucket->lock_held = false;
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

  histogram->tree_mtx = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));
  pthread_mutex_init(histogram->tree_mtx, NULL);

  histogram->generation_mtx =
      (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));
  pthread_mutex_init(histogram->generation_mtx, NULL);

  return histogram;
}

void dhist_destroy(struct dhist *histogram) {
  // TODO free the buckets.
  free(histogram->pow_table);
  free(histogram->tree_mtx);
  free(histogram->generation_mtx);
  free(histogram);
}

void dhist_insert(
    struct dhist *histogram, double observation, int mp_flag) {
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

    if (!bucket->below) {
      // Always insert the observaiton into this bucket.
      boundary = observation;
    } else {
      boundary = compute_bound(histogram, bucket->below, bucket);
    }

    if (mp_flag & DHIST_SINGLE_THREADED) {
      // Single threaded, handles both cases.
      add_succeeded = perform_add(
          histogram, (boundary <= observation) ? bucket : bucket->below,
          observation, false, false, mp_flag);
    } else if (observation < boundary) {
      // We need to insert into bucket->below. That bucket exists because we
      // would have created an artificial boundary if it didn't.
      if (trylock_boundary_succeeded(bucket->below, mp_flag)) {
        add_succeeded = perform_add(
            histogram, bucket->below, observation, true, false, mp_flag);
        unlock_boundary(bucket->below, mp_flag);
      } else {
        add_succeeded = false;
      }
    } else if (bucket->above) {
      // We need to insert into this bucket, and since we have a neighbor to
      // the right, we need to extend our lock.
      lock_boundary(bucket->above, mp_flag);
      add_succeeded = perform_add(
          histogram, bucket, observation, false, true, mp_flag);
      unlock_boundary(bucket->above, mp_flag);
    } else {
      // No bucket exists to the right, and since we already have the boundary
      // mutex for this bucket, nothing else can add a bucket to the right.

      // XXX split and delete_bucket should lock far left and far right buckets
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

// The caller needs to call free on the result.
char * dhist_get_json(
    struct dhist *histogram, const char *title, const char *xlabel,
    int mp_flag) {
  int num_chars;
  struct dhist_info info;
  char *buffer;

  extract_info(histogram, mp_flag, &info);
  num_chars = snprint_histogram(NULL, 0, &info, title, xlabel);
  buffer = (char *)malloc(sizeof(char) * (size_t)(num_chars + 1));
  snprint_histogram(
      buffer, (size_t)(num_chars + 1), &info, title, xlabel);
  clean_info(&info);

  return buffer;
}

void dhist_set_target_buckets(struct dhist *histogram, int target_buckets) {
  assert(false);
}

void dhist_set_alpha(struct dhist *histogram, double alpha) {
  assert(false);
}

double dhist_Jaccard_distance(
    struct dhist *hist1, struct dhist *hist2, int mp_flag) {
  int idx;
  struct dhist_info info1, info2, info_union, info1_redist, info2_redist;
  double distance, union_area, intersection_area, width;

  extract_info(hist1, mp_flag, &info1);
  extract_info(hist2, mp_flag, &info2);
  union_of_boundaries(&info1, &info2, &info_union);
  redistribute(&info1, &info_union, &info1_redist);
  redistribute(&info2, &info_union, &info2_redist);

  union_area = intersection_area = 0.0;
  for (idx = 0; idx < info_union.num_buckets; idx++) {
    width = info_union.boundaries[idx + 1] - info_union.boundaries[idx];
    union_area +=
        width * MAX(info1_redist.weights[idx], info2_redist.weights[idx]);
    intersection_area +=
        width * MIN(info1_redist.weights[idx], info2_redist.weights[idx]);
  }
  distance = 1.0 - intersection_area / union_area;

  clean_info(&info1);
  clean_info(&info2);
  clean_info(&info_union);
  clean_info(&info1_redist);
  clean_info(&info2_redist);

  return distance;
}

double dhist_Kolmogorov_Smirnov_statistic(
    struct dhist *hist1, struct dhist *hist2, int mp_flag) {
  int idx;
  struct dhist_info info1, info2, info_union, info1_redist, info2_redist;
  double distance, width, cdf[2];

  extract_info(hist1, mp_flag, &info1);
  extract_info(hist2, mp_flag, &info2);
  union_of_boundaries(&info1, &info2, &info_union);
  redistribute(&info1, &info_union, &info1_redist);
  redistribute(&info2, &info_union, &info2_redist);

  distance = cdf[0] = cdf[1] = 0.0;
  for (idx = 0; idx < info_union.num_buckets; idx++) {
    width = info_union.boundaries[idx + 1] - info_union.boundaries[idx];
    cdf[0] += width * info1_redist.weights[idx];
    cdf[1] += width * info2_redist.weights[idx];
    if (ABS(cdf[1] - cdf[0]) > distance)
      distance = ABS(cdf[1] - cdf[0]);
  }

  clean_info(&info1);
  clean_info(&info2);
  clean_info(&info_union);
  clean_info(&info1_redist);
  clean_info(&info2_redist);

  return distance;
}

double dhist_earth_movers_distance(
    struct dhist *hist1, struct dhist *hist2, int mp_flag) {
  int idx;
  struct dhist_info info1, info2, info_union, info1_redist, info2_redist;
  double distance, burden, step;

  extract_info(hist1, mp_flag, &info1);
  extract_info(hist2, mp_flag, &info2);
  union_of_boundaries(&info1, &info2, &info_union);
  redistribute(&info1, &info_union, &info1_redist);
  redistribute(&info2, &info_union, &info2_redist);

  distance = burden = step = 0.0;
  for (idx = 0; idx < info_union.num_buckets - 1; idx++) {
    burden += info2_redist.weights[idx] - info1_redist.weights[idx];
    // step is distance between above bucket's mean and this bucket's mean.
    step = (
        (info_union.boundaries[idx + 2] +
         info_union.boundaries[idx + 1]) / 2.0 -
        (info_union.boundaries[idx + 1] +
         info_union.boundaries[idx]) / 2.0);
    distance += ABS(burden) * step;
  }

  clean_info(&info1);
  clean_info(&info2);
  clean_info(&info_union);
  clean_info(&info1_redist);
  clean_info(&info2_redist);

  return distance;
}

static double
ipow(double coefficient, uint64_t power) {
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

static struct bucket *
init_bucket(
    struct dhist *histogram, int mp_flag) {
  struct bucket *bucket = NULL;
  uint32_t idx;

  for (idx = 0; idx < histogram->max_num_buckets; idx++) {
    bucket = &histogram->bucket_list[idx];
    if (bucket->is_enabled)
      continue;
    if (!trylock_boundary_succeeded(bucket, mp_flag))
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
  bucket->height = 1;
  unlock_boundary(bucket, mp_flag);

  return bucket;
}

/*
 * Returns true if the observation is between the lower bucket mu and the
 * current bucket mu.
 * obs \in [bucket->below->mu, bucket->mu)
 */
static bool
is_target_boundary(struct bucket *bucket, double observation) {
  double left_mu, mu;
  struct bucket *below;

  if ((below = bucket->below) != NULL) {
    left_mu = below->mu;
  } else {
    left_mu = observation - 1.0;
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

static int
compute_balance(struct bucket *bucket) {
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

// This is only safe if histogram->tree_mtx is held.
static void
obtain_write_lock(struct bucket *bucket, int mp_flag) {
  struct bucket *first, *second, *swap;

  if (mp_flag & DHIST_SINGLE_THREADED)
    return;

  if (bucket->above == NULL) {
    lock_boundary(bucket, mp_flag);
    return;
  }

  first = bucket;
  second = bucket->above;
  while (1) {
    lock_boundary(first, mp_flag);
    if (trylock_boundary_succeeded(second, mp_flag)) {
      return;
    }
    unlock_boundary(first, mp_flag);
    swap = first;
    first = second;
    second = swap;
  }
}

static void
release_write_lock(struct bucket *bucket, int mp_flag) {
  if (mp_flag & DHIST_SINGLE_THREADED)
    return;
  unlock_boundary(bucket, mp_flag);
  if (bucket->above)
    unlock_boundary(bucket->above, mp_flag);
}

static void
lock_boundary(struct bucket *bucket, int mp_flag) {
  if (mp_flag & DHIST_MULTI_THREADED) {
    pthread_mutex_lock(bucket->boundary_mtx);
    bucket->lock_held = true;
  }
}

static bool
trylock_boundary_succeeded(struct bucket *bucket, int mp_flag) {
  if (mp_flag & DHIST_SINGLE_THREADED) {
    return true;
  } else if (pthread_mutex_trylock(bucket->boundary_mtx) == 0) {
    bucket->lock_held = true;
    return true;
  } else {
    return false;
  }
}

static void
unlock_boundary(struct bucket *bucket, int mp_flag) {
  if (mp_flag & DHIST_MULTI_THREADED) {
    bucket->lock_held = false;
    pthread_mutex_unlock(bucket->boundary_mtx);
  }
}

static void
_print_tree(struct bucket *bucket, int depth) {
  int i;

  printf("%.2lf", bucket->mu);

  if (bucket->children[0]) {
    printf("\t  ");
    _print_tree(bucket->children[0], depth + 1);
  }

  if (bucket->children[1]) {
    printf("\n");
    for (i = 0; i < depth; i++)
      printf("\t");
    printf("\t\\ ");
    _print_tree(bucket->children[1], depth + 1);
  }
}

static void
print_tree(struct bucket *bucket) {
  _print_tree(bucket, 0);
  printf("\n");
}

static void
fix_height(struct bucket *bucket) {
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

static void
set_child(struct bucket *root, struct bucket *child, int dir) {
  root->children[dir] = child;
  if (child != NULL)
    child->parent = root;
  fix_height(root);
}

static int
count_buckets_in_tree(struct bucket *root) {
  if (root) {
    return 1 +
        count_buckets_in_tree(root->children[0]) +
        count_buckets_in_tree(root->children[1]);
  } else {
    return 0;
  }
}

static int
assert_invariant(struct bucket *root) {
  int left, right, dir;

  if (root == NULL) {
    return 0;
  } else if (root->children[0] == NULL && root->children[1] == NULL) {
    if (root->height != 1) {
      print_tree(root);
      assert(root->height == 1);
    }
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

    for (dir = 0; dir <= 1; dir++) {
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

static double
get_decay(
    struct dhist *histogram, uint64_t missed_generations) {
  if (missed_generations < (uint32_t)histogram->max_num_buckets)
    return histogram->pow_table[missed_generations];
  else
    return ipow(1.0 - histogram->alpha, missed_generations);
}

static uint64_t
get_generation(
    struct dhist *histogram, bool increment, int mp_flag) {
  uint64_t generation;

  if (mp_flag & DHIST_MULTI_THREADED) {
    pthread_mutex_lock(histogram->generation_mtx);
    if (increment)
      ++histogram->generation;
    generation = histogram->generation;
    pthread_mutex_unlock(histogram->generation_mtx);
  } else {
    if (increment)
      ++histogram->generation;
    generation = histogram->generation;
  }

  return generation;
}

// We can't increment generation until we finally make our update, but
// since the boundaries don't move upon decay, we can use any generation
// value to compute the bounds. We'll use the largest generation so that
// we don't need to recompute the decay for that bucket.
static double
compute_bound(
    struct dhist *histogram, struct bucket *left, struct bucket *right) {
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

static double compute_count(
    struct dhist *histogram, struct bucket *bucket, uint64_t generation) {
  return (
      bucket->count *
      get_decay(histogram, generation - bucket->update_generation));
}

static double compute_total_count(struct dhist *histogram, uint64_t generation) {
  return (1 - get_decay(histogram, generation)) /
         histogram->alpha;
}

static bool
perform_add(
    struct dhist *histogram, struct bucket *bucket, double observation,
    bool recheck_left_boundary, bool recheck_right_boundary, int mp_flag) {
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

/*
 * This returns the bucket with the greatest mu less than observation (which is
 * not necessarily the bucket where observation will be inserted) unless the
 * target bucket is the leftmost bucket, in which case observation may be less
 * than bucket->mu.
 * If mp_flag & DHIST_MULTI_THREADED, the boundary_mtx for this bucket and its
 * right hand neighbor will be locked.
 */
static struct bucket * find_bucket(
    struct dhist *histogram, double observation, int mp_flag) {
  struct bucket *bucket, *other;

  do {
    bucket = histogram->root;
    do {
      // The first two branches are quick decisions that avoid computing
      // the region where the boundary can exist.
      if ((other = bucket->below) && observation < other->mu) {
        bucket = bucket->children[0];
      } else if ((other = bucket->above) && other->mu <= observation) {
        bucket = bucket->children[1];
      } else if (is_target_boundary(bucket, observation)) {
        // It appears that we've found the target boundary, but we haven't
        // locked anything, so lock the bucket and make sure it's still the
        // boundary we wanted.
        if (mp_flag & DHIST_SINGLE_THREADED) {
          // There is nothing to race against.
          return bucket;
        } else {
          lock_boundary(bucket, mp_flag);

          if (!is_target_boundary(bucket, observation) || !bucket->is_enabled) {
            // We raced, so restart.
            unlock_boundary(bucket, mp_flag);
            bucket = histogram->root;
          } else {
            return bucket;
          }
        }
      } else {
        bucket = bucket->children[bucket->mu <= observation];
      }
    } while (bucket != NULL);
    // We can fail to find a bucket if the bucket boundaries moved around
    // during the search. Start over.
  } while (1);
  return NULL;
}

static void
fix_balance(
    struct dhist *histogram, struct bucket *bucket, bool upward) {
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

static struct bucket *
rotate_single(struct dhist *histogram, struct bucket *root, int dir) {
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

static struct bucket *
rotate_double(struct dhist *histogram, struct bucket *root, int dir) {
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

static void
split_bucket(struct dhist *histogram, struct bucket *bucket, int mp_flag) {
  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_lock(histogram->tree_mtx);

  if (bucket->is_enabled)
    _split_bucket(histogram, bucket, mp_flag);

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_unlock(histogram->tree_mtx);
}

static void
_split_bucket(struct dhist *histogram, struct bucket *bucket, int mp_flag) {
  double lower_bound, upper_bound, diameter, median;
  struct bucket *new_bucket, *cursor, *memo;

  // If we would exceed max_num_buckets, we need to first reclaim some old
  // buckets.
  if (histogram->num_buckets >= histogram->max_num_buckets) {
    cursor = histogram->root;
    while (cursor->children[0])
      cursor = cursor->children[0];

    while (cursor) {
      memo = cursor->above;
      if (cursor->count < histogram->delete_bucket_threshold &&
          histogram->num_buckets > 1) {
        _delete_bucket(histogram, cursor, mp_flag);
      }
      cursor = memo;
    }
  }

  obtain_write_lock(bucket, mp_flag);
  if (bucket->count <= histogram->split_bucket_threshold) {
    // We raced and lost.
    release_write_lock(bucket, mp_flag);
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
    upper_bound = compute_bound(histogram, bucket, bucket->above);
    if (!new_bucket->below) {
      lower_bound = bucket->mu - (upper_bound - bucket->mu);
    } else {
      lower_bound = compute_bound(histogram, new_bucket->below, bucket);
    }
    diameter = upper_bound - lower_bound;
    median = lower_bound + diameter / 2.0;
    new_bucket->mu = median - diameter / 6.0;
    bucket->mu = median + diameter / 6.0;
    //printf("??? lower_bound: %lf, upper_bound: %lf, "
    //       "diameter: %lf, median: %lf, new_bucket->mu: %lf, "
    //       "bucket->mu: %lf\n",
    //       lower_bound, upper_bound,
    //       diameter, median, new_bucket->mu,
    //       bucket->mu);
  }
  release_write_lock(bucket, mp_flag);

  new_bucket->above = bucket;
  if (new_bucket->below) {
    obtain_write_lock(new_bucket->below, mp_flag);
    memo = new_bucket->below->above;
    new_bucket->below->above = new_bucket;
    unlock_boundary(new_bucket->below, mp_flag);
    unlock_boundary(memo, mp_flag);
  }

  set_child(new_bucket, bucket->children[0], 0);
  set_child(bucket, new_bucket, 0);

  fix_balance(histogram, new_bucket, true);
}

static void
delete_bucket(struct dhist *histogram, struct bucket *bucket, int mp_flag) {
  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_lock(histogram->tree_mtx);

  if (bucket->is_enabled)
    _delete_bucket(histogram, bucket, mp_flag);

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_unlock(histogram->tree_mtx);
}

static void
_delete_bucket(struct dhist *histogram, struct bucket *bucket, int mp_flag) {
  struct bucket *attach, *lucky_bucket, *other;
  struct bucket *far_buckets[2];
  int dir_from_parent, dir, balance, idx;
  uint64_t generation;
  double below_count, above_count;

 restart:
  obtain_write_lock(bucket, mp_flag);
  if (bucket->count >= histogram->delete_bucket_threshold ||
      histogram->num_buckets == 1) {
    // We raced and lost.
    release_write_lock(bucket, mp_flag);
    return;
  }

  generation = get_generation(histogram, false, mp_flag);

  // Decide whether to merge the dying bucket above or below.
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

  // Now we need to extend our write lock to also cover the two adjacent
  // buckets.
  far_buckets[0] = bucket->below;
  if (bucket->above)
    far_buckets[1] = bucket->above->above;
  else
    far_buckets[1] = NULL;

  if (mp_flag & DHIST_MULTI_THREADED) {
    for (idx = 0; idx <= 1; idx++) {
      if (far_buckets[idx] != NULL &&
          !trylock_boundary_succeeded(far_buckets[idx], mp_flag)) {
        // Nuts. It isn't available. Let's restart once it is available.
        if (idx == 1 && far_buckets[0])
          unlock_boundary(far_buckets[0], mp_flag);
        release_write_lock(bucket, mp_flag);
        lock_boundary(far_buckets[idx], mp_flag);
        unlock_boundary(far_buckets[idx], mp_flag);
        goto restart;
      }
    }
  }

  // Update the lucky bucket data.
  lucky_bucket->mu =
      (lucky_bucket->mu * lucky_bucket->count +
       bucket->mu * bucket->count) /
      (lucky_bucket->count + bucket->count);
  lucky_bucket->count += bucket->count;

  // Remove bucket from the linked list.
  if (bucket->below)
    bucket->below->above = bucket->above;
  if (bucket->above)
    bucket->above->below = bucket->below;

  // We no longer need the outside boundaries.
  if (mp_flag & DHIST_MULTI_THREADED) {
    for (idx = 0; idx <= 1; idx++) {
      if (far_buckets[idx])
        unlock_boundary(far_buckets[idx], mp_flag);
    }
  }

  // Remove bucket from the tree.

  // If we move something to the attach point, we might need to rebalance.
  attach = NULL;
  if (bucket->children[0] == NULL || bucket->children[1] == NULL) {
    // Promote this bucket's child (no more than one exists) to the bucket's
    // spot.
    dir = (bucket->children[0] == NULL);
    if (bucket->parent == NULL) {
      // We're dropping from two buckets back down to one bucket. Reset the
      // root.
      histogram->root = bucket->children[dir];
      bucket->children[dir]->parent = NULL;
      fix_height(histogram->root);
    } else {
      dir_from_parent = (bucket->parent->children[1] == bucket);
      set_child(bucket->parent, bucket->children[dir], dir_from_parent);
      fix_height(bucket->parent);
      attach = bucket->parent;
    }
  } else {
    // We need to promote one of the two children and then attach the other
    // child. After that, we'll need to rebalance.
    balance = compute_balance(bucket);
    dir = (balance <= 0);
    attach = bucket->children[dir];
    while (attach->children[!dir] != NULL)
      attach = attach->children[!dir];
    set_child(attach, bucket->children[!dir], !dir);

    if (bucket->parent == NULL) {
      // We need a new root.
      histogram->root = bucket->children[dir];
      histogram->root->parent = NULL;
      fix_height(histogram->root);
    } else {
      // Move one of the two children up.
      dir_from_parent = (bucket->parent->children[1] == bucket);
      set_child(bucket->parent, bucket->children[dir], dir_from_parent);
    }
  }

  if (attach)
    fix_balance(histogram, attach, true);

  // Return bucket to the pool.
  other = bucket->above;
  bucket->below = bucket->above = NULL;
  bucket->children[0] = bucket->children[1] = NULL;
  bucket->is_enabled = false;

  if (mp_flag & DHIST_MULTI_THREADED) {
    unlock_boundary(bucket, mp_flag);
    if (other)
      unlock_boundary(other, mp_flag);
  }

  --histogram->num_buckets;
}

/*
 * This can efficiently apply the decay for multiple generations. This function
 * will not read histogram->generation since that may be updated by another
 * thread.
 */
static void
decay(struct dhist *histogram, struct bucket *bucket, uint64_t generation) {
  if (bucket == NULL || bucket->update_generation == generation)
    return;

  bucket->count *= get_decay(
      histogram, generation - bucket->update_generation);
  bucket->update_generation = generation;
}

static void extract_info(
    struct dhist *histogram, int mp_flag, struct dhist_info *info) {
  struct bucket *bucket;
  double total_count;
  int idx;

  info->histogram = histogram;

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_lock(histogram->tree_mtx);

  // Allocate enough space for both the weights and the boundaries.
  info->weights = (double *)malloc(
      sizeof(double) * 2 * histogram->max_num_buckets + 1);
  info->boundaries = info->weights + histogram->max_num_buckets;

  info->generation = get_generation(histogram, false, mp_flag);

  bucket = histogram->root;
  while (bucket->below)
    bucket = bucket->below;

  // First, populate weights with the raw counts and identify bucket
  // boundaries. Then compute the total count of all observations.
  // Finally, convert the raw counts into the weights.

  idx = 0;
  lock_boundary(bucket, mp_flag);
  info->boundaries[idx] = compute_bound(histogram, NULL, bucket);
  while (bucket) {
    if (bucket->above)
      lock_boundary(bucket->above, mp_flag);
    info->boundaries[idx + 1] = compute_bound(
        histogram, bucket, bucket->above);
    if (bucket->update_generation < info->generation)
      decay(histogram, bucket, info->generation);
    info->weights[idx] = bucket->count;

    unlock_boundary(bucket, mp_flag);
    bucket = bucket->above;
    idx++;
  }

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_unlock(histogram->tree_mtx);

  info->num_buckets = idx;
  info->num_boundaries = info->num_buckets + 1;

  total_count = 0.0;
  for (idx = 0; idx < info->num_buckets; idx++)
    total_count += info->weights[idx];

  for (idx = 0; idx < info->num_buckets; idx++) {
    info->weights[idx] /=
        (info->boundaries[idx + 1] - info->boundaries[idx]) * total_count;
  }
}

static void
clean_info(struct dhist_info *info) {
  free(info->weights);
  info->weights = NULL;
  info->boundaries = NULL;
  info->num_buckets = 0;
  info->num_boundaries = 0;
}

// Call this with n == 0 to get a count of how many characters (not including
// the \0) needed to print the histogram. Call again with an allocated buffer
// to construct the string.
static int
snprint_histogram(
    char *s_buffer, size_t n, struct dhist_info *info, const char *title,
    const char *xlabel) {
  int idx, num_chars = 0;

  num_chars += snprintf(s_buffer + num_chars, n, "{");
  if (title) {
    num_chars +=
        snprintf(s_buffer + num_chars, n, "\"title\": \"%s\", ", title);
  }
  if (xlabel) {
    num_chars +=
        snprintf(s_buffer + num_chars, n, "\"xlabel\": \"%s\", ", xlabel);
  }

  num_chars += snprintf(
      s_buffer + num_chars, n, "\"generation\": %lu, ", info->generation);
  num_chars += snprintf(
      s_buffer + num_chars, n, "\"id\": %lu, ", (uint64_t)info->histogram);

  num_chars += snprintf(s_buffer + num_chars, n, "\"weights\": [");
  for (idx = 0; idx < info->num_buckets - 1; idx++) {
    num_chars += snprintf(
        s_buffer + num_chars, n, "%lf, ", info->weights[idx]);
  }
  num_chars += snprintf(s_buffer + num_chars, n, "%lf], ", info->weights[idx]);

  num_chars += snprintf(s_buffer + num_chars, n, "\"boundaries\": [");
  for (idx = 0; idx < info->num_buckets; idx++) {
    num_chars += snprintf(
        s_buffer + num_chars, n, "%lf, ", info->boundaries[idx]);
  }
  num_chars += snprintf(
      s_buffer + num_chars, n, "%lf]}", info->boundaries[idx]);

  return num_chars;
}

static void
assert_consistent(struct dhist *histogram) {
  uint32_t num_buckets_seen = 0;
  double upper_bound, lower_bound;
  struct bucket *cursor = histogram->root, *cursor_two;
  while (cursor->children[0])
    cursor = cursor->children[0];

  while (cursor != NULL) {
    num_buckets_seen++;
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

    if (cursor->above) {
      assert(cursor->above->below == cursor);
    }

    if (cursor->below) {
      assert(cursor->below->above == cursor);
    }

    if (cursor->children[0]) {
      assert(cursor->children[0]->parent == cursor);
      assert(cursor->below);
      assert(cursor->below->mu <= cursor->mu);
    }

    if (cursor->children[1]) {
      assert(cursor->children[1]->parent == cursor);
      assert(cursor->above);
      assert(cursor->mu <= cursor->above->mu);
    }

    // Make sure we can find this bucket in the tree.
    cursor_two = histogram->root;
    while (true) {
      if (cursor_two == NULL) {
        printf("Could not find bucket with mu %lf in the tree.\n", cursor->mu);
        assert(false);
      } else if (cursor_two->mu == cursor->mu) {
        break;
      } else if (cursor->mu < cursor_two->mu) {
        cursor_two = cursor_two->children[0];
      } else {
        cursor_two = cursor_two->children[1];
      }
    }

    cursor = cursor->above;
  }

  assert(num_buckets_seen == histogram->num_buckets);
  assert(
      (uint64_t)count_buckets_in_tree(histogram->root) ==
      histogram->num_buckets);
  assert_invariant(histogram->root);
}

/*
 * Note: This does NOT initialize the weights, although it does reserve
 * enough space that the weights can be calculated later.
 */
static void
union_of_boundaries(
    struct dhist_info *info1, struct dhist_info *info2,
    struct dhist_info *union_info) {
  int idx1, idx2, union_idx;

  union_info->histogram = NULL;
  union_info->generation = 0;
  union_info->weights = (double *)malloc(sizeof(double) * (size_t)(
      info1->num_buckets + info1->num_boundaries +
      info2->num_buckets + info2->num_boundaries));
  union_info->boundaries =
      union_info->weights + (info1->num_buckets + info2->num_buckets);

  idx1 = idx2 = union_idx = 0;
  while (idx1 < info1->num_boundaries || idx2 < info2->num_boundaries) {
    if (idx1 < info1->num_boundaries && idx2 < info2->num_boundaries &&
        info1->boundaries[idx1] == info2->boundaries[idx2]) {
      union_info->boundaries[union_idx] = info2->boundaries[idx1];
      idx1++;
      idx2++;
    } else if (idx2 == info2->num_boundaries ||
               ((idx1 < info1->num_boundaries) &&
                (info1->boundaries[idx1] < info2->boundaries[idx2]))) {
      union_info->boundaries[union_idx] = info1->boundaries[idx1];
      idx1++;
    } else {
      union_info->boundaries[union_idx] = info2->boundaries[idx2];
      idx2++;
    }
    union_idx++;
  }
  union_info->num_boundaries = union_idx;
  union_info->num_buckets = union_idx - 1;
}

static void
redistribute(
    struct dhist_info *info_orig, struct dhist_info *info_union,
    struct dhist_info *info_redist) {
  int orig_idx, redist_idx;

  info_redist->histogram = NULL;
  info_redist->weights = (double *)malloc(sizeof(double) *
      (size_t)(info_union->num_buckets + info_union->num_boundaries));
  info_redist->boundaries = info_redist->weights + info_union->num_buckets;
  info_redist->num_buckets = info_union->num_buckets;
  info_redist->num_boundaries = info_union->num_boundaries;

  orig_idx = 0;
  for (redist_idx = 0; redist_idx < info_union->num_buckets; redist_idx++) {
    info_redist->boundaries[redist_idx] = info_union->boundaries[redist_idx];
    if (info_orig->boundaries[orig_idx + 1] ==
        info_union->boundaries[redist_idx]) {
      // Shift the orig boundary that we're considering.
      orig_idx++;
    }

    if (orig_idx == info_orig->num_buckets) {
      // We've exhausted the distribution. Fill in the rest and exit.
      for (; redist_idx < info_union->num_buckets; redist_idx++) {
        info_redist->boundaries[redist_idx] =
            info_union->boundaries[redist_idx];
        info_redist->weights[redist_idx] = 0.0;
      }
      break;
    }

    if (info_union->boundaries[redist_idx] < info_orig->boundaries[0]) {
      info_redist->weights[redist_idx] = 0.0;
    } else {
      info_redist->weights[redist_idx] = info_orig->weights[orig_idx];
    }
  }
  info_redist->boundaries[info_redist->num_boundaries - 1] =
      info_union->boundaries[info_redist->num_boundaries - 1];
}

