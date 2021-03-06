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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "dhist_fixedptc.h"

#ifndef ABS
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#endif

#define NUM_PRECOMPUTED_POWERS 200


// Interestingly, splitting bucket_data out from bucket provides a 4% speedup.
// This is likely due to the data that rarely changes (pointers only change on
// a split_bucket or a delete_bucket event) from the data the changes
// frequently. This means that the cache lines that hold the pointers are
// invalidated far less frequently than the cache lines for high-churn data.
struct bucket_data {
  // The values for count, mu, update_generation and the below and above
  // pointers are protected by bucket->boundary_mtx and
  // bucket->above->boundary_mtx. If this thread holds either of those
  // mutexes, these values cannot be modified by any other thread.
  fixedpt count;
  fixedpt mu;
  uint64_t update_generation;
  int height;
  // The is_enabled field exists because the bucket lookup is not protected
  // with any atomics, so it is possible for a thread to obtain a reference to
  // this bucket while it is enabled, for another thread to recycle the bucket,
  // and then for the first thread to attempt to update the bucket.
  bool is_enabled;
  bool lock_held;
  bool propogate_rebalance;
  bool __padding;
  struct bucket *fix_balance_stack_next;
};

struct bucket {
  struct bucket_data *data;
  // The parent and children pointers, height, and is_enabled fields are
  // protected by tree_mtx.
  struct bucket *children[2];
  struct bucket *parent;
  struct bucket *below;
  struct bucket *above;
  pthread_mutex_t *boundary_mtx;  /* lower boundary */
};

struct dhist_info {
  struct dhist *histogram;
  union {
    fixedpt *weights;
    fixedpt *CDF;
  };
  fixedpt *boundaries;
  uint64_t generation;
  int num_buckets;
  int num_boundaries;
};

struct thread_info {
  struct thread_info *prev;
  struct thread_info *next;
  struct bucket *bucket_to_free;
};

static fixedpt ipow(fixedpt coefficient, uint64_t power);

static struct bucket * bucket_init(int mp_flag);
static void bucket_destroy(struct bucket *bucket);

void thread_info_init_fields(
    struct dhist *histogram, struct thread_info *info, int mp_flag);
static void thread_info_finalize(
    struct dhist *histogram, struct thread_info *info, int mp_flag);
static void schedule_bucket_destruction(
    struct dhist *histogram, struct bucket *bucket, int mp_flag);

static bool is_target_boundary(struct bucket *bucket, fixedpt observation);
static int compute_balance(struct bucket *bucket);
static void obtain_write_lock(struct bucket *bucket, int mp_flag);
static void release_write_lock(struct bucket *bucket, int mp_flag);
static void lock_boundary(struct bucket *bucket, int mp_flag);
static bool trylock_boundary(struct bucket *bucket, int mp_flag);
static void unlock_boundary(struct bucket *bucket, int mp_flag);
static void fix_height(struct bucket *bucket);
static void set_child(struct bucket *root, struct bucket *child, int dir);

static fixedpt split_threshold(struct dhist *histogram);
static fixedpt get_decay(struct dhist *histogram, uint64_t missed_generations);
static uint64_t get_next_generation(
    struct dhist *histogram, bool increment, int mp_flag);
static fixedpt compute_bound(
    struct dhist *histogram, struct bucket *left, struct bucket *right);
static fixedpt compute_count(
    struct dhist *histogram, struct bucket *bucket, uint64_t generation);
static bool perform_add(
    struct dhist *histogram, struct bucket *bucket,
    fixedpt observation, bool recheck_left_boundary,
    bool recheck_right_boundary, int mp_flag);
static struct bucket * find_and_lock_bucket(
    struct dhist *histogram, fixedpt observation, int mp_flag);
static void fix_balance_enstack(
    struct dhist *histogram, struct bucket *bucket);
static void fix_balance_process_item(struct dhist *histogram);
static void fix_balance(struct dhist *histogram, struct bucket *bucket);
static struct bucket * rotate_single(
    struct dhist *histogram, struct bucket *root, int dir);
static struct bucket * rotate_fixedpt(
    struct dhist *histogram, struct bucket *root, int dir);
static void handle_bucket_split_and_deletes(
    struct dhist *histogram, struct bucket *bucket, int mp_flag);
static void handle_bucket_deletes(struct dhist *histogram, int mp_flag);
static void split_bucket(
    struct dhist *histogram, struct bucket *bucket, int mp_flag);
static void delete_bucket(
    struct dhist *histogram, struct bucket *bucket, int mp_flag);
static void decay(
    struct dhist *histogram, struct bucket *bucket, uint64_t generation);
static int snprint_histogram(
    char *s_buffer, size_t n, struct dhist_info *info, const char *title,
    const char *xlabel);
static inline fixedpt roundoff_error_in_bounds(
    fixedpt value, fixedpt lower, fixedpt higher);

static void extract_info(
    struct dhist *histogram, int mp_flag, struct dhist_info *info);
static void clean_info(struct dhist_info *info);


const int DHIST_SINGLE_THREADED = 1 << 0;
const int DHIST_MULTI_THREADED = 1 << 1;


struct dhist *
dhist_init(uint32_t target_buckets, fixedpt decay_rate) {
  struct dhist *histogram;
  uint64_t idx;

  histogram = (struct dhist *)malloc(sizeof(struct dhist));

  histogram->decay_rate = decay_rate;
  // total_count should approach 1.0 / (1.0 - decay_rate) since it's a
  // geometric series.
  histogram->total_count = 0.0;

  histogram->target_num_buckets = target_buckets;
  histogram->fix_balance_stack = NULL;

  // We're single threaded in initialization.
  histogram->root = bucket_init(DHIST_SINGLE_THREADED);
  histogram->num_buckets = 1;
  histogram->generation = 0;
  // The pow_table is a cache of the commonly used compound decay factors.
  histogram->num_precomputed_powers = NUM_PRECOMPUTED_POWERS;
  histogram->pow_table = (fixedpt *)malloc(
      histogram->num_precomputed_powers * sizeof(fixedpt));
  for (idx = 0; idx < histogram->num_precomputed_powers; idx++)
    histogram->pow_table[idx] = ipow(decay_rate, idx);

  histogram->thread_info_head = NULL;

  histogram->tree_mtx = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));
  histogram->generation_mtx =
      (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));
  histogram->thread_info_mtx =
      (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));

  pthread_mutex_init(histogram->tree_mtx, NULL);
  pthread_mutex_init(histogram->generation_mtx, NULL);
  pthread_mutex_init(histogram->thread_info_mtx, NULL);

  return histogram;
}

void dhist_destroy(struct dhist *histogram) {
  struct bucket *bucket, *memo;
  struct thread_info *info, *info_next;

  // Select leftmost bucket.
  bucket = histogram->root;
  while (bucket->children[0])
    bucket = bucket->children[0];

  while (bucket) {
    memo = bucket;
    bucket = bucket->above;
    bucket_destroy(memo);
  }

  info = histogram->thread_info_head;
  while (info) {
    info_next = info->next;
    free(info->bucket_to_free);
    free(info);
    info = info_next;
  }

  free(histogram->pow_table);
  free(histogram->tree_mtx);
  free(histogram->generation_mtx);
  free(histogram->thread_info_mtx);
  free(histogram);
}

void dhist_insert(
    struct dhist *histogram, uint64_t observation, int mp_flag) {
  fixedpt obs;
  struct bucket *bucket;
  bool add_succeeded;
  fixedpt boundary;
  struct thread_info info;

  obs = 

  thread_info_init_fields(histogram, &info, mp_flag);

  do {
    bucket = find_and_lock_bucket(histogram, observation, mp_flag);

    // We have a read lock on the bucket that we need to update, but we don't
    // have a write lock until we have both boundary mutexes. If we can't grab
    // the other lock, we might be in danger of a deadlock. To avoid that, if
    // we fail to take the other lock and that lock is below->boundary_mtx,
    // we'll drop our lock and try again.

    if (!bucket->below) {
      // observation is less than leftmost_bucket->mu; always insert the
      // observaiton into this bucket.
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
      if (trylock_boundary(bucket->below, mp_flag)) {
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
      add_succeeded = perform_add(
          histogram, bucket, observation, false, false, mp_flag);
    }

    unlock_boundary(bucket, mp_flag);
  } while (!add_succeeded);

  if (bucket->data->count > split_threshold(histogram))
    handle_bucket_split_and_deletes(histogram, bucket, mp_flag);

  thread_info_finalize(histogram, &info, mp_flag);
}

// The caller is expected to manage the memory associated with s_buffer. The
// return value is one fewer than the number of bytes needed to write the data.
int dhist_snprint_histogram(
    char *s_buffer, size_t n, struct dhist *histogram, const char *title,
    const char *xlabel, int mp_flag) {
  struct dhist_info info;

  extract_info(histogram, mp_flag, &info);
  n = (size_t)snprint_histogram(s_buffer, n, &info, title, xlabel);
  clean_info(&info);
  return (int)n;
}

void dhist_set_num_buckets(
    struct dhist *histogram, uint32_t target_buckets, int mp_flag) {
  uint32_t old_num_buckets;
  struct thread_info info;

  thread_info_init_fields(histogram, &info, mp_flag);

  old_num_buckets = histogram->num_buckets;
  histogram->target_num_buckets = target_buckets;

  if (old_num_buckets > target_buckets) {
    if (mp_flag & DHIST_MULTI_THREADED)
      pthread_mutex_lock(histogram->tree_mtx);

    handle_bucket_deletes(histogram, mp_flag);

    if (mp_flag & DHIST_MULTI_THREADED)
      pthread_mutex_unlock(histogram->tree_mtx);
  }

  thread_info_finalize(histogram, &info, mp_flag);
}

uint32_t dhist_get_num_buckets(
    struct dhist *histogram, bool get_actual_instead_of_target) {
  if (get_actual_instead_of_target)
    return histogram->num_buckets;
  else
    return histogram->target_num_buckets;
}

// This is only safe if called from a single-threaded context.
void dhist_set_decay_rate(struct dhist *histogram, fixedpt decay_rate) {
  struct bucket *cursor;
  uint64_t generation, idx;
  fixedpt count, old_total_count, new_total_count, fudge_factor;

  // Can't do it.
  if (decay_rate <= 0.0 || 1.0 <= decay_rate)
    return;

  generation = get_next_generation(histogram, false, DHIST_SINGLE_THREADED);
  old_total_count = histogram->total_count;
  new_total_count =
      1.0 + decay_rate * (1.0 - ipow(decay_rate, generation - 1)) /
      (1.0 - decay_rate);
  // fudge_factor: The ratio of the right answer to the wrong answer.
  fudge_factor = new_total_count / old_total_count;

  for (idx = 0; idx < histogram->num_precomputed_powers; idx++)
    histogram->pow_table[idx] = ipow(decay_rate, idx);

  cursor = histogram->root;
  while (cursor->children[0])
    cursor = cursor->children[0];

  count = 0.0;
  while (cursor) {
    decay(histogram, cursor, generation);
    cursor->data->count *= fudge_factor;
    count += cursor->data->count;
    cursor = cursor->above;
  }

  histogram->total_count = count;
  histogram->decay_rate = decay_rate;
}

fixedpt dhist_get_decay_rate(struct dhist *histogram) {
  return histogram->decay_rate;
}

static fixedpt
ipow(fixedpt coefficient, uint64_t power) {
  fixedpt result;

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
bucket_init(int mp_flag) {
  struct bucket *bucket;

  bucket = (struct bucket *)malloc(sizeof(struct bucket));
  bucket->data = (struct bucket_data *)malloc(sizeof(struct bucket));

  bucket->boundary_mtx =
    (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));
  pthread_mutex_init(bucket->boundary_mtx, NULL);
  bucket->data->fix_balance_stack_next = NULL;

  bucket->data->is_enabled = true;
  bucket->data->lock_held = false;
  bucket->data->count = 0.0;
  bucket->data->mu = 0.0;
  bucket->data->update_generation = 0;
  bucket->data->height = 1;
  bucket->data->propogate_rebalance = false;
  bucket->data->fix_balance_stack_next = NULL;

  bucket->below = NULL;
  bucket->above = NULL;
  bucket->parent = NULL;
  bucket->children[0] = NULL;
  bucket->children[1] = NULL;
  unlock_boundary(bucket, mp_flag);

  return bucket;
}

void
bucket_destroy(struct bucket *bucket) {
  free(bucket->data);
  free(bucket->boundary_mtx);
  free(bucket);
}

// This trio of functions (along with thread_info_finalize and
// schedule_bucket_destruction) allow the program to remove buckets from the
// tree in a safe way. This is done by requiring that every thread register
// itself as it enters dhist_insert and unregister itself when it leaves; if
// it's the only thread that could have had a reference to a stale bucket, it
// is then required to free that bucket. This requires two additional locked
// regions.
//
// The time required to take the locks appears to be about 15% of the entire
// execution time. However, removing this feature should provide a 30% speedup.
// This disparity is likely due to cache misses, but that's only a hypothesis.
//
// If the unit of time is discrete (every dhist_insert increments the
// generation counter), it would be preferable to create a pool of buckets
// since the maximum number of buckets is known at initialization time.
// However, this doesn't work if the unit of time is continuous.
//
// It's desirable to permit a continuous unit of time because it's possible to
// get a consistent clock across several processors. This will make it possible
// to create a dhist for each processor. Even if it's not possible to make each
// of those be single threaded, doing this will eliminate many of the cache
// misses.
//
// The process of merging multiple continuous time histograms together doesn't
// lose much information, which makes this route desirable.
void
thread_info_init_fields(
    struct dhist *histogram, struct thread_info *info, int mp_flag) {
  info->next = NULL;
  info->bucket_to_free = NULL;

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_lock(histogram->thread_info_mtx);

  info->prev = NULL;
  info->next = histogram->thread_info_head;
  if (histogram->thread_info_head)
    histogram->thread_info_head->prev = info;
  histogram->thread_info_head = info;

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_unlock(histogram->thread_info_mtx);
}

void
thread_info_finalize(
    struct dhist *histogram, struct thread_info *info, int mp_flag) {
  struct thread_info *prior, *post, *memo;

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_lock(histogram->thread_info_mtx);

  memo = info;
  prior = info->prev;
  post = info->next;
  if (post == NULL) {
    while (prior && prior->bucket_to_free) {
      // If we have work to do, pull out that chain.
      memo = prior;
      prior = memo->prev;
    }
  }

  if (histogram->thread_info_head == memo)
    histogram->thread_info_head = post;

  memo->prev = NULL;
  if (prior)
    prior->next = post;
  if (post)
    post->prev = prior;

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_unlock(histogram->thread_info_mtx);

  // Delete any excess buckets and their info tokens. Don't free info since
  // it's stack allocated.
  memo = info->prev;
  while (memo) {
    if (memo->bucket_to_free)
      bucket_destroy(memo->bucket_to_free);
    prior = memo->prev;
    free(memo);
    memo = prior;
  }
}

// Schedule the bucket for destruction (this will happen as some thread
// processes thread_info_finalize). This does not unlock the boundaries, but
// it does set all relevant fields and flags. The caller needs to hold
// tree_mtx.
static void schedule_bucket_destruction(
    struct dhist *histogram, struct bucket *bucket, int mp_flag) {
  struct thread_info *info;

  // Prepare and schedule bucket for destruction.
  bucket->data->is_enabled = false;

  bucket->below = bucket->above = NULL;
  bucket->children[0] = bucket->children[1] = NULL;
  info = (struct thread_info *)malloc(sizeof(struct thread_info));
  info->prev = NULL;
  info->bucket_to_free = bucket;

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_lock(histogram->thread_info_mtx);

  info->next = histogram->thread_info_head;
  histogram->thread_info_head->prev = info;
  histogram->thread_info_head = info;

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_unlock(histogram->thread_info_mtx);
}

/*
 * Returns true if the observation is between the lower bucket mu (inclusive)
 * and the current bucket mu (exclusive).
 * obs \in [bucket->below->data->mu, bucket->data->mu)
 */
static bool
is_target_boundary(struct bucket *bucket, fixedpt observation) {
  fixedpt left_mu, mu;
  struct bucket *below;

  if ((below = bucket->below) != NULL) {
    left_mu = below->data->mu;
  } else {
    left_mu = observation - 1.0;
  }

  mu = bucket->data->mu;
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

// Identifies the difference in heights of the two sub-trees beneath this node.
// A possitive return means the right sub-branch is taller.
static int
compute_balance(struct bucket *bucket) {
  int dir;
  int heights[2];

  for (dir = 0; dir <= 1; dir++) {
    if (bucket->children[dir] == NULL)
      heights[dir] = 0;
    else
      heights[dir] = bucket->children[dir]->data->height;
  }

  return heights[1] - heights[0];
}

// Obtains both locks associated with a bucket.
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
    if (trylock_boundary(second, mp_flag)) {
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
    bucket->data->lock_held = true;
  }
}

static bool
trylock_boundary(struct bucket *bucket, int mp_flag) {
  if (mp_flag & DHIST_SINGLE_THREADED) {
    return true;
  } else if (pthread_mutex_trylock(bucket->boundary_mtx) == 0) {
    bucket->data->lock_held = true;
    return true;
  } else {
    return false;
  }
}

static void
unlock_boundary(struct bucket *bucket, int mp_flag) {
  if (mp_flag & DHIST_MULTI_THREADED) {
    bucket->data->lock_held = false;
    pthread_mutex_unlock(bucket->boundary_mtx);
  }
}

// Rebalances this node in the AVL tree.
static void
fix_height(struct bucket *bucket) {
  struct bucket *left, *right;
  int prior_height;

  while (bucket != NULL) {
    left = bucket->children[0];
    right = bucket->children[1];

    prior_height = bucket->data->height;
    if (left == NULL && right == NULL) {
      bucket->data->height = 1;
    } else if (left == NULL) {
      bucket->data->height = 1 + right->data->height;
    } else if (right == NULL) {
      bucket->data->height = 1 + left->data->height;
    } else {
      bucket->data->height = 1 +
          (left->data->height > right->data->height ?
           left->data->height : right->data->height);
    }

    bucket = bucket->parent;
  }
}

static void set_child(struct bucket *root, struct bucket *child, int dir) {
  root->children[dir] = child;
  if (child != NULL)
    child->parent = root;
  fix_height(root);
}

// Calculate the decay factor that should be applied to a count that hasn't
// bene updated in missed_generations updates.
static fixedpt get_decay(struct dhist *histogram, uint64_t missed_generations) {
  if (missed_generations < histogram->num_precomputed_powers)
    return histogram->pow_table[missed_generations];
  else
    return ipow(histogram->decay_rate, missed_generations);
}

static fixedpt split_threshold(struct dhist *histogram) {
  //return
  //    2 * (1.0 / (1.0 - histogram->decay_rate)) /
  //    histogram->target_num_buckets;
  return 2 * (histogram->total_count / histogram->target_num_buckets);
}

// Get the next generation count and increment the global counter.
static uint64_t get_next_generation(
    struct dhist *histogram, bool increment, int mp_flag) {
  uint64_t generation;

  if (mp_flag & DHIST_MULTI_THREADED) {
    pthread_mutex_lock(histogram->generation_mtx);
    if (increment) {
      ++histogram->generation;
      histogram->total_count =
          1.0 + histogram->total_count * histogram->decay_rate;
    }
    generation = histogram->generation;
    pthread_mutex_unlock(histogram->generation_mtx);
  } else {
    if (increment) {
      ++histogram->generation;
      histogram->total_count =
          1.0 + histogram->total_count * histogram->decay_rate;
    }
    generation = histogram->generation;
  }

  return generation;
}

// Compute the boundary between two buckets.
static fixedpt
compute_bound(
    struct dhist *histogram, struct bucket *left, struct bucket *right) {
  uint64_t generation;
  fixedpt count_left, mu_left, count_right, mu_right, ret;

  // We can't increment generation until we finally make our update, but since
  // the boundaries don't move upon decay, we can use any generation value to
  // compute the bounds. We'll use the largest generation so that we don't need
  // to recompute the decay for one of the buckets.
  if (!left) {
    return right->data->mu - 0.5;
  } else if (!right) {
    return left->data->mu + 0.5;
  }

  if (left->data->update_generation > right->data->update_generation) {
    generation = left->data->update_generation;
    mu_left = left->data->mu;
    mu_right = right->data->mu;
    count_left = left->data->count;
    count_right = compute_count(histogram, right, generation);
  } else {
    generation = right->data->update_generation;
    mu_left = left->data->mu;
    mu_right = right->data->mu;
    count_left = compute_count(histogram, left, generation);
    count_right = right->data->count;
  }

  ret = (
      (mu_left * count_left + mu_right * count_right) /
      (count_left + count_right));
  return roundoff_error_in_bounds(ret, mu_left, mu_right);
}

// Compute the count in a bucket.
static fixedpt compute_count(
    struct dhist *histogram, struct bucket *bucket, uint64_t generation) {
  return (
      bucket->data->count *
      get_decay(histogram, generation - bucket->data->update_generation));
}

// The caller should have a write lock on the bucket. This can fail if one of
// the recheck_* flags is true and the boundary is not correct.
static bool
perform_add(
    struct dhist *histogram, struct bucket *bucket, fixedpt observation,
    bool recheck_left_boundary, bool recheck_right_boundary, int mp_flag) {
  fixedpt boundary, mu, mu_lower_bound, mu_upper_bound;
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

  update_generation = get_next_generation(histogram, true, mp_flag);
  bucket->data->count = 1.0 + compute_count(
      histogram, bucket, update_generation);
  bucket->data->update_generation = update_generation;

  mu = (bucket->data->count * bucket->data->mu + observation) /
       (bucket->data->count + 1);
  mu_lower_bound = bucket->below ? bucket->below->data->mu : mu;
  mu_upper_bound = bucket->above ? bucket->above->data->mu : mu;
  bucket->data->mu =
      roundoff_error_in_bounds(mu, mu_lower_bound, mu_upper_bound);

  return true;
}

/*
 * Return the bucket with the greatest mu less than the observation (which is
 * not necessarily the bucket where observation will be inserted) unless the
 * target bucket is the leftmost bucket, in which case observation may be less
 * than bucket->data->mu.
 * If mp_flag & DHIST_MULTI_THREADED, the boundary_mtx for this bucket and its
 * right hand neighbor will be locked.
 */
static struct bucket * find_and_lock_bucket(
    struct dhist *histogram, fixedpt observation, int mp_flag) {
  struct bucket *bucket, *other;

  do {
    bucket = histogram->root;
    do {
      // The first two branches are quick decisions that avoid computing
      // the region where the boundary can exist.
      if ((other = bucket->below) && observation < other->data->mu) {
        bucket = bucket->children[0];
      } else if ((other = bucket->above) && other->data->mu <= observation) {
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

          if (!is_target_boundary(bucket, observation) ||
              !bucket->data->is_enabled) {
            // We raced, so restart.
            unlock_boundary(bucket, mp_flag);
            bucket = histogram->root;
          } else {
            return bucket;
          }
        }
      } else {
        bucket = bucket->children[bucket->data->mu <= observation];
      }
    } while (bucket != NULL);
    // We can fail to find a bucket if the bucket boundaries moved around
    // during the search. Start over.
  } while (1);
  return NULL;
}

// Add a bucket to the AVL tree rebalance stack.
static void
fix_balance_enstack(struct dhist *histogram, struct bucket *bucket) {
  if (bucket) {
    bucket->data->fix_balance_stack_next = histogram->fix_balance_stack;
    histogram->fix_balance_stack = bucket;
  }
}

static void
fix_balance_process_item(struct dhist *histogram) {
  int balance;
  int prior_height;
  int dir;
  struct bucket *bucket;

  bucket = histogram->fix_balance_stack;
  if (bucket == NULL)
    return;

  fix_height(bucket);
  prior_height = bucket->data->height;
  balance = compute_balance(bucket);
  if (ABS(compute_balance(bucket)) >= 2) {
    balance = compute_balance(bucket);
    dir = (balance < 0);

    if (bucket->children[!dir]->children[!dir] == NULL) {
      bucket = rotate_fixedpt(histogram, bucket, dir);
    } else if (bucket->children[!dir]->children[dir] == NULL) {
      bucket = rotate_single(histogram, bucket, dir);
    } else if (bucket->children[!dir]->children[dir]->data->height >
        bucket->children[!dir]->children[!dir]->data->height) {
      bucket = rotate_fixedpt(histogram, bucket, dir);
    } else {
      bucket = rotate_single(histogram, bucket, dir);
    }
  }
  fix_height(bucket);

  histogram->fix_balance_stack = bucket->data->fix_balance_stack_next;
  bucket->data->fix_balance_stack_next = NULL;

  if (bucket->data->propogate_rebalance &&
      prior_height == bucket->data->height) {
    bucket->data->propogate_rebalance = false;
    bucket = bucket->parent;
    if (bucket)
      bucket->data->propogate_rebalance = true;
    fix_balance_enstack(histogram, bucket);
  }
}

static void
fix_balance(
    struct dhist *histogram, struct bucket *bucket) {
  bucket->data->propogate_rebalance = true;
  fix_balance_enstack(histogram, bucket);

  while (histogram->fix_balance_stack) {
    fix_balance_process_item(histogram);
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

  fix_balance_enstack(histogram, new_root);
  fix_balance_enstack(histogram, new_root->children[0]);
  fix_balance_enstack(histogram, new_root->children[1]);

  return new_root;
}

static struct bucket *
rotate_fixedpt(struct dhist *histogram, struct bucket *root, int dir) {
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

  fix_balance_enstack(histogram, new_root);
  fix_balance_enstack(histogram, new_root->children[0]);
  fix_balance_enstack(histogram, new_root->children[1]);

  return new_root;
}

// Split the bucket. If we grew past target_num_buckets, find the least
// populated bucket and schedule it for destruction.
void
handle_bucket_split_and_deletes(
    struct dhist *histogram, struct bucket *bucket, int mp_flag) {
  if (mp_flag & DHIST_MULTI_THREADED) {
    pthread_mutex_lock(histogram->tree_mtx);

    // Did we race?
    if (bucket->data->is_enabled == false ||
        bucket->data->count <= split_threshold(histogram)) {
      pthread_mutex_unlock(histogram->tree_mtx);
      return;
    }
  }

  split_bucket(histogram, bucket, mp_flag);
  handle_bucket_deletes(histogram, mp_flag);

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_unlock(histogram->tree_mtx);
}

// This should be called from inside a context that already holds the tree_mtx.
void handle_bucket_deletes(struct dhist *histogram, int mp_flag) {
  struct bucket *cursor, *min_bucket;

  while (histogram->num_buckets > histogram->target_num_buckets) {
    cursor = histogram->root;
    while (cursor->children[0])
      cursor = cursor->children[0];

    min_bucket = NULL;
    while (cursor) {
      if (cursor->data->is_enabled &&
          (min_bucket == NULL ||
           cursor->data->count < min_bucket->data->count)) {
        min_bucket = cursor;
      }
      cursor = cursor->above;
    }

    delete_bucket(histogram, min_bucket, mp_flag);
  }
}

static void
split_bucket(struct dhist *histogram, struct bucket *bucket, int mp_flag) {
  fixedpt lower_bound, upper_bound;
  struct bucket *new_bucket, *memo;

  obtain_write_lock(bucket, mp_flag);

  new_bucket = bucket_init(mp_flag);
  ++histogram->num_buckets;
  new_bucket->parent = bucket;

  bucket->data->count /= 2.0;
  new_bucket->data->count = bucket->data->count;
  new_bucket->data->update_generation = bucket->data->update_generation;
  new_bucket->below = bucket->below;
  bucket->below = new_bucket;

  if (histogram->num_buckets == 2) {
    // This is an awkward split because we don't have enough information to
    // compute a diameter. Instead, we can set mu to the same value in both
    // buckets, which will also make the boundary between them be mu. As soon
    // as we observe a non-mu entry, the split between the buckets
    // will be more sensical. However, we will need to take care to not split
    // a bucket again until this happens.
    new_bucket->data->mu = bucket->data->mu;
  } else {
    if (!new_bucket->below) {
      upper_bound = compute_bound(histogram, bucket, bucket->above);
      // If no lower bound, construct one such that mu is in the exact center.
      lower_bound = bucket->data->mu - (upper_bound - bucket->data->mu);
    } else if (!bucket->above) {
      lower_bound = compute_bound(histogram, new_bucket->below, bucket);
      upper_bound = bucket->data->mu + (bucket->data->mu - lower_bound);
    } else {
      lower_bound = compute_bound(histogram, new_bucket->below, bucket);
      upper_bound = compute_bound(histogram, bucket, bucket->above);
    }

    new_bucket->data->mu = (lower_bound + bucket->data->mu) / 2.0;
    bucket->data->mu = (bucket->data->mu + upper_bound) / 2.0;
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

  fix_balance(histogram, new_bucket);
}

static void
delete_bucket(struct dhist *histogram, struct bucket *bucket, int mp_flag) {
  struct bucket *attach, *lucky_bucket, *other;
  struct bucket *far_buckets[2];
  int dir_from_parent, dir, balance, idx;
  uint64_t generation;
  fixedpt below_count, above_count, mu, mu_lower_bound, mu_upper_bound;

 restart:
  obtain_write_lock(bucket, mp_flag);
  generation = get_next_generation(histogram, false, mp_flag);

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
          !trylock_boundary(far_buckets[idx], mp_flag)) {
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
  decay(histogram, bucket, generation);
  decay(histogram, lucky_bucket, generation);

  // Remove bucket from the linked list.
  if (bucket->below)
    bucket->below->above = bucket->above;
  if (bucket->above)
    bucket->above->below = bucket->below;

  // Update the lucky_bucket's stats.
  mu = (lucky_bucket->data->mu * lucky_bucket->data->count +
        bucket->data->mu * bucket->data->count) /
       (lucky_bucket->data->count + bucket->data->count);
  mu_lower_bound = lucky_bucket->below ? lucky_bucket->below->data->mu : mu;
  mu_upper_bound = lucky_bucket->above ? lucky_bucket->above->data->mu : mu;
  lucky_bucket->data->mu =
      roundoff_error_in_bounds(mu, mu_lower_bound, mu_upper_bound);

  lucky_bucket->data->count += bucket->data->count;

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
    fix_balance(histogram, attach);

  other = bucket->above;
  schedule_bucket_destruction(histogram, bucket, mp_flag);
  if (mp_flag & DHIST_MULTI_THREADED) {
    unlock_boundary(bucket, mp_flag);
    if (other)
      unlock_boundary(other, mp_flag);
  }

  --histogram->num_buckets;
}

// Update the bucket count to reflect the amount of decay that occurs up
// until the supplied generation.
static void
decay(struct dhist *histogram, struct bucket *bucket, uint64_t generation) {
  if (bucket == NULL || bucket->data->update_generation == generation)
    return;

  bucket->data->count *= get_decay(
      histogram, generation - bucket->data->update_generation);
  bucket->data->update_generation = generation;
}

// Scan through the histogram and record the information associated with the
// boundaries and counts.
static void extract_info(
    struct dhist *histogram, int mp_flag, struct dhist_info *info) {
  struct bucket *bucket;
  fixedpt total_count;
  int idx;

  info->histogram = histogram;

  if (mp_flag & DHIST_MULTI_THREADED)
    pthread_mutex_lock(histogram->tree_mtx);

  // Allocate enough space for both the weights and the boundaries. We can
  // have one more than the target number of buckets for a short period of
  // time, but we can't grow past that.
  info->weights = (fixedpt *)malloc(
      sizeof(fixedpt) * 2 * (histogram->target_num_buckets + 1) + 1);
  info->boundaries = info->weights + histogram->target_num_buckets + 1;

  info->generation = get_next_generation(histogram, false, mp_flag);

  bucket = histogram->root;
  // Astonishingly, this is possibly actually slower than traversing
  // the bucket->below linked list.
  while (bucket->children[0])
    bucket = bucket->children[0];

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
    if (bucket->data->update_generation < info->generation)
      decay(histogram, bucket, info->generation);
    info->weights[idx] = bucket->data->count;

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

// Free the internal memory associated with a dhist_info object.
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
  int idx, num_chars, num_for_call;
  size_t remaining_chars;

#define dhist_snprint_helper() do {               \
    num_chars += num_for_call;                    \
    if (num_for_call > (int)remaining_chars) {    \
      remaining_chars = 0;                        \
      s_buffer = NULL;                            \
    } else {                                      \
      s_buffer = s_buffer + num_for_call;         \
      remaining_chars -= (size_t)num_for_call;    \
    }                                             \
  } while (0)

  num_chars = 0;
  remaining_chars = n;

  num_for_call = snprintf(s_buffer, remaining_chars, "{");
  dhist_snprint_helper();

  if (title) {
    num_for_call =
        snprintf(s_buffer, remaining_chars, "\"title\": \"%s\", ", title);
    dhist_snprint_helper();
  }
  if (xlabel) {
    num_for_call =
        snprintf(s_buffer, remaining_chars, "\"xlabel\": \"%s\", ", xlabel);
    dhist_snprint_helper();
  }

  num_for_call = snprintf(
      s_buffer, remaining_chars, "\"generation\": %lu, ", info->generation);
  dhist_snprint_helper();
  num_for_call = snprintf(
      s_buffer, remaining_chars, "\"id\": %lu, ", (uint64_t)info->histogram);
  dhist_snprint_helper();

  num_for_call = snprintf(s_buffer, remaining_chars, "\"weights\": [");
  dhist_snprint_helper();
  for (idx = 0; idx < info->num_buckets - 1; idx++) {
    num_for_call = snprintf(
        s_buffer, remaining_chars, "%lf, ", info->weights[idx]);
    dhist_snprint_helper();
  }
  num_for_call = snprintf(
      s_buffer, remaining_chars, "%lf], ", info->weights[idx]);
  dhist_snprint_helper();

  num_for_call = snprintf(s_buffer, remaining_chars, "\"boundaries\": [");
  dhist_snprint_helper();
  for (idx = 0; idx < info->num_buckets; idx++) {
    num_for_call = snprintf(
        s_buffer, remaining_chars, "%lf, ", info->boundaries[idx]);
    dhist_snprint_helper();
  }
  num_for_call = snprintf(
      s_buffer, remaining_chars, "%lf]}", info->boundaries[idx]);
  dhist_snprint_helper();

#undef dhist_snprint_helper

  return num_chars;
}

// Some operations can produce roundoff error the makes a value shift outside
// of sensible bounds. For example, the mu for a bucket could be smaller than
// its left neighbor due to roundoff error. If this phenomenon exists, replace
// the value with the acceptible boundary value.
fixedpt roundoff_error_in_bounds(fixedpt value, fixedpt lower, fixedpt upper) {
  if (lower <= value && value <= upper) {
    return value;
  } else if (value < lower) {
    return lower;
  } else {
    return upper;
  }
}

