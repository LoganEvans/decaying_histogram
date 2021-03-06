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
#include <stdio.h>
#include <iostream>
#include <random>

#include "src/dhist.c"
#include "gtest/gtest.h"

using std::cout;
using std::endl;

// TODO: Move the struct definitions in src/dhist.c into a header so that
// these functions can be offloaded to a separate compilation unit.
static void _print_tree(struct bucket *bucket, int depth);
static void print_tree(struct bucket *bucket);
static void assert_consistent(struct dhist *histogram);
static int assert_invariant(struct bucket *root);
static int count_buckets_in_tree(struct bucket *root);

class HistogramTest : public testing::Test {
 protected:
  virtual void SetUp() {
    decay_rate_ = 0.999;
    target_buckets_ = 40;
    histogram_ = dhist_init(target_buckets_, decay_rate_);
    generator_.seed(0);
  }

  virtual void TearDown() {
    dhist_destroy(histogram_);
  }

  struct dhist *histogram_;
  double decay_rate_;
  int target_buckets_;
  std::default_random_engine generator_;
  std::normal_distribution<double> normal_;
};

TEST_F(HistogramTest, Stress) {
  for (int idx = 0; idx < 100000; idx++) {
    dhist_insert(histogram_, normal_(generator_), DHIST_SINGLE_THREADED);
    assert_consistent(histogram_);
  }
}

TEST_F(HistogramTest, PositiveDensity) {
  double observation;

  for (int idx = 0; idx < 100000; idx++) {
    observation = normal_(generator_);
    dhist_insert(histogram_, observation, DHIST_SINGLE_THREADED);
  }

  struct bucket *cursor = histogram_->root;
  while (cursor->children[0])
    cursor = cursor->children[0];

  while (cursor) {
    EXPECT_GT(cursor->data->count, 0.0);
    cursor = cursor->above;
  }
}

TEST_F(HistogramTest, TotalCount) {
  double count;
  int iterations = 100000;

  count = 0.0;
  for (int idx = 0; idx < iterations; idx++) {
    // The total count should match this update pattern.
    count *= decay_rate_;
    count += 1.0;
    dhist_insert(histogram_, normal_(generator_), DHIST_SINGLE_THREADED);
    ASSERT_NEAR(count, histogram_->total_count, 1e-10);
  }

  EXPECT_LT(0, histogram_->num_buckets);

  // The buckets should sum to the total count.

  struct bucket *cursor = histogram_->root;
  while (cursor->children[0])
    cursor = cursor->children[0];

  int idx = 0;
  count = 0.0;
  while (cursor) {
    idx++;
    count += compute_count(histogram_, cursor, histogram_->generation);
    cursor = cursor->above;
  }

  ASSERT_NEAR(count, histogram_->total_count, 1e-10);

  ASSERT_NEAR(
      histogram_->total_count,
      1.0 + decay_rate_ * (1.0 - ipow(decay_rate_, iterations - 1)) /
        (1.0 - decay_rate_),
      1e-10);
}

TEST_F(HistogramTest, InsertConstant) {
  // This should not freeze if a constant value is inserted.
  for (int idx = 0; idx < 100000; idx++)
    dhist_insert(histogram_, 1, DHIST_SINGLE_THREADED);

  // We can also stress it a bit by first making sure several buckets exist.
  for (int idx = 0; idx < 100000; idx++)
    dhist_insert(histogram_, normal_(generator_), DHIST_SINGLE_THREADED);
  for (int idx = 0; idx < 100000; idx++)
    dhist_insert(histogram_, 1, DHIST_SINGLE_THREADED);

  // Let's also be obnoxious and use a painfully quick decay rate.
  dhist_destroy(histogram_);
  histogram_ = dhist_init(50, 0.9);
  for (int idx = 0; idx < 10000; idx++)
    dhist_insert(histogram_, normal_(generator_), DHIST_SINGLE_THREADED);
  for (int idx = 0; idx < 10000; idx++)
    dhist_insert(histogram_, 1, DHIST_SINGLE_THREADED);
}

TEST_F(HistogramTest, Trace1) {
  char buffer[100];
  int buffer_idx = 0;
  FILE *fin;
  char ch;

  fin = fopen("trace1.txt", "r");
  if (!fin)
    fin = fopen("tests/trace1.txt", "r");
  ASSERT_TRUE(fin);

  dhist_destroy(histogram_);
  histogram_ = dhist_init(20, 0.9);

  do {
    ch = static_cast<char>(fgetc(fin));
    if (ch == EOF || ch == ',') {
      buffer[buffer_idx] = '\0';
      if (buffer_idx > 0)
        dhist_insert(histogram_, atoi(buffer), DHIST_SINGLE_THREADED);
      assert_consistent(histogram_);
      buffer_idx = 0;
    } else {
      buffer[buffer_idx] = ch;
      buffer_idx++;
    }
  } while (ch != EOF);

  fclose(fin);
}

TEST_F(HistogramTest, SetNumBuckets) {
  // Burn in the histogram.
  for (int idx = 0; idx < 100000; idx++)
    dhist_insert(histogram_, normal_(generator_), DHIST_SINGLE_THREADED);

  ASSERT_EQ(target_buckets_, histogram_->target_num_buckets);
  ASSERT_EQ(target_buckets_, dhist_get_num_buckets(histogram_, true));
  ASSERT_EQ(target_buckets_, dhist_get_num_buckets(histogram_, false));

  // Now lower the number of buckets.
  dhist_set_num_buckets(
      histogram_, target_buckets_ - 5, DHIST_SINGLE_THREADED);

  // Buckets should immediately be deleted.
  ASSERT_EQ(target_buckets_ - 5, dhist_get_num_buckets(histogram_, true));
  ASSERT_EQ(target_buckets_ - 5, dhist_get_num_buckets(histogram_, false));

  // Raise the number of buckets.
  dhist_set_num_buckets(
      histogram_, target_buckets_ + 5, DHIST_SINGLE_THREADED);

  // We haven't had a chance to split buckets yet, so the actual and target
  // number of buckets will have a mismatch.
  ASSERT_EQ(target_buckets_ - 5, dhist_get_num_buckets(histogram_, true));
  ASSERT_EQ(target_buckets_ + 5, dhist_get_num_buckets(histogram_, false));

  // Burn in the new number of buckets. Shift the distribution so that we
  // actually need to split buckets.
  for (int idx = 0; idx < 10000; idx++)
    dhist_insert(histogram_, 9.0 + normal_(generator_), DHIST_SINGLE_THREADED);

  // Now the target and actual number of buckets should match.
  ASSERT_EQ(target_buckets_ + 5, dhist_get_num_buckets(histogram_, true));
  ASSERT_EQ(target_buckets_ + 5, dhist_get_num_buckets(histogram_, false));
}

TEST_F(HistogramTest, SetDecayRate) {
  double uut_decay_rate;
  ASSERT_EQ(histogram_->decay_rate, dhist_get_decay_rate(histogram_));

  for (int idx = 0; idx < 100000; idx++)
    dhist_insert(histogram_, normal_(generator_), DHIST_SINGLE_THREADED);

  uut_decay_rate = decay_rate_ / 2.0;
  dhist_set_decay_rate(histogram_, uut_decay_rate);
  ASSERT_NEAR(dhist_get_decay_rate(histogram_), uut_decay_rate, 1e-10);

  // Make sure we didn't break anything obvious.
  assert_consistent(histogram_);
  for (int idx = 0; idx < 1000; idx++) {
    dhist_insert(histogram_, normal_(generator_), DHIST_SINGLE_THREADED);
    assert_consistent(histogram_);
  }

  uut_decay_rate = (1.0 + decay_rate_) / 2.0;
  dhist_set_decay_rate(histogram_, uut_decay_rate);
  ASSERT_NEAR(dhist_get_decay_rate(histogram_), uut_decay_rate, 1e-10);

  for (int idx = 0; idx < 1000; idx++) {
    dhist_insert(histogram_, normal_(generator_), DHIST_SINGLE_THREADED);
    assert_consistent(histogram_);
  }
}

TEST_F(HistogramTest, RoundoffErrorInBounds) {
  double count_left, mu_left, count_right, mu_right, mu;

  count_left = 1.9415677466378186e-09;
  mu_left = 1279999.4617896655;
  count_right = 9.9539473717148734;
  mu_right = 1279999.9999999998;

  mu = ((count_left * mu_left) + (count_right * mu_right)) /
       (count_left + count_right);

  // Make sure that this DOES actually go the wrong way... the mean is larger
  // than either value.
  ASSERT_LT(mu_right, mu);
  ASSERT_EQ(
      mu_right,
      roundoff_error_in_bounds(mu, mu_left, mu_right));
}

static int count_buckets_in_tree(struct bucket *root) {
  if (root) {
    return 1 +
        count_buckets_in_tree(root->children[0]) +
        count_buckets_in_tree(root->children[1]);
  } else {
    return 0;
  }
}

static void _print_tree(struct bucket *bucket, int depth) {
  int i;

  printf("%.2lf", bucket->data->mu);

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

static void print_tree(struct bucket *bucket) {
  _print_tree(bucket, 0);
  printf("\n");
}

static void assert_consistent(struct dhist *histogram) {
  uint32_t num_buckets_seen = 0;
  double upper_bound, lower_bound, count = 0.0;
  struct bucket *cursor = histogram->root, *cursor_two;
  while (cursor->children[0])
    cursor = cursor->children[0];

  while (cursor != NULL) {
    count += compute_count(histogram, cursor, histogram->generation);
    num_buckets_seen++;
    upper_bound = compute_bound(histogram, cursor, cursor->above);
    lower_bound = compute_bound(histogram, cursor->below, cursor);
    if (cursor->above && cursor->data->mu > cursor->above->data->mu) {
      printf("ERROR: cursor->data->mu(%.19lf) > cursor->above->data->mu(%.19lf)\n",
          cursor->data->mu, cursor->above->data->mu);
      print_tree(histogram->root);
      assert(0);
    }
    if (upper_bound < lower_bound) {
      printf("ERROR: upper_bound(%.11lf) < lower_bound(%.11lf)\n",
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
      assert(cursor->below->data->mu <= cursor->data->mu);
    }

    if (cursor->children[1]) {
      assert(cursor->children[1]->parent == cursor);
      assert(cursor->above);
      assert(cursor->data->mu <= cursor->above->data->mu);
    }

    // Make sure we can find this bucket in the tree.
    cursor_two = histogram->root;
    while (true) {
      if (cursor_two == NULL) {
        printf("Could not find bucket with mu %lf in the tree.\n",
            cursor->data->mu);
        assert(false);
      } else if (cursor_two->data->mu == cursor->data->mu) {
        break;
      } else if (cursor->data->mu < cursor_two->data->mu) {
        cursor_two = cursor_two->children[0];
      } else {
        cursor_two = cursor_two->children[1];
      }
    }

    cursor = cursor->above;
  }

  if (ABS(count - histogram->total_count) >= 1e-9) {
    printf(
        "ERROR: Sum of bucket->data->count (%.11lf) does not equal "
        "total_count (%.11lf). Difference: %.11lf\n",
        count, histogram->total_count, ABS(count - histogram->total_count));
    assert(false);
  }
  assert(num_buckets_seen == histogram->num_buckets);
  assert(
      (uint64_t)count_buckets_in_tree(histogram->root) ==
      histogram->num_buckets);
  assert_invariant(histogram->root);
}

static int assert_invariant(struct bucket *root) {
  int left, right, dir;

  if (root == NULL) {
    return 0;
  } else if (root->children[0] == NULL && root->children[1] == NULL) {
    if (root->data->height != 1) {
      print_tree(root);
      assert(root->data->height == 1);
    }
    return 1;
  } else {
    if (root->children[0] && root->data->mu < root->children[0]->data->mu) {
      printf("ORDER ERROR(0): root->data->mu: %lf "
             "< root->children[0]->data->mu: %lf ...\n",
          root->data->mu, root->children[0]->data->mu);
      assert(false);
    }

    if (root->children[1] && root->data->mu > root->children[1]->data->mu) {
      printf("ORDER ERROR(1): root->data->mu: %lf "
             "> root->children[1]->data->mu: %lf ...\n",
          root->data->mu, root->children[1]->data->mu);
      assert(false);
    }

    left = assert_invariant(root->children[0]);
    if (root->children[0] == NULL && left != 0) {
      printf("ERROR(1): root->children[0] == NULL && left(%d) != 0\n", left);
      assert(false);
    } else if (root->children[0] &&
               (root->children[0]->data->height != left)) {
      printf("ERROR(2): root->children[0]->hieght(%d) != left(%d)\n",
          root->children[0]->data->height, left);
      assert(false);
    }

    right = assert_invariant(root->children[1]);
    if (root->children[1] == NULL && right != 0) {
      printf("ERROR(3): root->children[1] == NULL && right(%d) != 0\n", right);
      assert(false);
    } else if (root->children[1] &&
               (root->children[1]->data->height != right)) {
      printf("ERROR(4): root->children[1]->hieght(%d) != right(%d)\n",
          root->children[1]->data->height, right);
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

    if ((root->data->height != left + 1 && root->data->height != right + 1) ||
        (root->data->height <= left || root->data->height <= right)) {
      printf(
        "root height is not correct. heights -- root: %d left: %d right: %d\n",
        root->data->height, left, right);
      assert(false);
    }

    if (ABS(right - left) > 1) {
      assert(false);
    }

    return 1 + (left > right ? left : right);
  }
}

