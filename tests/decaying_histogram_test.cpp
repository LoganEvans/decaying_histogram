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

#include <iostream>
#include <random>

#include "decaying_histogram.h"
#include "gtest/gtest.h"


using std::cout;
using std::endl;

static std::default_random_engine g_generator;
static std::normal_distribution<double> g_normal;

class HistogramTest : public testing::Test {
 protected:
  virtual void SetUp() {
    alpha_ = 0.001;
    target_buckets_ = 40;
    histogram_ = new struct decaying_histogram;
    init_decaying_histogram(histogram_, target_buckets_, alpha_);
    g_generator.seed(0);
  }

  virtual void TearDown() {
    clean_decaying_histogram(histogram_);
    delete histogram_;
  }

  struct decaying_histogram *histogram_;
  double alpha_;
  int target_buckets_;
};

TEST_F(HistogramTest, TargetBoundary) {
  int num_buckets;
  double dens, lower_bound, upper_bound;
  struct bucket *bucket;
  bool target_boundary_seen;

#if 0
  for (num_buckets = 20; num_buckets <= 20;
      num_buckets++) {
    histogram_->num_buckets = num_buckets;
    histogram_->generation = 9001;
    for (int idx = 0; idx < num_buckets; idx++) {
      histogram_->bucket_list[idx].count = 1.0 / (histogram_->alpha * num_buckets);
      histogram_->bucket_list[idx].mu = idx;
      histogram_->bucket_list[idx].update_generation = 9001;
    }
    for (int idx = 0; idx < num_buckets - 1; idx++) {
      histogram_->bucket_list[idx].above = &histogram_->bucket_list[idx + 1];
    }
    for (int idx = 1; idx < num_buckets; idx++) {
      histogram_->bucket_list[idx].below = &histogram_->bucket_list[idx - 1];
    }
    histogram_->bucket_list[0].below = NULL;
    histogram_->bucket_list[num_buckets - 1].above = NULL;

    if (num_buckets == 2) {
      histogram_->bucket_list[0].above = &histogram_->bucket_list[1];
      histogram_->bucket_list[1].below = &histogram_->bucket_list[0];
    }

    for (int offset_num = -1; offset_num <= 1; offset_num++) {
      double offset = offset_num / 4.0;
      for (double center = 1; center < num_buckets - 1; center += 1.0) {
        target_boundary_seen = false;
        for (int idx = 0; idx < num_buckets; idx++) {
          bucket = &histogram_->bucket_list[idx];
          if (is_target_boundary(bucket, center + offset)) {
            EXPECT_FALSE(target_boundary_seen)
                << "buckets: " << num_buckets << " idx: " << idx << " loc: " << center + offset << endl;
            if (target_boundary_seen)
              assert(false);
            target_boundary_seen = true;
          }
        }
      }
      ASSERT_TRUE(target_boundary_seen);
    }
  }
#endif
}

TEST_F(HistogramTest, PositiveDensity) {
  double observation;

  for (int idx = 0; idx < 1000; idx++) {
    observation = g_normal(g_generator);
    dh_insert(histogram_, observation, DHIST_SINGLE_THREADED);
    //if (histogram_->num_buckets > 1)
    //  FAIL() << histogram_->bucket_list[0].lower_bound << " | " << histogram_->bucket_list[0].upper_bound
    //         << " ; " << histogram_->bucket_list[1].lower_bound << " | " << histogram_->bucket_list[1].upper_bound << endl;
  }

  for (int idx = 0; idx < histogram_->num_buckets; idx++) {
    ASSERT_GT(
        weight(histogram_, histogram_->bucket_list[idx], 0, NULL, NULL),
        0.0);
  }
}

TEST_F(HistogramTest, DensitySumsToOne) {
  double observation, acc, dens, lower_bound, upper_bound;

  for (int idx = 0; idx < 1000000; idx++)
    dh_insert(histogram_, g_normal(g_generator), DHIST_SINGLE_THREADED);

  for (int idx = 0; idx < histogram_->num_buckets; idx++) {
    decay(histogram_, histogram_->bucket_list[idx], histogram_->generation);
  }

  acc = 0.0;
  for (int idx = 0; idx < histogram_->num_buckets; idx++) {
    dens = weight(
        histogram_, histogram_->bucket_list[idx], 0,
        &lower_bound, &upper_bound);
    acc += dens * (upper_bound - lower_bound);
  }
  EXPECT_NEAR(1.0, acc, 0.00004);
}

TEST_F(HistogramTest, TotalCount) {
  double count;

  count = 0.0;
  for (int idx = 0; idx < 1000; idx++) {
    count *= (1.0 - alpha_);
    count += 1.0;
    dh_insert(histogram_, g_normal(g_generator), DHIST_SINGLE_THREADED);
  }

  EXPECT_LT(0, histogram_->num_buckets);
  EXPECT_NEAR(count, total_count(histogram_, histogram_->generation), 0.00004);
}

TEST_F(HistogramTest, BucketCountsSumToTotalCount) {
  double count;

  count = 0.0;
  for (int idx = 0; idx < histogram_->num_buckets; idx++) {
    count += histogram_->bucket_list[idx]->count;
  }

  EXPECT_NEAR(count, total_count(histogram_, histogram_->generation), 0.00004);
}

TEST_F(HistogramTest, CheckSplitBucket) {
  struct bucket *cursor;
  std::uniform_int_distribution<int> distribution(0, 3);

  histogram_->root = init_bucket();
  histogram_->root->height = 1;

  for (int idx = 0; idx < 900; idx++) {
    cursor = histogram_->root;
    for (int scramble = 0; scramble < 90; scramble++) {
      int choice = distribution(g_generator);
      switch (choice) {
      case 0:
      case 1:
        if (cursor->children[choice])
          cursor = cursor->children[choice];
        break;
      default:
        if (cursor->parent)
          cursor = cursor->parent;
        break;
      }
    }

    split_bucket(histogram_, cursor);
    assert_invariant(histogram_->root);
  }
  //print_tree(histogram_->root);
  assert_invariant(histogram_->root);
}

TEST_F(HistogramTest, CheckDeleteBucket) {
  struct bucket *cursor;
  std::uniform_int_distribution<int> distribution(0, 3);
  int iterations = 100;
  int count = 99;

  while (iterations--) {
    histogram_->root = init_bucket();
    histogram_->root->height = 1;

    for (int idx = 0; idx < count; idx++) {
      cursor = histogram_->root;
      for (int scramble = 0; scramble < 90; scramble++) {
        int choice = distribution(g_generator);
        switch (choice) {
        case 0:
        case 1:
          if (cursor->children[choice]) {
            cursor = cursor->children[choice];
          }
          break;
        default:
          if (cursor->parent) {
            cursor = cursor->parent;
          }
          break;
        }
      }

      split_bucket(histogram_, cursor);
      assert_invariant(histogram_->root);
    }
    assert_invariant(histogram_->root);

    for (int idx = 1; idx < count; idx++) {
      cursor = histogram_->root;
      for (int scramble = 0; scramble < count; scramble++) {
        int choice = distribution(g_generator);
        switch (choice) {
        case 0:
        case 1:
          if (cursor->children[choice])
            cursor = cursor->children[choice];
          break;
        default:
          if (cursor->parent)
            cursor = cursor->parent;
          break;
        }
      }

      delete_bucket(histogram_, cursor);
      assert_invariant(histogram_->root);
    }
  }
}

