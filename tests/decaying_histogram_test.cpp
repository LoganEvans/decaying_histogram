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
    destroy_buckets(histogram_);
    delete histogram_;
  }

  struct decaying_histogram *histogram_;
  double alpha_;
  int target_buckets_;
};

TEST_F(HistogramTest, FindBucket) {
  int num_buckets = histogram_->max_num_buckets;
  histogram_->num_buckets = num_buckets;
  for (int idx = 0; idx < num_buckets; idx++) {
    histogram_->bucket_list[idx].lower_bound = idx + 0.5;
    histogram_->bucket_list[idx].upper_bound = idx + 1.5;
    histogram_->bucket_list[idx].count = 1.0;
    histogram_->bucket_list[idx].mu = idx;
  }
  for (int idx = 1; idx < num_buckets - 1; idx++) {
    histogram_->bucket_list[idx].below = &histogram_->bucket_list[idx - 1];
    histogram_->bucket_list[idx].above = &histogram_->bucket_list[idx + 1];
  }
  histogram_->bucket_list[0].below = NULL;
  histogram_->bucket_list[num_buckets - 1].above = NULL;

  for (int idx = 0; idx < num_buckets; idx++) {
    EXPECT_EQ(idx, find_bucket_idx(histogram_, idx + 1.0));
  }

  // If we only have one bucket, everything goes into that bucket.
  histogram_->num_buckets = 1;
  histogram_->bucket_list[0].below = NULL;
  histogram_->bucket_list[0].above = NULL;
  EXPECT_EQ(0, find_bucket_idx(histogram_, -1000));
}

TEST_F(HistogramTest, Bootstrap) {
  EXPECT_EQ(0, 0);
}

TEST_F(HistogramTest, PositiveDensity) {
  double observation;

  for (int idx = 0; idx < 1000; idx++) {
    observation = g_normal(g_generator);
    add_observation(histogram_, observation);
    //if (histogram_->num_buckets > 1)
    //  FAIL() << histogram_->bucket_list[0].lower_bound << " | " << histogram_->bucket_list[0].upper_bound
    //         << " ; " << histogram_->bucket_list[1].lower_bound << " | " << histogram_->bucket_list[1].upper_bound << endl;
  }

  for (int idx = 0; idx < histogram_->num_buckets; idx++) {
    ASSERT_GT(density(histogram_, &histogram_->bucket_list[idx]), 0.0);
  }
}

TEST_F(HistogramTest, DensitySumsToOne) {
  double observation, acc;

  for (int idx = 0; idx < 1000000; idx++)
    add_observation(histogram_, g_normal(g_generator));

  acc = 0.0;
  for (int idx = 0; idx < histogram_->num_buckets; idx++) {
    acc += density(histogram_, &histogram_->bucket_list[idx]) *
        (histogram_->bucket_list[idx].upper_bound -
         histogram_->bucket_list[idx].lower_bound);
  }
  EXPECT_EQ(1.0, acc);
}

TEST_F(HistogramTest, TotalCount) {
  double count;

  count = 0.0;
  for (int idx = 0; idx < 1000; idx++) {
    count *= (1.0 - alpha_);
    count += 1.0;
    add_observation(histogram_, g_normal(g_generator));
  }

  EXPECT_LT(0, histogram_->num_buckets);
  EXPECT_NEAR(count, total_count(histogram_), 0.0000001);

  cout << count << endl;
  cout << histogram_->num_buckets << endl;
  cout << histogram_->max_num_buckets << endl;
  sprint_histogram_new(histogram_, NULL);
}

TEST(HistogramBucketTest, BucketInit) {
    struct bucket buckets[3];
    init_bucket(&buckets[0], NULL, NULL, 0.01);
    init_bucket(&buckets[1], &buckets[0], &buckets[2], 0.01);
    init_bucket(&buckets[2], NULL, &buckets[2], 0.01);
    for (int idx = 0; idx < 3; idx++) {
        EXPECT_EQ(0.01, buckets[idx].alpha);
        EXPECT_EQ(0.0, buckets[idx].count);
        EXPECT_EQ(0.0, buckets[idx].mu);
        EXPECT_EQ(0.0, buckets[idx].lower_bound);
        EXPECT_EQ(0.0, buckets[idx].upper_bound);
        EXPECT_EQ(0, buckets[idx].last_decay_generation);
    }
    EXPECT_EQ(NULL, buckets[0].below);
    EXPECT_EQ(NULL, buckets[0].above);
    EXPECT_EQ(&buckets[0], buckets[1].below);
    EXPECT_EQ(&buckets[2], buckets[1].above);
    EXPECT_EQ(NULL, buckets[2].below);
    EXPECT_EQ(&buckets[2], buckets[2].above);
}

