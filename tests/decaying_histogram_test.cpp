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

TEST_F(HistogramTest, FindBucket) {
  int num_buckets;
  double dens, lower_bound, upper_bound;
  struct bucket *bucket;

  for (num_buckets = 1; num_buckets <= histogram_->max_num_buckets;
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

    for (int idx = 0; idx < num_buckets; idx++) {
      bucket = find_bucket(histogram_, idx - 0.25, DHIST_SINGLE_THREADED);
      ASSERT_EQ(&histogram_->bucket_list[idx], bucket)
          << "idx: " << idx << " num_buckets: " << num_buckets << endl;
    }
  }

  // If we only have one bucket, everything goes into that bucket.
  histogram_->num_buckets = 1;
  histogram_->bucket_list[0].below = NULL;
  histogram_->bucket_list[0].above = NULL;
  EXPECT_EQ(
      &histogram_->bucket_list[0],
      find_bucket(histogram_, -1000, DHIST_SINGLE_THREADED));
}

TEST_F(HistogramTest, TargetBoundary) {
  int num_buckets;
  double dens, lower_bound, upper_bound;
  struct bucket *bucket;
  bool target_boundary_seen;

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
}

TEST_F(HistogramTest, PositiveDensity) {
  double observation;

  for (int idx = 0; idx < 1000; idx++) {
    observation = g_normal(g_generator);
    add_observation(histogram_, observation, DHIST_SINGLE_THREADED);
    //if (histogram_->num_buckets > 1)
    //  FAIL() << histogram_->bucket_list[0].lower_bound << " | " << histogram_->bucket_list[0].upper_bound
    //         << " ; " << histogram_->bucket_list[1].lower_bound << " | " << histogram_->bucket_list[1].upper_bound << endl;
  }

  for (int idx = 0; idx < histogram_->num_buckets; idx++) {
    ASSERT_GT(
        density(histogram_, &histogram_->bucket_list[idx], 0, NULL, NULL),
        0.0);
  }
}

TEST_F(HistogramTest, DensitySumsToOne) {
  double observation, acc, dens, lower_bound, upper_bound;

  for (int idx = 0; idx < 1000000; idx++)
    add_observation(histogram_, g_normal(g_generator), DHIST_SINGLE_THREADED);

  for (int idx = 0; idx < histogram_->num_buckets; idx++) {
    decay(histogram_, &histogram_->bucket_list[idx], histogram_->generation);
  }

  acc = 0.0;
  for (int idx = 0; idx < histogram_->num_buckets; idx++) {
    dens = density(
        histogram_, &histogram_->bucket_list[idx], 0,
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
    add_observation(histogram_, g_normal(g_generator), DHIST_SINGLE_THREADED);
  }

  EXPECT_LT(0, histogram_->num_buckets);
  EXPECT_NEAR(count, total_count(histogram_, histogram_->generation), 0.00004);
}

TEST_F(HistogramTest, BucketCountsSumToTotalCount) {
  double count;

  count = 0.0;
  full_refresh(histogram_, DHIST_SINGLE_THREADED);
  for (int idx = 0; idx < histogram_->num_buckets; idx++) {
    count += histogram_->bucket_list[idx].count;
  }

  EXPECT_NEAR(count, total_count(histogram_, histogram_->generation), 0.00004);
}

TEST_F(HistogramTest, OddBug) {
  struct decaying_histogram *histogram;
  histogram = new decaying_histogram;
  init_decaying_histogram(histogram, 50, 0.000001);
  histogram->generation = 1998485954;
  histogram->num_buckets = 47;

  for (int idx = 0; idx < histogram->num_buckets - 1; idx++) {
    histogram->bucket_list[idx].above = &histogram->bucket_list[idx + 1];
    histogram->bucket_list[idx + 1].below = &histogram->bucket_list[idx];
  }
  histogram->bucket_list[0].below = NULL;
  histogram->bucket_list[histogram->num_buckets - 1].below = NULL;

  histogram->bucket_list[0].count = 6961.3413249679315;
  histogram->bucket_list[0].mu = 9.0585741371445749;
  histogram->bucket_list[0].update_generation = 1998485954;

  histogram->bucket_list[1].count = 10519.276825604944;
  histogram->bucket_list[1].mu = 9.1476029277719419;
  histogram->bucket_list[1].update_generation = 1998485954;

  histogram->bucket_list[2].count = 14082.687659465486;
  histogram->bucket_list[2].mu = 9.2030289543231003;
  histogram->bucket_list[2].update_generation = 1998485954;

  histogram->bucket_list[3].count = 13772.257107927997;
  histogram->bucket_list[3].mu = 9.2437798680184571;
  histogram->bucket_list[3].update_generation = 1998485954;

  histogram->bucket_list[4].count = 15679.291501800875;
  histogram->bucket_list[4].mu = 9.2747806360000347;
  histogram->bucket_list[4].update_generation = 1998485954;

  histogram->bucket_list[5].count = 16277.035393519884;
  histogram->bucket_list[5].mu = 9.302418345043181;
  histogram->bucket_list[5].update_generation = 1998485954;

  histogram->bucket_list[6].count = 15991.270860627799;
  histogram->bucket_list[6].mu = 9.3268128675525617;
  histogram->bucket_list[6].update_generation = 1998485954;

  histogram->bucket_list[7].count = 17431.803581745786;
  histogram->bucket_list[7].mu = 9.3476386594087941;
  histogram->bucket_list[7].update_generation = 1998485954;

  histogram->bucket_list[8].count = 16616.661529967387;
  histogram->bucket_list[8].mu = 9.3679510449588879;
  histogram->bucket_list[8].update_generation = 1998485954;

  histogram->bucket_list[9].count = 19357.041846887576;
  histogram->bucket_list[9].mu = 9.3879324816191652;
  histogram->bucket_list[9].update_generation = 1998485954;

  histogram->bucket_list[10].count = 19933.23635818695;
  histogram->bucket_list[10].mu = 9.4073165147007654;
  histogram->bucket_list[10].update_generation = 1998485954;

  histogram->bucket_list[11].count = 17226.322162945609;
  histogram->bucket_list[11].mu = 9.4259209733903067;
  histogram->bucket_list[11].update_generation = 1998485954;

  histogram->bucket_list[12].count = 15308.779225360091;
  histogram->bucket_list[12].mu = 9.4417783634730021;
  histogram->bucket_list[12].update_generation = 1998485954;

  histogram->bucket_list[13].count = 17047.252602391902;
  histogram->bucket_list[13].mu = 9.4567087149321818;
  histogram->bucket_list[13].update_generation = 1998485954;

  histogram->bucket_list[14].count = 16892.87213652255;
  histogram->bucket_list[14].mu = 9.4719764978239418;
  histogram->bucket_list[14].update_generation = 1998485954;

  histogram->bucket_list[15].count = 13516.610651886003;
  histogram->bucket_list[15].mu = 9.4871109016366297;
  histogram->bucket_list[15].update_generation = 1998485954;

  histogram->bucket_list[16].count = 18509.590650480943;
  histogram->bucket_list[16].mu = 9.5021941213488574;
  histogram->bucket_list[16].update_generation = 1998485954;

  histogram->bucket_list[17].count = 18967.582805759837;
  histogram->bucket_list[17].mu = 9.5197098849084547;
  histogram->bucket_list[17].update_generation = 1998485954;

  histogram->bucket_list[18].count = 22352.516500903967;
  histogram->bucket_list[18].mu = 9.5406914933287084;
  histogram->bucket_list[18].update_generation = 1998485954;

  histogram->bucket_list[19].count = 21457.542632850687;
  histogram->bucket_list[19].mu = 9.5640584535576103;
  histogram->bucket_list[19].update_generation = 1998485954;

  histogram->bucket_list[20].count = 18473.260540208597;
  histogram->bucket_list[20].mu = 9.5872903770552096;
  histogram->bucket_list[20].update_generation = 1998485954;

  histogram->bucket_list[21].count = 20144.672577434478;
  histogram->bucket_list[21].mu = 9.6128554154635495;
  histogram->bucket_list[21].update_generation = 1998485954;

  histogram->bucket_list[22].count = 20270.892559220072;
  histogram->bucket_list[22].mu = 9.6422806651710733;
  histogram->bucket_list[22].update_generation = 1998485954;

  histogram->bucket_list[23].count = 20862.947598877476;
  histogram->bucket_list[23].mu = 9.6750160715512123;
  histogram->bucket_list[23].update_generation = 1998485954;

  histogram->bucket_list[24].count = 23058.340054507382;
  histogram->bucket_list[24].mu = 9.7114605534858391;
  histogram->bucket_list[24].update_generation = 1998485954;

  histogram->bucket_list[25].count = 25973.050304038239;
  histogram->bucket_list[25].mu = 9.7512649026562972;
  histogram->bucket_list[25].update_generation = 1998485954;

  histogram->bucket_list[26].count = 29517.040832196002;
  histogram->bucket_list[26].mu = 9.7938758612504149;
  histogram->bucket_list[26].update_generation = 1998485954;

  histogram->bucket_list[27].count = 33323.19328332704;
  histogram->bucket_list[27].mu = 9.8407810634997173;
  histogram->bucket_list[27].update_generation = 1998485954;

  histogram->bucket_list[28].count = 35003.821654246021;
  histogram->bucket_list[28].mu = 9.8915305603600263;
  histogram->bucket_list[28].update_generation = 1998485954;

  histogram->bucket_list[29].count = 34447.767953428178;
  histogram->bucket_list[29].mu = 9.9486715979972633;
  histogram->bucket_list[29].update_generation = 1998485954;

  histogram->bucket_list[30].count = 30006.4990400235;
  histogram->bucket_list[30].mu = 10.013032298260068;
  histogram->bucket_list[30].update_generation = 1998485954;

  histogram->bucket_list[31].count = 23554.156646552394;
  histogram->bucket_list[31].mu = 10.086572273730033;
  histogram->bucket_list[31].update_generation = 1998485954;

  histogram->bucket_list[32].count = 19411.366850554998;
  histogram->bucket_list[32].mu = 10.166330167423325;
  histogram->bucket_list[32].update_generation = 1998485954;

  histogram->bucket_list[33].count = 19055.093343485169;
  histogram->bucket_list[33].mu = 10.24228329779112;
  histogram->bucket_list[33].update_generation = 1998485954;

  histogram->bucket_list[34].count = 20623.442034099193;
  histogram->bucket_list[34].mu = 10.312455032284648;
  histogram->bucket_list[34].update_generation = 1998485954;

  histogram->bucket_list[35].count = 23226.744372891324;
  histogram->bucket_list[35].mu = 10.379877432094174;
  histogram->bucket_list[35].update_generation = 1998485954;

  histogram->bucket_list[36].count = 26868.068941156045;
  histogram->bucket_list[36].mu = 10.448248274196965;
  histogram->bucket_list[36].update_generation = 1998485954;

  histogram->bucket_list[37].count = 30296.6838122962;
  histogram->bucket_list[37].mu = 10.521956250880839;
  histogram->bucket_list[37].update_generation = 1998485954;

  histogram->bucket_list[38].count = 31786.386339147281;
  histogram->bucket_list[38].mu = 10.606920415278219;
  histogram->bucket_list[38].update_generation = 1998485954;

  histogram->bucket_list[39].count = 29275.837688718573;
  histogram->bucket_list[39].mu = 10.708818262801133;
  histogram->bucket_list[39].update_generation = 1998485954;

  histogram->bucket_list[40].count = 25989.799592068361;
  histogram->bucket_list[40].mu = 10.829352052285435;
  histogram->bucket_list[40].update_generation = 1998485954;

  histogram->bucket_list[41].count = 26587.812668665061;
  histogram->bucket_list[41].mu = 10.957072489293742;
  histogram->bucket_list[41].update_generation = 1998485954;

  histogram->bucket_list[42].count = 27106.796043975275;
  histogram->bucket_list[42].mu = 11.091263869618814;
  histogram->bucket_list[42].update_generation = 1998485954;

  histogram->bucket_list[43].count = 24662.875960578283;
  histogram->bucket_list[43].mu = 11.252438478177837;
  histogram->bucket_list[43].update_generation = 1998485954;

  histogram->bucket_list[44].count = 19149.160078622746;
  histogram->bucket_list[44].mu = 11.499108094372476;
  histogram->bucket_list[44].update_generation = 1998485954;

  histogram->bucket_list[45].count = 8733.1823486090398;
  histogram->bucket_list[45].mu = 11.919118174005325;
  histogram->bucket_list[45].update_generation = 1998485954;

  histogram->bucket_list[46].count = 24690.13349128791;
  histogram->bucket_list[46].mu = 15.957688766645205;
  histogram->bucket_list[46].update_generation = 1998485954;

  add_observation(histogram, 10.604553229079038, DHIST_MULTI_THREADED);
  clean_decaying_histogram(histogram);
  delete histogram;
}

