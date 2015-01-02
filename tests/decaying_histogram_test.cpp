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


#include "decaying_histogram.h"
#include "gtest/gtest.h"

TEST(TestCase1, Test1) {
    struct bucket buckets[3];
    init_bucket(&buckets[0], NULL, NULL);
    init_bucket(&buckets[1], &buckets[0], &buckets[2]);
    init_bucket(&buckets[2], NULL, &buckets[2]);
    for (int idx = 0; idx < 3; idx++) {
        EXPECT_EQ(0.0, buckets[idx].alpha_mu);
        EXPECT_EQ(0.0, buckets[idx].alpha_count);
        EXPECT_EQ(0.0, buckets[idx].count);
        EXPECT_EQ(0.0, buckets[idx].mu);
        EXPECT_EQ(0.0, buckets[idx].min);
        EXPECT_EQ(0.0, buckets[idx].max);
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

