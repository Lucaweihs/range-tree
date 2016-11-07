/***
 * MIT License
 *
 * Copyright (c) 2016 Luca Weihs
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <vector>
#include "gtest/gtest.h"
#include "RangeTree.h"

namespace RT = RangeTree;

TEST(point_test, has_input_values)
{
    std::vector<double> vals = {1.0, 2.0, -1.0, 3.0};

    RT::Point<int> p(vals, 5);

    EXPECT_EQ(p.count(), 1);
    p.increaseCountBy(1);
    p.increaseCountBy(3);
    EXPECT_EQ(p.count(), 5);
    EXPECT_THROW(p.increaseCountBy(-1), std::logic_error);

	EXPECT_EQ(p[0], 1.0);
    EXPECT_EQ(p[1], 2.0);
    EXPECT_EQ(p[2], -1.0);
    EXPECT_EQ(p[3], 3.0);
    EXPECT_THROW(p[5], std::out_of_range);
    EXPECT_THROW(p[-1], std::out_of_range);

    EXPECT_EQ(p.dim(), 4);

    EXPECT_EQ(p.value(), 5);
}

TEST(point_test, compares_correctly)
{
    RT::PointOrdering<int> ptOrd0(0);
    RT::PointOrdering<int> ptOrd1(1);
    RT::PointOrdering<int> ptOrd2(2);
    RT::PointOrdering<int> ptOrd12(12);

    std::vector<double> vals1 = {1.0, 4.0, -1.0, 4.0};
    std::vector<double> vals2 = {2.2, 3.3, -1.0, 3.0};
    std::vector<double> vals3 = {1.0, 2.0, -1.0, 3.0, 4.0};

    RT::Point<int> p1(vals1, 5);
    RT::Point<int> p1yes(vals1, 5);
    RT::Point<int> p1no(vals1, 6);
    RT::Point<int> p2(vals2, 2);
    RT::Point<int> p3(vals3, 2);

    EXPECT_TRUE(p1 == p1yes);
    EXPECT_TRUE(ptOrd0.equals(p1, p1yes));
    EXPECT_TRUE(ptOrd1.equals(p1, p1yes));
    EXPECT_TRUE(ptOrd2.equals(p1, p1yes));
    EXPECT_FALSE(p1 == p1no);
    EXPECT_TRUE(ptOrd0.equals(p1, p1no));
    EXPECT_TRUE(ptOrd1.equals(p1, p1no));
    EXPECT_TRUE(ptOrd2.equals(p1, p1no));

    EXPECT_TRUE(ptOrd0.less(p1, p2));
    EXPECT_TRUE(ptOrd0.lessOrEq(p1, p2));
    EXPECT_TRUE(ptOrd0.greaterOrEq(p2, p1));
    EXPECT_TRUE(ptOrd0.greater(p2, p1));
    EXPECT_FALSE(ptOrd0.equals(p1, p2));

    EXPECT_TRUE(ptOrd1.greater(p1, p2));
    EXPECT_TRUE(ptOrd1.greaterOrEq(p1, p2));
    EXPECT_TRUE(ptOrd1.lessOrEq(p2, p1));
    EXPECT_TRUE(ptOrd1.less(p2, p1));
    EXPECT_FALSE(ptOrd1.equals(p1, p2));

    EXPECT_TRUE(ptOrd2.greater(p1, p2));
    EXPECT_TRUE(ptOrd2.greaterOrEq(p1, p2));
    EXPECT_TRUE(ptOrd2.lessOrEq(p2, p1));
    EXPECT_TRUE(ptOrd2.less(p2, p1));
    EXPECT_FALSE(ptOrd2.equals(p1, p2));

    EXPECT_THROW(ptOrd12.less(p1, p2), std::logic_error);

    EXPECT_THROW(ptOrd0(p1, p3), std::logic_error);
}