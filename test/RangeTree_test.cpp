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

TEST(range_tree_test, returns_correct_for_one_dim)
{
    std::vector<double> values = {1.0,1.0,2.0,2.0,3.0,3.0,3.0,4.0,4.0,5.0};
    std::vector<RT::Point<int>> points = {};

    auto f = [](double a) { std::vector<double> b = {a}; return b;};
    for (int i = 0; i < values.size(); i++) {
        RT::Point<int> a(f(values[i]), 0);
        points.push_back(a);
    }

    RT::RangeTree<int> rtree(points);

    auto g = [](bool a) { std::vector<bool> b = {a}; return b;};
    EXPECT_EQ(rtree.countInRange(f(-12.0), f(30.0), g(true), g(true)), 10);
    EXPECT_EQ(rtree.countInRange(f(2.0), f(4.0), g(true), g(true)), 7);
    EXPECT_EQ(rtree.countInRange(f(2.0), f(4.0), g(false), g(true)), 5);
    EXPECT_EQ(rtree.countInRange(f(2.0), f(4.0), g(true), g(false)), 5);
    EXPECT_EQ(rtree.countInRange(f(2.0), f(4.0), g(false), g(false)), 3);
    EXPECT_EQ(rtree.countInRange(f(3.0), f(3.0), g(false), g(false)), 0);
    EXPECT_EQ(rtree.countInRange(f(3.0), f(4.0), g(false), g(false)), 0);
}

TEST(range_tree_test, returns_correct_for_two_dim)
{
    std::vector<double> v1 = {3.0, 3.0, 3.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 1.0, 3.1};
    std::vector<double> v2 = {1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 3.2};
    std::vector<RT::Point<int>> points = {};

    auto f = [](double a, double b) { std::vector<double> c = {a, b}; return c;};
    for (int i = 0; i < v1.size(); i++) {
        RT::Point<int> a(f(v1[i], v2[i]), 0);
        points.push_back(a);
    }

    RT::RangeTree<int> rtree(points);

    auto g = [](bool a, bool b) { std::vector<bool> c = {a, b}; return c;};
    // Selecting 2-dim regions with boundary
    EXPECT_EQ(rtree.countInRange(f(0.0, 0.0), f(4.0, 4.0),
                                 g(true, true), g(true, true)), 11);
    EXPECT_EQ(rtree.countInRange(f(0.0, 0.0), f(1.0, 1.0),
                                 g(true, true), g(true, true)), 2);
    EXPECT_EQ(rtree.countInRange(f(0.0, 0.0), f(1.0, 3.0),
                                 g(true, true), g(true, true)), 4);
    EXPECT_EQ(rtree.countInRange(f(1.9, 1.9), f(2.1, 2.1),
                                 g(true, true), g(true, true)), 1);
    EXPECT_EQ(rtree.countInRange(f(-1000.0, -1000.0), f(1000.0, 0.9999999),
                                 g(true, true), g(true, true)), 0);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
                                 g(true, true), g(true, true)), 10);

    // Selecting 2-dim regions without (some) boundaries
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
                                 g(false, true), g(true, true)), 6);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
                                 g(true, false), g(true, true)), 6);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
                                 g(true, true), g(false, true)), 7);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
                                 g(true, true), g(true, false)), 7);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
                                 g(false, false), g(true, true)), 4);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
                                 g(true, true), g(false, false)), 5);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
                                 g(true, false), g(true, false)), 3);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
                                 g(false, true), g(false, true)), 3);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(3.0, 3.0),
                                 g(false, false), g(false, false)), 1);

    // Selecting 0/1-dim regions with boundary
    EXPECT_EQ(rtree.countInRange(f(1.0, 0.0), f(1.0, 2.0),
                                 g(true, true), g(true, true)), 3);
    EXPECT_EQ(rtree.countInRange(f(0.0, 2.0), f(3.0, 2.0),
                                 g(true, true), g(true, true)), 3);
    EXPECT_EQ(rtree.countInRange(f(0.0, 2.0), f(3.0, 2.0),
                                 g(true, true), g(true, true)), 3);
    EXPECT_EQ(rtree.countInRange(f(1.3, 3.0), f(2.3, 3.0),
                                 g(true, true), g(true, true)), 1);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 1.0),
                                 g(true, true), g(true, true)), 2);
    EXPECT_EQ(rtree.countInRange(f(0.0, 0.0), f(0.0, 1.0),
                                 g(true, true), g(true, true)), 0);

    // Selecting 0/1-dim regions without (some) boundaries
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 2.0),
                                 g(false, true), g(true, true)), 0);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 2.0),
                                 g(true, false), g(true, true)), 1);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 2.0),
                                 g(true, true), g(false, true)), 0);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 2.0),
                                 g(true, true), g(true, false)), 2);

    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 1.0),
                                 g(false, true), g(true, true)), 0);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 1.0),
                                 g(true, false), g(true, true)), 0);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 1.0),
                                 g(true, true), g(false, true)), 0);
    EXPECT_EQ(rtree.countInRange(f(1.0, 1.0), f(1.0, 1.0),
                                 g(true, true), g(true, false)), 0);
}

TEST(range_tree_test, returns_correct_for_three_dim)
{
    auto f = [](double a, double b, double c) { std::vector<double> d = {a, b, c}; return d;};
    std::vector<RT::Point<int>> points = {};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                RT::Point<int> a(f(i, j, k), 0);
                points.push_back(a);
                if (i == 1 && j == 1 && k == 1) {
                    points.push_back(a);
                }
            }
        }
    }
    RT::RangeTree<int> rtree(points);
    RT::NaiveRangeCounter<int> nrc(points);

    auto g = [](bool a, bool b, bool c) { std::vector<bool> d = {a, b, c}; return d;};
    std::vector<std::vector<bool>> boolVectors = {};
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                std::vector<bool> a = {};
                a.push_back(i == 1);
                a.push_back(j == 1);
                a.push_back(k == 1);
                boolVectors.push_back(a);
            }
        }
    }

    // Selecting 3 dim region with all possible boundaries
    auto lower = {0.0, 0.0, 0.0};
    auto upper = {2.0, 2.0, 2.0};

    for (int i = 0; i < boolVectors.size(); i++) {
        for (int j = 0; j < boolVectors.size(); j++) {
            auto withLower = boolVectors[i];
            auto withUpper = boolVectors[j];

            EXPECT_EQ(
                    nrc.countInRange(lower, upper, withLower, withUpper),
                    rtree.countInRange(lower, upper, withLower, withUpper)
            );
        }
    }

    // Selecting 0,1,2,3 dim region with all possible boundaries
    auto lower3d = {0.0, 0.0, 0.0};
    auto upper3d = {2.0, 2.0, 2.0};

    auto lower2d = {1.0, 0.0, 0.0};
    auto upper2d = {1.0, 2.0, 2.0};

    auto lower1d = {1.0, 1.0, 0.0};
    auto upper1d = {1.0, 1.0, 2.0};

    auto lower0d = {1.0, 1.0, 1.0};
    auto upper0d = {1.0, 1.0, 1.0};

    for (int i = 0; i < boolVectors.size(); i++) {
        for (int j = 0; j < boolVectors.size(); j++) {
            auto withLower = boolVectors[i];
            auto withUpper = boolVectors[j];

            EXPECT_EQ(nrc.countInRange(lower3d, upper3d, withLower, withUpper),
                      rtree.countInRange(lower3d, upper3d, withLower, withUpper));
            EXPECT_EQ(nrc.countInRange(lower2d, upper2d, withLower, withUpper),
                      rtree.countInRange(lower2d, upper2d, withLower, withUpper));
            EXPECT_EQ(nrc.countInRange(lower1d, upper1d, withLower, withUpper),
                      rtree.countInRange(lower1d, upper1d, withLower, withUpper));
            EXPECT_EQ(nrc.countInRange(lower0d, upper0d, withLower, withUpper),
                      rtree.countInRange(lower0d, upper0d, withLower, withUpper));
        }
    }
}