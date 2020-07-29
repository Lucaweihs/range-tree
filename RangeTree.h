/*
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

/**
* Implements the RangeTree datastructure.
*
* See the documentation of the class RangeTree for more details regarding the purpose
* of RangeTree's.
*
* This file contains two main abstractions meant for use:
* 1. The RangeTree class.
* 2. A Point class which captures the idea of a d-dimensional euclidean point.
*/

#ifndef RANGETREE_H
#define RANGETREE_H

#include <vector>
#include <iostream>
#include <sstream>
#include <numeric>
#include <type_traits>
#include <deque>
#include <cmath>
#include <algorithm>

namespace RangeTree {

    /**
    * A point in euclidean space.
    *
    * A class that represents a multi-dimensional euclidean point
    * with some associated value. We allow for each point to have an
    * associated value so that some more information can be stored with
    * each point. Points can also have a multiplicity/count, this corresponds
    * to having several duplicates of the same point.
    */
    template<typename T, class S>
    class Point {
        static_assert(std::is_arithmetic<T>::value, "Type T must be numeric");
    private:
        std::vector<T> vec;
        S val;
        int multiplicity;

    public:
        /**
        * Constructs an empty point.
        *
        * Creates a point in 0 dimensional euclidean space. This constructor
        * is provided only to make certain edge cases easier to handle.
        */
        Point() : multiplicity(0) {}

        /**
        * Constructs a point.
        *
        * Creates a point with its position in euclidean space defined by vec,
        * value defined by val, and a multiplicity/count of 1.
        *
        * @param vec the position in euclidean space.
        * @param val the value associated with the point.
        */
        Point(const std::vector<T>& vec, const S& val): val(val), vec(vec), multiplicity(1) {}

        /**
        * Constructs a point.
        *
        * Copies a point.
        *
        * @param vec the position in euclidean space.
        * @param val the value associated with the point.
        */
        Point(const Point<T,S>& p): val(p.val), vec(p.vec), multiplicity(p.count()) {}


        /**
        * Euclidean position of the point.
        *
        * @return the euclidean position of the point as a std::vector.
        */
        const std::vector<T>& asVector() const {
            return vec;
        }

        /**
        * The point's ambient dimension.
        *
        * @return the dimension of the space in which the point lives. I.e. a point of the
        *         form (1,2,3) lives in dimension 3.
        */
        unsigned long dim() const {
            return vec.size();
        }

        /**
        * The point's count/multiplicity.
        *
        * @return returns the count/multiplicity.
        */
        int count() const {
            return multiplicity;
        }

        /**
        * Increase the point's count/multiplicity.
        *
        * @param n amount to increase by.
        */
        void increaseCountBy(int n) {
            if (n < 0) {
                throw std::logic_error("Can't increase by a negative amount");
            }
            multiplicity += n;
        }

        /**
        * The point's value.
        *
        * @return the value stored in the point.
        */
        S value() const {
            return val;
        }

        /**
        * Index a point.
        *
        * Get the ith coordinate value of the point. I.e. if a point is of the form (4,5,6),
        * then its 0th coordinate value is 4 while its 2nd is 6.
        *
        * @param index the coordinate to index.
        * @return the coordinate value.
        */
        T operator[](int index) const {
            if(index < 0 || index >= dim()) {
                throw std::out_of_range("[] access index for point is out of range.");
            }
            return vec[index];
        }

        /**
        * Check for equality.
        *
        * Two points are considered equal if they are in the same spot, have the same
        * multiplicity/count, and store the same value.
        *
        * @param p some other point
        * @return true if \p equals the current point, otherwise false.
        */
        bool operator==(const Point<T,S>& p) const {
            return vec == p.vec && multiplicity == p.multiplicity && val == p.val;
        }

        /**
        * Check for inequality.
        *
        * The opposite of ==.
        *
        * @param p some other point.
        * @return false if \p equals the current point, otherwise true.
        */
        bool operator!=(const Point<T,S>& p) const {
            return !((*this) == p);
        }

        /**
        * Prints the point to standard out.
        *
        * As an example, a point with euclidean location (3,4,5) and with a
        * multiplicity/count of 4 will be printed as
        *
        * (3, 4, 5) : 4
        *
        * @param withCount whether or not to display the points count/multiplicity.
        */
        void print(bool withCount=true) const {
            std::cout << "(";
            for (int i = 0; i < dim() - 1; i++) {
                std::cout << (*this)[i] << ", ";
            }
            if (withCount) {
                std::cout << (*this)[dim() - 1] << ") : " << count() << std::endl;
            } else {
                std::cout << (*this)[dim() - 1] << ") : " << std::endl;
            }
        }
    };

    /**
    * A class that totally orders Point<T,S>'s in euclidean space.
    *
    * A total order of Points is required in the RangeTree. This is an implementation
    * detail that can be ignored. Given a start index \compareStartIndex, this class
    * orders points so that, for two points p_1 = (p_{11}, p_{12},...,p_{1n}) and
    * p_2 = (p_{21}, p_{22},...,p_{2n}) we have that p_1 < p_2 if and only if
    *
    * (p_{1\compareStartInd},...,p_{1n}) < (p_{2\compareStartInd},...,p_{2n})
    *
    * using the usual lexicographic order, or
    *
    * (p_{1\compareStartInd},...,p_{1n}) == (p_{2\compareStartInd},...,p_{2n}) and
    * (p_{11},...,p_{1(\compareStartInd-1)}) < (p_{21},...,p_{2(\compareStartInd-1)})
    *
    * again using the usual lexicographic order.
    */
    template <typename T, class S>
    class PointOrdering {
        static_assert(std::is_arithmetic<T>::value, "Type T must be numeric");
    private:
        int compareStartIndex;

    public:
        PointOrdering(int compareStartIndex): compareStartIndex(compareStartIndex) {
            if (compareStartIndex < 0) {
                throw new std::logic_error("Cannot have comparison start index <0.");
            }
        }

        static bool equals(const Point<T,S>& p1, const Point<T,S>& p2) {
            return p1.asVector() == p2.asVector();
        }

        int getCompareStartIndex() const {
            return compareStartIndex;
        }

        bool less(const Point<T,S>& p1, const Point<T,S>& p2) const {
            if (p1.dim() != p2.dim()) {
                throw std::logic_error("Points are incomparable (differing dims).");
            }
            if (compareStartIndex >= p1.dim()) {
                throw std::logic_error("Cannot compare points, compare start index >= point dimension.");
            }
            for (int i = compareStartIndex; i < p1.dim(); i++) {
                if (p1[i] < p2[i]) {
                    return true;
                } else if (p1[i] > p2[i]) {
                    return false;
                }
            }
            for (int i = 0; i < compareStartIndex; i++) {
                if (p1[i] < p2[i]) {
                    return true;
                } else if (p1[i] > p2[i]) {
                    return false;
                }
            }
            return false;
        }

        bool lessOrEq(const Point<T,S>& p1, const Point<T,S>& p2) const {
            return less(p1, p2) || equals(p1, p2);
        }

        bool greater(const Point<T,S>& p1, const Point<T,S>& p2) const {
            return less(p2, p1);
        }

        bool greaterOrEq(const Point<T,S>& p1, const Point<T,S>& p2) const {
            return greater(p1, p2) || equals(p1,p2);
        }

        bool operator()(const Point<T,S>& p1, const Point<T,S>& p2) const {
            return this->less(p1, p2);
        }
    };

    /**
    * A matrix that keeps a collection of points sorted on each coordinate
    */
    template<typename T, class S>
    class SortedPointMatrix {
        static_assert(std::is_arithmetic<T>::value, "Type T must be numeric");
    private:
        std::vector<Point<T,S>* > pointsSortedByCurrentDim;
        std::deque<std::vector<int> > redirectionTable;
        int currentDim;
        int dim;
        static const int MAX_POINTS_BEFORE_SWITCH = 1000;

        std::vector<int> sortOrder(const std::vector<Point<T,S>* >& points, int onDim) {
            std::vector<int> order(points.size());
            for (int i = 0; i < points.size(); i++) { order[i] = i; }
            PointOrdering<T,S> pointOrdering(onDim);
            std::sort(order.begin(), order.end(),
                      [pointOrdering, points](int i, int j) {
                          return pointOrdering.less(*(points[i]), *(points[j]));
                      });
            return order;
        }

        void sort(std::vector<Point<T,S>* >& points, int onDim) {
            PointOrdering<T,S> pointOrdering(onDim);
            std::sort(points.begin(), points.end(),
                      [pointOrdering](Point<T,S>* pt0, Point<T,S>* pt1) {
                          return pointOrdering.less(*(pt0), *(pt1));
                      });
        }

        void rearrangeGivenOrder(std::vector<Point<T,S>* >& points,
                                  const std::vector<int>& order) {
            std::vector<Point<T,S>* > tmp = points;
            for (int i = 0; i < points.size(); i++) {
                points[i] = tmp[order[i]];
            }
        }

        SortedPointMatrix(const std::vector<Point<T,S>* >& pointsSortedByCurrentDim,
                          const std::deque<std::vector<int> >& redirectionTable,
                          int currentDim, int dim) : pointsSortedByCurrentDim(pointsSortedByCurrentDim), redirectionTable(redirectionTable),
                                        currentDim(currentDim), dim(dim) {}

    public:
        /**
        * Constructs a sorted point matrix
        */
        SortedPointMatrix(std::vector<Point<T,S>* >& points): currentDim(0) {
            if (points.size() == 0) {
                throw std::range_error("Cannot construct a SortedPointMatrix with 0 points.");
            } else {
                dim = points[0]->dim();
                for (int i = 1; i < points.size(); i++) {
                    if (points[i]->dim() != dim) {
                        throw std::logic_error("Input points to SortedPointMatrix must all"
                                                       " have the same dimension.");
                    }
                }

                int sortDimension = (points.size() > MAX_POINTS_BEFORE_SWITCH) ? dim - 1 : 0;
                PointOrdering<T,S> pointOrdering(sortDimension);
                std::sort(points.begin(), points.end(),
                          [pointOrdering](Point<T,S>* p1, Point<T,S>* p2) {
                              return pointOrdering.less(*p1, *p2);
                          });

                pointsSortedByCurrentDim.push_back(points[0]);
                int k = 0;
                for (int i = 1; i < points.size(); i++) {
                    if (pointOrdering.equals(*(pointsSortedByCurrentDim[k]), *points[i])) {
                        if (pointsSortedByCurrentDim[k]->value() != points[i]->value()) {
                            throw std::logic_error("Input points have same position but different values");
                        }
                        pointsSortedByCurrentDim[k]->increaseCountBy(points[i]->count());
                    } else {
                        pointsSortedByCurrentDim.push_back(points[i]);
                        k++;
                    }
                }

                if (pointsSortedByCurrentDim.size() > MAX_POINTS_BEFORE_SWITCH) {
                    for (int i = dim - 2; i >= currentDim; i--) {
                        std::vector<int> order = sortOrder(pointsSortedByCurrentDim, i);
                        redirectionTable.push_front(order);
                        rearrangeGivenOrder(pointsSortedByCurrentDim, order);
                    }
                }
            }
        }

        void moveToNextDimension() {
            if (currentDim == dim - 1) {
                throw std::logic_error("Already at max dimension, cannot move to next.");
            }
            currentDim++;
            if (pointsSortedByCurrentDim.size() > MAX_POINTS_BEFORE_SWITCH) {
                std::vector<Point<T,S>* > tmp = pointsSortedByCurrentDim;
                for (int i = 0; i < pointsSortedByCurrentDim.size(); i++) {
                    pointsSortedByCurrentDim[redirectionTable[0][i]] = tmp[i];
                }
                redirectionTable.pop_front();
            } else {
                sort(pointsSortedByCurrentDim, currentDim);
            }
        }

        Point<T,S>* getMidPoint() {
            int mid = (numUniquePoints() - 1) / 2;
            return pointsSortedByCurrentDim[mid];
        }

        int numUniquePoints() {
            return pointsSortedByCurrentDim.size();
        }

        int getCurrentDim() {
            return currentDim;
        }

        std::vector<Point<T,S>* > getSortedPointsAtCurrentDim() {
            return pointsSortedByCurrentDim;
        }

        /**
        * Constructs two sorted point matrices after splitting on the current midpoint
        */
        std::pair<SortedPointMatrix, SortedPointMatrix> splitOnMid() {
            int n = numUniquePoints();
            if (n == 1) {
                throw std::logic_error("Cannot split on mid when there is only one point.");
            }

            int mid = (n - 1) / 2;
            std::vector<Point<T, S> *> sortedPointsLeft(mid + 1), sortedPointsRight(n - mid - 1);
            for (int i = 0; i < mid + 1; i++) {
                sortedPointsLeft[i] = pointsSortedByCurrentDim[i];
            }
            for (int i = mid + 1; i < n; i++) {
                sortedPointsRight[i - mid - 1] = pointsSortedByCurrentDim[i];
            }

            if (n <= MAX_POINTS_BEFORE_SWITCH) {
                std::deque<std::vector<int> > redirectionTableLeft, redirectionTableRight;
                return std::make_pair(
                        SortedPointMatrix(sortedPointsLeft, redirectionTableLeft, currentDim, dim),
                        SortedPointMatrix(sortedPointsRight, redirectionTableRight, currentDim, dim));
            } else {
                std::vector<bool> onLeft(n);
                for (int i = 0; i < n; i++) {
                    onLeft[i] = i <= mid;
                }

                std::deque<std::vector<int> > redirectionTableLeft(redirectionTable.size(),
                                                                   std::vector<int>(mid + 1));
                std::deque<std::vector<int> > redirectionTableRight(redirectionTable.size(),
                                                                    std::vector<int>(n - mid - 1));
                for (int i = 0; i < redirectionTable.size(); i++) {
                    std::vector<bool> lastOnLeft = onLeft;

                    for (int j = 0; j < numUniquePoints(); j++) {
                        onLeft[redirectionTable[i][j]] = lastOnLeft[j];
                    }

                    std::vector<int> newRedirect(numUniquePoints());
                    int kLeft = 0, kRight = 0;
                    for (int j = 0; j < numUniquePoints(); j++) {
                        if (onLeft[j]) {
                            newRedirect[j] = kLeft;
                            kLeft++;
                        } else {
                            newRedirect[j] = kRight;
                            kRight++;
                        }
                    }

                    kLeft = 0, kRight = 0;
                    for (int j = 0; j < numUniquePoints(); j++) {
                        if (lastOnLeft[j]) {
                            redirectionTableLeft[i][kLeft] = newRedirect[redirectionTable[i][j]];
                            kLeft++;
                        } else {
                            redirectionTableRight[i][kRight] = newRedirect[redirectionTable[i][j]];
                            kRight++;
                        }
                    }
                }
                return std::make_pair(
                        SortedPointMatrix(sortedPointsLeft, redirectionTableLeft, currentDim, dim),
                        SortedPointMatrix(sortedPointsRight, redirectionTableRight, currentDim, dim));
            }
        }
    };

    /**
    * A class representing a single node in a RangeTree. These should not be
    * constructed directly, instead use the RangeTree class.
    */
    template <typename T, class S>
    class RangeTreeNode {
        static_assert(std::is_arithmetic<T>::value, "Type T must be numeric");
    private:
        std::shared_ptr<RangeTreeNode<T,S> > left; /**< Contains points <= the comparison point **/
        std::shared_ptr<RangeTreeNode<T,S> > right; /**< Contains points > the comparison point **/
        std::shared_ptr<RangeTreeNode<T,S> > treeOnNextDim; /**< Tree on the next dimension **/
        Point<T,S>* point; /**< The comparison point **/
        bool isLeaf; /**< Whether or not the point is a leaf **/
        int pointCountSum; /**< Total number of points, counting multiplicities, at leaves of the tree **/
        PointOrdering<T,S> pointOrdering; /**< Helper to totally order input points **/

        // For fractional cascading
        std::vector<T> pointsLastDimSorted;
        std::vector<Point<T,S>* > allPointsSorted;
        std::vector<int> pointerToGeqLeft;
        std::vector<int> pointerToLeqLeft;
        std::vector<int> pointerToGeqRight;
        std::vector<int> pointerToLeqRight;
        std::vector<int> cumuCountPoints;

    public:
        /**
        * Construct a range tree structure from points.
        *
        * Creates a range tree structure on the input collection \allPoints using the lexicographic order
        * starting at \compareStartInd.
        *
        * @param uniquePoints a collection of points.
        * @return a range tree structure
        */
        RangeTreeNode(SortedPointMatrix<T,S>& spm,
                      bool onLeftEdge = true,
                      bool onRightEdge = true): pointOrdering(spm.getCurrentDim()) {
            point = spm.getMidPoint();

            if (spm.numUniquePoints() == 1) {
                isLeaf = true;
                pointCountSum = point->count();
                pointsLastDimSorted.push_back((*point)[point->dim() - 1]);
                if (spm.getCurrentDim() == point->dim() - 2) {
                    spm.moveToNextDimension();
                }
            } else {
                auto spmPair = spm.splitOnMid();
                left = std::shared_ptr<RangeTreeNode<T,S> >(
                        new RangeTreeNode<T,S>(spmPair.first, onLeftEdge, false));
                right = std::shared_ptr<RangeTreeNode<T,S> >(
                        new RangeTreeNode<T,S>(spmPair.second, false, onRightEdge));
                pointCountSum = left->totalPoints() + right->totalPoints();

                int dim = point->dim();
                if (spm.getCurrentDim() + 2 == dim) {
                    spm.moveToNextDimension();

                    allPointsSorted = spm.getSortedPointsAtCurrentDim();
                    cumuCountPoints.push_back(0);
                    for (int i = 0; i < allPointsSorted.size(); i++) {
                        pointsLastDimSorted.push_back((*allPointsSorted[i])[dim - 1]);
                        cumuCountPoints.push_back(cumuCountPoints.back() + allPointsSorted[i]->count());
                    }
                    const auto& leftSorted = left->pointsLastDimSorted;
                    const auto& rightSorted = right->pointsLastDimSorted;

                    pointerToGeqLeft = createGeqPointers(pointsLastDimSorted, leftSorted);
                    pointerToGeqRight = createGeqPointers(pointsLastDimSorted, rightSorted);
                    pointerToLeqLeft = createLeqPointers(pointsLastDimSorted, leftSorted);
                    pointerToLeqRight = createLeqPointers(pointsLastDimSorted, rightSorted);
                }  else if (!onLeftEdge && !onRightEdge && spm.getCurrentDim() + 1 != point->dim()) {
                        spm.moveToNextDimension();
                        treeOnNextDim = std::shared_ptr<RangeTreeNode>(new RangeTreeNode(spm));
                }
                isLeaf = false;
            }
        }

        std::vector<int> createGeqPointers(const std::vector<T>& vec,
                                           const std::vector<T>& subVec) {
            std::vector<int> grePointers(vec.size());
            int k = 0;
            for (int i = 0; i < vec.size(); i++) {
                while (k < subVec.size() && subVec[k] < vec[i]) {
                    k++;
                }
                grePointers[i] = k;
            }
            return grePointers;
        }

        std::vector<int> createLeqPointers(const std::vector<T>& vec,
                                           const std::vector<T>& subVec) {
            std::vector<int> leqPointers(vec.size());
            int k = subVec.size() - 1;
            for (int i = vec.size() - 1; i >= 0; i--) {
                while (k >= 0 && subVec[k] > vec[i]) {
                    k--;
                }
                leqPointers[i] = k;
            }
            return leqPointers;
        }

        int binarySearchFirstGeq(T needle, int left, int right) const {
            if (left == right) {
                if (needle <= pointsLastDimSorted[left]) {
                    return left;
                } else {
                    return left + 1;
                }
            }
            int mid = (left + right) / 2;
            if (needle <= pointsLastDimSorted[mid]) {
                return binarySearchFirstGeq(needle, left, mid);
            } else {
                return binarySearchFirstGeq(needle, mid + 1, right);
            }
        }

        int binarySearchFirstLeq(T needle, int left, int right) const {
            if (left == right) {
                if (needle >= pointsLastDimSorted[left]) {
                    return left;
                } else {
                    return left - 1;
                }
            }
            int mid = (left + right + 1) / 2;
            if (needle >= pointsLastDimSorted[mid]) {
                return binarySearchFirstLeq(needle, mid, right);
            } else {
                return binarySearchFirstLeq(needle, left, mid - 1);
            }
        }

        /**
        * Construct a RangeTreeNode representing a leaf.
        *
        * @param pointAtLeaf the point to use at the leaf.
        * @param compareStartInd the index defining the lexicographic ordering.
        * @return
        */
        RangeTreeNode(Point<T,S>* pointAtLeaf, int compareStartInd) :
                point(pointAtLeaf), isLeaf(true), pointCountSum(pointAtLeaf->count()), pointOrdering(compareStartInd) {}

        /**
        * Total count of points at the leaves of the range tree rooted at this node.
        *
        * The total count returned INCLUDES the multiplicities/count of the points at the leaves.
        *
        * @return the total count.
        */
        int totalPoints() const {
            return pointCountSum;
        }

        /**
        * Return all points at the leaves of the range tree rooted at this node.
        * @return all the points.
        */
        std::vector<Point<T,S> > getAllPoints() const {
            if (isLeaf) {
                std::vector<Point<T,S> > vec;
                vec.push_back(*point);
                return vec;
            }
            auto allPointsLeft = left->getAllPoints();
            auto allPointsRight = right->getAllPoints();

            allPointsLeft.insert(allPointsLeft.end(), allPointsRight.begin(), allPointsRight.end());
            return allPointsLeft;
        }

        /**
        * Check if point is in a euclidean box.
        *
        * Determines whether or not a point is in a euclidean box. That is, if p_1 = (p_{11},...,p_{1n}) is an
        * n-dimensional point. Then this function returns true if, for all 1 <= i <= n we have
        *
        * lower[i] <= p_{1i} <= upper[i] if withLower[i] == true and withUpper[i] == true, or
        * lower[i] < p_{1i} <= upper[i] if withLower[i] == false and withUpper[i] == true, or
        * lower[i] <= p_{1i} < upper[i] if withLower[i] == true and withUpper[i] == false, or
        * lower[i] < p_{1i} < upper[i] if withLower[i] == false and withUpper[i] == false.
        *
        * @param point the point to check.
        * @param lower the lower points of the rectangle.
        * @param upper the upper bounds of the rectangle.
        * @param withLower whether to use strict (<) or not strict (<=) inequalities at certain coordiantes of p_1
        *                  for the lower bounds.
        * @param withUpper as for \withLower but for the upper bounds.
        * @return true if the point is in the rectangle, false otherwise.
        */
        bool pointInRange(const Point<T,S>& point,
                          const std::vector<T>& lower,
                          const std::vector<T>& upper) const {
            for (int i = 0; i < point.dim(); i++) {
                if (point[i] < lower[i]) {
                    return false;
                }
                if (point[i] > upper[i]) {
                    return false;
                }
            }
            return true;
        }

        /**
        * Count the number of points at leaves of tree rooted at the current node that are within the given bounds.
        *
        * @param lower see the pointInRange(...) function.
        * @param upper
        * @return the count.
        */
        int countInRange(const std::vector<T>& lower,
                         const std::vector<T>& upper) const {
            if (isLeaf) {
                if (pointInRange(*point, lower, upper)) {
                    return totalPoints();
                } else {
                    return 0;
                }
            }
            int compareInd = pointOrdering.getCompareStartIndex();

            if ((*point)[compareInd] > upper[compareInd]) {
                return left->countInRange(lower, upper);
            }
            if ((*point)[compareInd] < lower[compareInd]) {
                return right->countInRange(lower, upper);
            }

            int dim = point->dim();
            if (compareInd + 2 == dim) {
                int n = pointsLastDimSorted.size();
                int geqInd = binarySearchFirstGeq(lower.back(), 0, n - 1);
                int leqInd = binarySearchFirstLeq(upper.back(), 0, n - 1);

                if (geqInd > leqInd) {
                    return 0;
                }
                std::vector<RangeTreeNode<T, S>* > nodes;
                std::vector<std::pair<int,int> > inds;
                left->leftFractionalCascade(lower,
                                            pointerToGeqLeft[geqInd],
                                            pointerToLeqLeft[leqInd],
                                            nodes,
                                            inds);
                right->rightFractionalCascade(upper,
                                              pointerToGeqRight[geqInd],
                                              pointerToLeqRight[leqInd],
                                              nodes,
                                              inds);
                int sum = 0;
                for (int i = 0; i < nodes.size(); i++) {
                    if (nodes[i]->isLeaf) {
                        sum += nodes[i]->totalPoints();
                    } else {
                        sum += nodes[i]->cumuCountPoints[inds[i].second + 1] -
                                nodes[i]->cumuCountPoints[inds[i].first];
                    }
                }
                return sum;
            } else {
                std::vector<std::shared_ptr<RangeTreeNode<T, S> > > canonicalNodes;

                if (left->isLeaf) {
                    canonicalNodes.push_back(left);
                } else {
                    left->leftCanonicalNodes(lower, canonicalNodes);
                }

                if (right->isLeaf) {
                    canonicalNodes.push_back(right);
                } else {
                    right->rightCanonicalNodes(upper, canonicalNodes);
                }

                int numPointsInRange = 0;
                for (int i = 0; i < canonicalNodes.size(); i++) {
                    std::shared_ptr<RangeTreeNode<T, S> > node = canonicalNodes[i];
                    if (node->isLeaf) {
                        if (pointInRange(*(node->point), lower, upper)) {
                            numPointsInRange += node->totalPoints();
                        }
                    } else if (compareInd + 1 == point->dim()) {
                        numPointsInRange += node->totalPoints();
                    } else {
                        numPointsInRange += node->treeOnNextDim->countInRange(lower, upper);
                    }
                }
                return numPointsInRange;
            }
        }

        /**
        * Return the points at leaves of tree rooted at the current node that are within the given bounds.
        *
        * @param lower see the pointInRange(...) function.
        * @param upper
        * @return a std::vector of the Points.
        */
        std::vector<Point<T,S> > pointsInRange(const std::vector<T>& lower,
                                               const std::vector<T>& upper) const {
            std::vector<Point<T,S> > pointsToReturn = {};
            if (isLeaf) {
                if (pointInRange(*point, lower, upper)) {
                    pointsToReturn.push_back(*point);
                }
                return pointsToReturn;
            }
            int compareInd = pointOrdering.getCompareStartIndex();

            if ((*point)[compareInd] > upper[compareInd]) {
                return left->pointsInRange(lower, upper);
            }
            if ((*point)[compareInd] < lower[compareInd]) {
                return right->pointsInRange(lower, upper);
            }

            int dim = point->dim();
            if (compareInd + 2 == dim) {
                int n = pointsLastDimSorted.size();
                int geqInd = binarySearchFirstGeq(lower.back(), 0, n - 1);
                int leqInd = binarySearchFirstLeq(upper.back(), 0, n - 1);

                if (geqInd > leqInd) {
                    return pointsToReturn;
                }
                std::vector<RangeTreeNode<T, S>* > nodes;
                std::vector<std::pair<int,int> > inds;
                left->leftFractionalCascade(lower,
                                            pointerToGeqLeft[geqInd],
                                            pointerToLeqLeft[leqInd],
                                            nodes,
                                            inds);
                right->rightFractionalCascade(upper,
                                              pointerToGeqRight[geqInd],
                                              pointerToLeqRight[leqInd],
                                              nodes,
                                              inds);
                for (int i = 0; i < nodes.size(); i++) {
                    if (nodes[i]->isLeaf) {
                        pointsToReturn.push_back(*(nodes[i]->point));
                    } else {
                        for (int j = inds[i].first; j <= inds[i].second; j++) {
                            pointsToReturn.push_back(*(nodes[i]->allPointsSorted[j]));
                        }
                    }
                }
                return pointsToReturn;
            } else {
                std::vector<std::shared_ptr<RangeTreeNode<T, S> > > canonicalNodes = {};

                if (left->isLeaf) {
                    canonicalNodes.push_back(left);
                } else {
                    left->leftCanonicalNodes(lower, canonicalNodes);
                }

                if (right->isLeaf) {
                    canonicalNodes.push_back(right);
                } else {
                    right->rightCanonicalNodes(upper, canonicalNodes);
                }

                for (int i = 0; i < canonicalNodes.size(); i++) {
                    std::shared_ptr<RangeTreeNode<T, S> > node = canonicalNodes[i];
                    if (node->isLeaf) {
                        if (pointInRange(*(node->point), lower, upper)) {
                            pointsToReturn.push_back(*(node->point));
                        }
                    } else if (compareInd + 1 == point->dim()) {
                        auto allPointsAtNode = node->getAllPoints();
                        pointsToReturn.insert(pointsToReturn.end(), allPointsAtNode.begin(), allPointsAtNode.end());
                    } else {
                        auto allPointsAtNode = node->treeOnNextDim->pointsInRange(lower, upper);
                        pointsToReturn.insert(pointsToReturn.end(), allPointsAtNode.begin(), allPointsAtNode.end());
                    }
                }
                return pointsToReturn;
            }
        }

        void leftFractionalCascade(const std::vector<T>& lower,
                                  int geqInd,
                                  int leqInd,
                                  std::vector<RangeTreeNode<T,S>* >& nodes,
                                  std::vector<std::pair<int,int> >& inds) {
            if (leqInd < geqInd) {
                return;
            }

            int compareInd = point->dim() - 2;

            if (lower[compareInd] <= (*point)[compareInd]) {
                if (isLeaf) {
                    nodes.push_back(this);
                    inds.push_back(std::pair<int,int>(0,0));
                    return;
                }

                int geqIndRight = pointerToGeqRight[geqInd];
                int leqIndRight = pointerToLeqRight[leqInd];
                if (leqIndRight >= geqIndRight) {
                    nodes.push_back(right.get());
                    if (right->isLeaf) {
                        inds.push_back(std::pair<int,int>(0,0));
                    } else {
                        inds.push_back(std::pair<int,int>(geqIndRight,leqIndRight));
                    }
                }

                left->leftFractionalCascade(lower,
                                            pointerToGeqLeft[geqInd],
                                            pointerToLeqLeft[leqInd],
                                            nodes,
                                            inds);
            } else {
                if (isLeaf) {
                    return;
                }
                right->leftFractionalCascade(lower,
                                             pointerToGeqRight[geqInd],
                                             pointerToLeqRight[leqInd],
                                             nodes,
                                             inds);
            }
        }

        void rightFractionalCascade(const std::vector<T>& upper,
                                   int geqInd,
                                   int leqInd,
                                   std::vector<RangeTreeNode<T,S>* >& nodes,
                                   std::vector<std::pair<int,int> >& inds) {
            if (leqInd < geqInd) {
                return;
            }

            int compareInd = point->dim() - 2;

            if ((*point)[compareInd] <= upper[compareInd]) {
                if (isLeaf) {
                    nodes.push_back(this);
                    inds.push_back(std::pair<int,int>(0,0));
                    return;
                }

                int geqIndLeft = pointerToGeqLeft[geqInd];
                int leqIndLeft = pointerToLeqLeft[leqInd];
                if (leqIndLeft >= geqIndLeft) {
                    nodes.push_back(left.get());
                    if (left->isLeaf) {
                        inds.push_back(std::pair<int,int>(0,0));
                    } else {
                        inds.push_back(std::pair<int,int>(geqIndLeft,leqIndLeft));
                    }
                }
                right->rightFractionalCascade(upper,
                                              pointerToGeqRight[geqInd],
                                              pointerToLeqRight[leqInd],
                                              nodes,
                                              inds);
            } else {
                if (isLeaf) {
                    return;
                }
                left->rightFractionalCascade(upper,
                                             pointerToGeqLeft[geqInd],
                                             pointerToLeqLeft[leqInd],
                                             nodes,
                                             inds);
            }
        }

        /**
        * Helper function for countInRange(...).
        * @param lower
        * @param withLower
        * @param nodes
        */
        void leftCanonicalNodes(const std::vector<T>& lower,
                                std::vector<std::shared_ptr<RangeTreeNode<T,S> > >& nodes) {
            if (isLeaf) {
                throw std::logic_error("Should never have a leaf deciding if its canonical.");
            }
            int compareInd = pointOrdering.getCompareStartIndex();
            int totalPoints = 0;
            if (lower[compareInd] <= (*point)[compareInd]) {
                nodes.push_back(right);
                if (left->isLeaf) {
                    nodes.push_back(left);
                } else {
                    left->leftCanonicalNodes(lower, nodes);
                }
            } else {
                if (right->isLeaf) {
                    nodes.push_back(right);
                } else {
                    right->leftCanonicalNodes(lower, nodes);
                }
            }
        }

        /**
        * Helper function for countInRange(...).
        * @param upper
        * @param nodes
        */
        void rightCanonicalNodes(const std::vector<T>& upper,
                                 std::vector<std::shared_ptr<RangeTreeNode<T,S> > >& nodes) {
            if (isLeaf) {
                throw std::logic_error("Should never have a leaf deciding if its canonical.");
            }
            int compareInd = pointOrdering.getCompareStartIndex();
            int totalPoints = 0;
            if (upper[compareInd] >= (*point)[compareInd]) {
                nodes.push_back(left);
                if (right->isLeaf) {
                    nodes.push_back(right);
                } else {
                    right->rightCanonicalNodes(upper, nodes);
                }
            } else {
                if (left->isLeaf) {
                    nodes.push_back(left);
                } else {
                    left->rightCanonicalNodes(upper, nodes);
                }
            }
        }

        /**
        * Print the structure of the tree rooted at the curret node.
        *
        * The printed structure does not reflect any subtrees for other coordinates.
        *
        * @param numIndents the number of indents to use before every line printed.
        */
        void print(int numIndents) {
            for (int i = 0; i < numIndents; i++) { std::cout << "\t"; }
            if (isLeaf) {
                point->print(true);
            } else {
                point->print(false);
                left->print(numIndents + 1);
                right->print(numIndents + 1);
            }
        }
    };

    /**
    * A class facilitating fast orthogonal range queries.
    *
    * A RangeTree allows for 'orthogonal range queries.' That is, given a collection of
    * points P = {p_1, ..., p_n} in euclidean d-dimensional space, a RangeTree can efficiently
    * answer questions of the form
    *
    * "How many points of p are in the box high dimensional rectangle
    * [l_1, u_1] x [l_2, u_2] x ... x [l_d, u_d]
    * where l_1 <= u_1, ..., l_n <= u_n?"
    *
    * It returns the number of such points in worst case
    * O(log(n)^d) time. It can also return the points that are in the rectangle in worst case
    * O(log(n)^d + k) time where k is the number of points that lie in the rectangle.
    *
    * The particular algorithm implemented here is described in Chapter 5 of the book
    *
    * Mark de Berg, Otfried Cheong, Marc van Kreveld, and Mark Overmars. 2008.
    * Computational Geometry: Algorithms and Applications (3rd ed. ed.). TELOS, Santa Clara, CA, USA.
    */
    template <typename T, class S>
    class RangeTree {
        static_assert(std::is_arithmetic<T>::value, "Type T must be numeric");
    private:
        std::shared_ptr<RangeTreeNode<T,S> > root;
        std::vector<std::shared_ptr<Point<T,S> > > savedPoints;
        std::vector<Point<T,S>* > savedPointsRaw;

        std::vector<std::shared_ptr<Point<T,S> > > copyPointsToHeap(const std::vector<Point<T,S> >& points) {
            std::vector<std::shared_ptr<Point<T,S> > > vecOfPointers;
            for (int i = 0; i < points.size(); i++) {
                vecOfPointers.push_back(std::shared_ptr<Point<T,S> >(new Point<T,S>(points[i])));
            }
            return vecOfPointers;
        }

        std::vector<Point<T,S>* > getRawPointers(std::vector<std::shared_ptr<Point<T,S> > >& points) {
            std::vector<Point<T,S>* > vecOfPointers;
            for (int i = 0; i < points.size(); i++) {
                vecOfPointers.push_back(points[i].get());
            }
            return vecOfPointers;
        }

        std::vector<T> getModifiedLower(const std::vector<T>& lower,
                         const std::vector<bool>& withLower) const {
            std::vector<T> newLower = lower;
            for (int i = 0; i < lower.size(); i++) {
                if (std::is_integral<T>::value) {
                    if (!withLower[i]) {
                        newLower[i]++;
                    }
                } else {
                    if (!withLower[i]) {
                        newLower[i] = std::nextafter(newLower[i], std::numeric_limits<T>::max());
                    }
                }
            }
            return newLower;
        }

        std::vector<T> getModifiedUpper(const std::vector<T>& upper,
                                        const std::vector<bool>& withUpper) const {
            std::vector<T> newUpper = upper;
            for (int i = 0; i < upper.size(); i++) {
                if (std::is_integral<T>::value) {
                    if (!withUpper[i]) {
                        newUpper[i]--;
                    }
                } else {
                    if (!withUpper[i]) {
                        newUpper[i] = std::nextafter(newUpper[i], std::numeric_limits<T>::lowest());
                    }
                }
            }
            return newUpper;
        }

    public:
        /**
        * Construct a new RangeTree from input points.
        *
        * Input points may have duplicates but if two points p_1 and p_2 are at the same euclidean position
        * then they are required to have the same value as all duplicate points will be accumulated into a
        * single point with multiplicity/count equal to the sum of the multiplicities/counts of all such duplicates.
        *
        * @param points the points from which to create a RangeTree
        */
        RangeTree(const std::vector<Point<T,S> >& points): savedPoints(copyPointsToHeap(points)),
                                                           savedPointsRaw(getRawPointers(savedPoints)) {
            SortedPointMatrix<T,S> spm(savedPointsRaw);
            root = std::shared_ptr<RangeTreeNode<T,S> >(new RangeTreeNode<T,S>(spm));
        }

        /**
        * The number of points within a high dimensional rectangle.
        *
        * The rectangle is defined by the input parameters. In particular, an n-dimensional point
        * p_1 = (p_{11},...,p_{1n}) is an is in the rectangle if, for all 1 <= i <= n, we have
        *
        * lower[i] <= p_{1i} <= upper[i] if withLower[i] == true and withUpper[i] == true, or
        * lower[i] < p_{1i} <= upper[i] if withLower[i] == false and withUpper[i] == true, or
        * lower[i] <= p_{1i} < upper[i] if withLower[i] == true and withUpper[i] == false, or
        * lower[i] < p_{1i} < upper[i] if withLower[i] == false and withUpper[i] == false.
        *
        * @param lower the lower bounds of the rectangle.
        * @param upper the upper bounds of the rectangle.
        * @param withLower whether to use strict (<) or not strict (<=) inequalities at certain coordiantes of p_1
        *                  for the lower bounds.
        * @param withUpper as for \withLower but for the upper bounds.
        * @return the number of points in the rectangle.
        */
        int countInRange(const std::vector<T>& lower,
                         const std::vector<T>& upper,
                         const std::vector<bool>& withLower,
                         const std::vector<bool>& withUpper) const {
            if (lower.size() != upper.size() || lower.size() != withLower.size() ||
                    lower.size() != withUpper.size()) {
                throw std::logic_error("All vectors inputted to countInRange must have the same length.");
            }
            for (int i = 0; i < lower.size(); i++) {
                if (((!withUpper[i] || !withLower[i]) && lower[i] >= upper[i]) ||
                    lower[i] > upper[i]) {
                    return 0;
                }
            }
            return root->countInRange(getModifiedLower(lower, withLower),
                                      getModifiedUpper(upper, withUpper));
        }

        /**
        * The number of points within a high dimensional rectangle.
        *
        * The rectangle is defined by the input parameters. In particular, an n-dimensional point
        * p_1 = (p_{11},...,p_{1n}) is an is in the rectangle if, for all 1 <= i <= n, we have
        *
        * lower[i] <= p_{1i} <= upper[i]
        *
        * @param lower the lower bounds of the rectangle.
        * @param upper the upper bounds of the rectangle.

        * @return the number of points in the rectangle.
        */
        int countInRange(const std::vector<T>& lower,
                         const std::vector<T>& upper) const {
            if (lower.size() != upper.size()) {
                throw std::logic_error("upper and lower in countInRange must have the same length.");
            }
            return root->countInRange(lower, upper);
        }

        /**
        * Return all points in range.
        *
        * Returns a std::vector of all points in the given rectangle. See \countInRange for how
        * this rectangle is specified. NOTE: these points may not be identical to those points
        * that were given as input to the RangeTree at construction time. This is because
        * duplicate points are merged together with appropriate incrementing of their multiplicity.
        * That is, two points at euclidean position (1,2,3) and multiplicities/counts of 2 and 3
        * respectively will be merged into a single Point with position (1,2,3) and multiplicity 5
        * (recall that all points with the same euclidean position are required to have the same
        * associated value so it is ok to merge in this way).
        *
        * @param lower the lower bounds of the rectangle.
        * @param upper the upper bounds of the rectangle.
        * @param withLower whether to use strict (<) or not strict (<=) inequalities at certain coordiantes of p_1
        *                  for the lower bounds.
        * @param withUpper as for \withLower but for the upper bounds.
        * @return the number of points in the rectangle.
        */
        std::vector<Point<T,S> > pointsInRange(const std::vector<T>& lower,
                                               const std::vector<T>& upper,
                                               const std::vector<bool>& withLower,
                                               const std::vector<bool>& withUpper) const {
            if (lower.size() != upper.size() || lower.size() != withLower.size() ||
                lower.size() != withUpper.size()) {
                throw std::logic_error("All vectors inputted to pointsInRange must have the same length.");
            }
            for (int i = 0; i < lower.size(); i++) {
                if (((!withUpper[i] || !withLower[i]) && lower[i] >= upper[i]) ||
                    lower[i] > upper[i]) {
                    return std::vector<Point<T,S> >();
                }
            }
            return root->pointsInRange(getModifiedLower(lower, withLower),
                                       getModifiedUpper(upper, withUpper));
        }

        void print() const {
            root->print(0);
        }
    };

/**
* A class which is used to naively count the number of points in a given rectangle. This class is used
* for testing an benchmarking, it should not be used in practice.
*/
    template <typename T, class S>
    class NaiveRangeCounter {
        static_assert(std::is_arithmetic<T>::value, "Type T must be numeric");
    private:
        std::vector<Point<T,S> > points;

        static bool pointInRange(const Point<T,S>& point,
                                 const std::vector<T>& lower,
                                 const std::vector<T>& upper,
                                 const std::vector<bool>& withLower,
                                 const std::vector<bool>& withUpper) {
            for (int i = 0; i < point.dim(); i++) {
                if (point[i] < lower[i] ||
                    (point[i] == lower[i] && !withLower[i])) {
                    return false;
                }
                if (point[i] > upper[i] ||
                    (point[i] == upper[i] && !withUpper[i])) {
                    return false;
                }
            }
            return true;
        }

    public:
        NaiveRangeCounter(std::vector<Point<T,S> > points): points(points) {}

        int countInRange(const std::vector<T>& lower,
                         const std::vector<T>& upper,
                         const std::vector<bool>& withLower,
                         const std::vector<bool>& withUpper) const {
            int count = 0;
            for (int i = 0; i < points.size(); i++) {
                if (pointInRange(points[i], lower, upper, withLower, withUpper)) {
                    count += points[i].count();
                }
            }
            return count;
        }

        std::vector<Point<T,S> > pointsInRange(const std::vector<T>& lower,
                                               const std::vector<T>& upper,
                                               const std::vector<bool>& withLower,
                                               const std::vector<bool>& withUpper) const {
            std::vector<Point<T,S> > selectedPoints = {};
            for (int i = 0; i < points.size(); i++) {
                if (pointInRange(points[i], lower, upper, withLower, withUpper)) {
                    selectedPoints.push_back(points[i]);
                }
            }
            return selectedPoints;
        }
    };

} // namespace

#endif //RANGETREE_H
