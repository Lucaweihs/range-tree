//
// Created by Luca Weihs on 11/3/16.
//

#ifndef RANGETREE_RANGETREE_H
#define RANGETREE_RANGETREE_H

#include <vector>
#include <iostream>
#include <sstream>

namespace RangeTree {

    /**
     * A class that represents a multi-dimensional point
     * with some value. This point can have a multiplicity
     * represented by its count.
     */
    template<class T>
    class Point {
    private:
        std::vector<double> vec;
        T val;
        int multiplicity;

    public:
        Point() : multiplicity(0) {}

        Point(std::vector<double> vec, T val): val(val), vec(vec), multiplicity(1) {}

        std::vector<double> asVector() const {
            return vec;
        }

        unsigned long dim() const {
            return vec.size();
        }

        int count() const {
            return multiplicity;
        }

        void increaseCountBy(int n) {
            if (n < 0) {
                throw std::logic_error("Can't increase by a negative amount");
            }
            multiplicity += n;
        }

        T value() const {
            return val;
        }

        double operator[](int index) const {
            if(index < 0 || index >= dim()) {
                throw std::out_of_range("[] access index for point is out of range.");
            }
            return vec[index];
        }

        double operator==(Point<T> p) const {
            return vec == p.vec && multiplicity == p.multiplicity && val == p.val;
        }

        double operator!=(Point<T> p) const {
            return !((*this) == p);
        }

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

    template <class T>
    class PointOrdering {
    private:
        int compareStartIndex;

    public:
        PointOrdering(int compareStartIndex): compareStartIndex(compareStartIndex) {
            if (compareStartIndex < 0) {
                throw new std::logic_error("Cannot have comparison start index <0.");
            }
        }

        static bool equals(const Point<T>& p1, const Point<T>& p2) {
            return p1.asVector() == p2.asVector();
        }

        int getCompareStartIndex() const {
            return compareStartIndex;
        }

        bool less(const Point<T>& p1, const Point<T>& p2) const {
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

        bool lessOrEq(const Point<T>& p1, const Point<T>& p2) const {
            return less(p1, p2) || equals(p1, p2);
        }

        bool greater(const Point<T>& p1, const Point<T>& p2) const {
            return less(p2, p1);
        }

        bool greaterOrEq(const Point<T>& p1, const Point<T>& p2) const {
            return greater(p1, p2) || equals(p1,p2);
        }

        bool operator()(const Point<T>& p1, const Point<T>& p2) const {
            return this->less(p1, p2);
        }
    };

    template <class T>
    class RangeTreeNode {
    private:
        std::shared_ptr<RangeTreeNode<T>> left;
        std::shared_ptr<RangeTreeNode<T>> right;
        std::shared_ptr<RangeTreeNode<T>> treeOnNextDim;
        Point<T> point;
        std::vector<Point<T>> allPoints;
        bool isLeaf;
        int pointCountSum;
        PointOrdering<T> pointOrdering;

        std::shared_ptr<RangeTreeNode<T>> sortedPointsToBinaryTree(std::vector<Point<T>> sortedUniquePoints,
                                                                   int first, int last, int compareStartInd) {
            if (first == last) {
                return std::shared_ptr<RangeTreeNode<T>>(
                        new RangeTreeNode(sortedUniquePoints[first], compareStartInd));
            } else {
                int mid = (first + last) / 2;
                auto left = sortedPointsToBinaryTree(sortedUniquePoints, first, mid, compareStartInd);
                auto right = sortedPointsToBinaryTree(sortedUniquePoints, mid + 1, last, compareStartInd);

                auto allPointsLeft = (*left).getAllPoints();
                auto allPointsRight = (*right).getAllPoints();

                allPointsLeft.insert(allPointsLeft.end(), allPointsRight.begin(), allPointsRight.end());

                return std::shared_ptr<RangeTreeNode<T>>(
                        new RangeTreeNode(left, right,
                                          allPointsLeft,
                                          sortedUniquePoints[mid],
                                          compareStartInd));
            }
        }

    public:
        RangeTreeNode(std::vector<Point<T>> allPoints, int compareStartInd): pointOrdering(compareStartInd) {
            if (allPoints.size() == 0) {
                throw std::range_error("Range tree requires input vector of points to not be empty.");
            }

            std::sort(allPoints.begin(), allPoints.end(), pointOrdering);
            std::vector<Point<T>> sortedUniquePoints = {allPoints[0]};
            int k = 0;
            for (int i = 1; i < allPoints.size(); i++) {
                if (pointOrdering.equals(sortedUniquePoints[k], allPoints[i])) {
                    if (sortedUniquePoints[k].value() != allPoints[i].value()) {
                        throw std::logic_error("Input points have same position but different values");
                    }
                    sortedUniquePoints[k].increaseCountBy(allPoints[i].count());
                } else {
                    sortedUniquePoints.push_back(allPoints[i]);
                    k++;
                }
            }

            int numLeafNodes = sortedUniquePoints.size();
            int mid = (numLeafNodes - 1) / 2;
            point = sortedUniquePoints[mid];

            if (numLeafNodes == 1) {
                isLeaf = true;
                pointCountSum = point.count();
            } else {
                left = sortedPointsToBinaryTree(sortedUniquePoints, 0, mid, compareStartInd);
                right = sortedPointsToBinaryTree(sortedUniquePoints, mid + 1, numLeafNodes - 1, compareStartInd);
                pointCountSum = left->totalPoints() + right->totalPoints();
                if (compareStartInd + 1 != point.dim()) {
                    treeOnNextDim = std::shared_ptr<RangeTreeNode>(
                            new RangeTreeNode(sortedUniquePoints, compareStartInd + 1));
                }
                isLeaf = false;
            }
        }

        RangeTreeNode(std::shared_ptr<RangeTreeNode<T>> left,
                      std::shared_ptr<RangeTreeNode<T>> right,
                      std::vector<Point<T>> allPoints,
                      Point<T> comparePoint,
                      int compareStartInd) :
                left(left), right(right), allPoints(allPoints),
                point(comparePoint), isLeaf(false),
                pointOrdering(compareStartInd) {
            if (compareStartInd + 1 != point.dim()) {
                treeOnNextDim = std::shared_ptr<RangeTreeNode<T>>(new RangeTreeNode<T>(allPoints, compareStartInd + 1));
            }
            pointCountSum = (*left).totalPoints() + (*right).totalPoints();
        }

        RangeTreeNode(Point<T> pointAtLeaf, int compareStartInd) :
                point(pointAtLeaf), isLeaf(true), pointCountSum(pointAtLeaf.count()), pointOrdering(compareStartInd) {
            allPoints.push_back(point);
        }

        int totalPoints() const {
            return pointCountSum;
        }

        std::vector<Point<T>> getAllPoints() const {
            return allPoints;
        }

        bool pointInRange(const Point<T>& point,
                          const std::vector<double>& lower,
                          const std::vector<double>& upper,
                          const std::vector<bool>& withLower,
                          const std::vector<bool>& withUpper) const {
            for (int i = 0; i < point.dim(); i++) {
                if (point[i] < lower[i] ||
                    (!withLower[i] && point[i] == lower[i])) {
                    return false;
                }
                if (point[i] > upper[i] ||
                    (!withUpper[i] && point[i] == upper[i])) {
                    return false;
                }
            }
            return true;
        }

        int countInRange(const std::vector<double>& lower,
                         const std::vector<double>& upper,
                         const std::vector<bool>& withLower,
                         const std::vector<bool>& withUpper) const {
            if (isLeaf) {
                if (pointInRange(point, lower, upper, withLower, withUpper)) {
                    return totalPoints();
                } else {
                    return 0;
                }
            }
            int compareInd = pointOrdering.getCompareStartIndex();

            if ((point[compareInd] > upper[compareInd]) ||
                    (point[compareInd] == upper[compareInd] && !withUpper[compareInd])) {
                return (*left).countInRange(lower, upper, withLower, withUpper);
            }
            if ((point[compareInd] < lower[compareInd]) ||
                    (point[compareInd] == lower[compareInd] && !withLower[compareInd])) {
                return (*right).countInRange(lower, upper, withLower, withUpper);
            }

            std::vector<std::shared_ptr<RangeTreeNode<T>>> canonicalNodes = {};

            if ((*left).isLeaf) {
                canonicalNodes.push_back(left);
            } else {
                (*left).leftCanonicalNodes(lower, withLower, canonicalNodes);
            }

            if ((*right).isLeaf) {
                canonicalNodes.push_back(right);
            } else {
                (*right).rightCanonicalNodes(upper, withUpper, canonicalNodes);
            }

            int numPointsInRange = 0;
            for (int i = 0; i < canonicalNodes.size(); i++) {
                std::shared_ptr<RangeTreeNode<T>> node = canonicalNodes[i];
                if ((*node).isLeaf) {
                    if (pointInRange((*node).point, lower, upper, withLower, withUpper)) {
                        numPointsInRange += (*node).totalPoints();
                    }
                } else if (compareInd + 1 == point.dim()) {
                    numPointsInRange += (*node).totalPoints();
                } else {
                    numPointsInRange += node->treeOnNextDim->countInRange(lower, upper, withLower, withUpper);
                }
            }

            return numPointsInRange;
        }

        void leftCanonicalNodes(const std::vector<double>& lower,
                                const std::vector<bool>& withLower,
                                std::vector<std::shared_ptr<RangeTreeNode<T>>>& nodes) {
            if (isLeaf) {
                throw std::logic_error("Should never have a leaf deciding if its canonical.");
            }
            int compareInd = pointOrdering.getCompareStartIndex();
            int totalPoints = 0;
            if (lower[compareInd] < point[compareInd] ||
                    (lower[compareInd] == point[compareInd] && withLower[compareInd])) {
                nodes.push_back(right);
                if ((*left).isLeaf) {
                    nodes.push_back(left);
                } else {
                    (*left).leftCanonicalNodes(lower, withLower, nodes);
                }
            } else {
                if ((*right).isLeaf) {
                    nodes.push_back(right);
                } else {
                    (*right).leftCanonicalNodes(lower, withLower, nodes);
                }
            }
        }

        void rightCanonicalNodes(const std::vector<double>& upper,
                                 const std::vector<bool>& withUpper,
                                 std::vector<std::shared_ptr<RangeTreeNode<T>>>& nodes) {
            if (isLeaf) {
                throw std::logic_error("Should never have a leaf deciding if its canonical.");
            }
            int compareInd = pointOrdering.getCompareStartIndex();
            int totalPoints = 0;
            if (upper[compareInd] > point[compareInd] ||
                    (upper[compareInd] == point[compareInd] && withUpper[compareInd])) {
                nodes.push_back(left);
                if ((*right).isLeaf) {
                    nodes.push_back(right);
                } else {
                    (*right).rightCanonicalNodes(upper, withUpper, nodes);
                }
            } else {
                if ((*left).isLeaf) {
                    nodes.push_back(left);
                } else {
                    (*left).rightCanonicalNodes(upper, withUpper, nodes);
                }
            }
        }

        void print(int numIndents) {
            for (int i = 0; i < numIndents; i++) { std::cout << "\t"; }
            if (isLeaf) {
                point.print(true);
            } else {
                point.print(false);
                left->print(numIndents + 1);
                right->print(numIndents + 1);
            }
        }
    };

    template <class T>
    class RangeTree {
    private:
        std::shared_ptr<RangeTreeNode<T>> root;

    public:
        RangeTree(std::vector<Point<T>> points) {
            root = std::shared_ptr<RangeTreeNode<T>>(new RangeTreeNode<T>(points, 0));
        }

        int countInRange(const std::vector<double>& lower,
                         const std::vector<double>& upper,
                         const std::vector<bool>& withLower,
                         const std::vector<bool>& withUpper) {
            return (*root).countInRange(lower, upper, withLower, withUpper);
        }

        void print() {
            root->print(0);
        }
    };

    template <class T>
    class NaiveRangeCounter {
    private:
        std::vector<Point<T>> points;

        static bool pointInRange(const Point<T>& point,
                          const std::vector<double>& lower,
                          const std::vector<double>& upper,
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
        NaiveRangeCounter(std::vector<Point<T>> points): points(points) {}

        int countInRange(const std::vector<double>& lower,
                         const std::vector<double>& upper,
                         const std::vector<bool>& withLower,
                         const std::vector<bool>& withUpper) {
            int count = 0;
            for (int i = 0; i < points.size(); i++) {
                if (pointInRange(points[i], lower, upper, withLower, withUpper)) {
                    count += points[i].count();
                }
            }
            return count;
        }
    };

} // namespace

#endif //RANGETREE_RANGETREE_H
