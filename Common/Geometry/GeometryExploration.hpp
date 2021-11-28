#pragma once
#include <Eigen/Dense>
#include <vector>
#include <set>
#include <iostream>
#include <functional>

namespace geometry {
    template <typename T>
    using Vector3 = Eigen::Matrix<T, 3, 1>;
    template <typename T>
    using Vector2 = Eigen::Matrix<T, 2, 1>;
    template <typename T>
    using VectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    template <typename T>
    class VectorCompare {
    public:
        bool operator () (const Vector3<T> &v1, const Vector3<T> &v2) { return v1(1) < v2(1); }
    };

    // TODO: HW5
    // part 2.1 implement nlog(n) method for 2d Pareto front
    // The input is a vector of 2d points
    // The output should be the 2d points on the Pareto front
    // sort your return vector for comparing with 
    // the reference result
    template <typename T>
    bool comp(Vector2<T> a,Vector2<T> b){
        if (a(0)<b(0)){
            return true;
        }else if(a(0)==b(0)&&a(1)<b(1)){
            return true;
        }else{
            return false;
        }
    }
    template <typename T>
    std::vector<Vector2<T>> ParetoFront2D(const std::vector<Vector2<T>> &points) {
        std::vector<Vector2<T>> results = points;
        std::sort(results.begin(),results.end(),comp<T>);
        T curry = results[0](1);
        std::vector<Vector2<T>> finalresults;
        finalresults.clear();
        finalresults.push_back(results[0]);
        for (int i = 1;i<results.size();i++){
            if(results[i](1)<curry){
                finalresults.push_back(results[i]);
                curry = results[i](1);
            }
        }
        return finalresults; // you should change this line as well
    }
    
}
