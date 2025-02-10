#include "utilities.h"

#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RSnapshotOptions.hxx>

using namespace ROOT::VecOps;

template <typename T>
RVec<float> DiffNext(RVec<T> v) {
    RVec<float> diff;
    for (int i=0; i<v.size()-1; i++) diff.push_back(v[i+1] - v[i]);
    diff.push_back(0);
    return diff;
}
template <typename T>
RVec<float> DiffPrev(RVec<T> v) {
    RVec<float> diff;
    diff.push_back(0);
    for (int i=1; i<v.size(); i++) diff.push_back(v[i] - v[i-1]);
    return diff;
}
template <typename T>
RVec<T> Acc(RVec<T> v) {
    RVec<T> acc;
    T sum = 0;
    for (T x : v) acc.push_back(sum += x);
    return acc; 
}



template <typename T>
float Cov(RVec<T> u, RVec<T> v) {
    return Mean(u*v) - Mean(u) * Mean(v);
}

// template <typename T>
// std::tuple<float, float, float> LinReg(RVec<T> x, RVec<T> y) {
//     float esp_x = Mean(x);
//     float esp_y = Mean(y);
//     float espsq_x = Mean(x*x);
//     float espsq_y = Mean(y*y);
//     float esp_xy = Mean(x*y);
//     float cov = Cov(x, y);

//     float m = (esp_xy - esp_x * esp_y) / (espsq_x - esp_x * esp_x);
//     float q = esp_y - m * esp_x;

//     return {m, q};
// } 

