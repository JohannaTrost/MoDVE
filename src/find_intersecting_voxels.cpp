#include <algorithm>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

//' Find all the voxels that a segment intersects in 3D space
//'
//' @param xyz_begin a numeric vector of length 3 with the x-y-z coordinates
//' of the start point of the segment
//' @param xyz_end a numeric vector of length 3 with the x-y-z coordinates
//' of the end point of the segment
//'
//' @returns a list of length 3 integer vectors, the 3D coordinates of all voxels
//' intersected by the segment
//'
//

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
List find_intersecting_voxels(NumericVector xyz_begin, NumericVector xyz_end) {

  //std::transform(xyz_begin.begin(), xyz_end.cend(), xyz_begin.begin(), std::ceil);
  IntegerVector trunc_begin, trunc_end;
  for (auto& el : xyz_begin) trunc_begin.push_back(std::ceil(el));
  for (auto& el : xyz_end) trunc_end.push_back(std::ceil(el));

  return List::create(trunc_begin, trunc_end);
}

// [[Rcpp::export]]
IntegerVector findInterval2(NumericVector x, NumericVector breaks) {
  IntegerVector out(x.size());

  NumericVector::iterator it, pos;
  IntegerVector::iterator out_it;

  for(it = x.begin(), out_it = out.begin(); it != x.end();
  ++it, ++out_it) {
    pos = std::upper_bound(breaks.begin(), breaks.end(), *it);
    *out_it = std::distance(breaks.begin(), pos);
  }

  return out;
}

