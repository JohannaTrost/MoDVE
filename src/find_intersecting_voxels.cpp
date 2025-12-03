#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

//' Find all the voxels that a segment intersects in 3D space
//'
//' The voxels are found using a 3D version of the Bresenham algorithm.
//'
//' @param coords_begin a numeric vector of length 3 with the x-y-z coordinates
//' of the start point of the segment
//' @param coords_end a numeric vector of length 3 with the x-y-z coordinates
//' of the end point of the segment
//'
//' @returns a list of length 3 integer vectors, the 3D coordinates of all voxels
//' intersected by the segment
//'
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List find_intersecting_voxels(NumericVector coords_begin, NumericVector coords_end) {

  // Find which cell coordinates the points belong to (i.e. ceiling)
  IntegerVector xyz_begin, xyz_end;
  for (auto& el : coords_begin) xyz_begin.push_back(std::ceil(el));
  for (auto& el : coords_end) xyz_end.push_back(std::ceil(el));

  // Bresendham 3D algorithm
  // code adapted from GeeksForGeeks contributor ishankhandelwals
  List intersctdVoxels = List::create(xyz_begin);
  int x1 = xyz_begin[0];
  int y1 = xyz_begin[1];
  int z1 = xyz_begin[2];
  int x2 = xyz_end[0];
  int y2 = xyz_end[1];
  int z2 = xyz_end[2];

  int dx = abs(x2 - x1);
  int dy = abs(y2 - y1);
  int dz = abs(z2 - z1);
  int x_incr = x2 > x1 ? 1 : -1;
  int y_incr = y2 > y1 ? 1 : -1;
  int z_incr = z2 > z1 ? 1 : -1;

  // Driving axis is X-axis"
  if (dx >= dy && dx >= dz) {
    int p1 = 2 * dy - dx;
    int p2 = 2 * dz - dx;
    while (x1 != x2) {
      x1 += x_incr;
      if (p1 >= 0) {
        y1 += y_incr;
        p1 -= 2 * dx;
      }
      if (p2 >= 0) {
        z1 += z_incr;
        p2 -= 2 * dx;
      }
      p1 += 2 * dy;
      p2 += 2 * dz;
      IntegerVector coords{x1, y1, z1};
      intersctdVoxels.push_back(coords);
    }
  }
  // Driving axis is Y-axis
  else if (dy >= dx && dy >= dz) {
    int p1 = 2 * dx - dy;
    int p2 = 2 * dz - dy;
    while (y1 != y2) {
      y1 += y_incr;
      if (p1 >= 0) {
        x1 += x_incr;
        p1 -= 2 * dy;
      }
      if (p2 >= 0) {
        z1 += z_incr;
        p2 -= 2 * dy;
      }
      p1 += 2 * dx;
      p2 += 2 * dz;
      IntegerVector coords{x1, y1, z1};
      intersctdVoxels.push_back(coords);
    }
  }
  // Driving axis is Z-axis
  else {
    int p1 = 2 * dy - dz;
    int p2 = 2 * dx - dz;
    while (z1 != z2) {
      z1 += z_incr;
      if (p1 >= 0) {
        y1 += y_incr;
        p1 -= 2 * dz;
      }
      if (p2 >= 0) {
        x1 += x_incr;
        p2 -= 2 * dz;
      }
      p1 += 2 * dy;
      p2 += 2 * dx;
      IntegerVector coords{x1, y1, z1};
      intersctdVoxels.push_back(coords);
    }
  }

  return intersctdVoxels;
}

