#ifndef CMathCorrelation_H
#define CMathCorrelation_H

#include <vector>
#include <cassert>
#include <cmath>

namespace CMathCorrelation {

using Values = std::vector<double>;

double mean(const Values &xv) {
  int nv = xv.size();

  if (nv == 0)
    return 0.0;

  double sum = 0.0;

  for (const auto &x : xv)
    sum += x;

  return sum/nv;
}

double calc(const Values &xv, const Values &yv) {
  assert(xv.size() == yv.size());

  double xmean = mean(xv);
  double ymean = mean(yv);

  int nxy = xv.size();

  if (nxy == 0)
    return 0.0;

  double sum_dxy = 0.0;
  double sum_dxx = 0.0;
  double sum_dyy = 0.0;

  for (int i = 0; i < nxy; ++i) {
    double dx = xv[i] - xmean;
    double dy = yv[i] - ymean;

    double dxy = dx*dy;

    sum_dxy += dxy;

    sum_dxx += dx*dx;
    sum_dyy += dy*dy;
  }

  double p = sum_dxy/(sqrt(sum_dxx)*sqrt(sum_dyy));

  return p;
}

}

#endif
