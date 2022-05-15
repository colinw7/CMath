#ifndef CMathCorrelation_H
#define CMathCorrelation_H

#include <vector>
#include <cassert>
#include <cmath>

namespace CMathCorrelation {

using Values = std::vector<double>;

inline double mean(const Values &xv) {
  auto nv = xv.size();
  if (nv == 0) return 0.0;

  double sum = 0.0;

  for (const auto &x : xv)
    sum += x;

  return sum/double(nv);
}

inline double calc(const Values &xv, const Values &yv) {
  assert(xv.size() == yv.size());

  // calc means
  double xmean = mean(xv);
  double ymean = mean(yv);

  // calc correlation
  auto nxy = xv.size();
  if (nxy == 0) return 0.0;

  double sum_dxy = 0.0;
  double sum_dxx = 0.0;
  double sum_dyy = 0.0;

  for (size_t i = 0; i < nxy; ++i) {
    double dx = xv[i] - xmean;
    double dy = yv[i] - ymean;

    double dxy = dx*dy;

    sum_dxy += dxy;

    sum_dxx += dx*dx;
    sum_dyy += dy*dy;
  }

  double d = sqrt(sum_dxx)*sqrt(sum_dyy);

  double p = (d > 0 ? sum_dxy/d : 0.0);

  return p;
}

inline double stddev(const Values &v) {
  // calc means
  double m = mean(v);

  // calc correlation
  auto n = v.size();
  if (n == 0) return 0.0;

  double sum2 = 0.0;

  for (size_t i = 0; i < n; ++i) {
    double d = v[i] - m;

    sum2 += d*d;
  }

  return sqrt(sum2/double(n));
}

}

//------

class CMathBivariate {
 public:
  using Values = std::vector<double>;

 public:
  CMathBivariate(const Values &xv, const Values &yv) :
   xv_(xv), yv_(yv) {
    init();
  }

  double calc(double x, double y) {
    double p1 = 1 - corr_*corr_;
    if (p1 <= 0.0) return 0.0;

    double xxstddev = xstddev_*xstddev_;
    double yystddev = ystddev_*ystddev_;
    double xystddev = xstddev_*ystddev_;

    double f = 2.0*M_PI*xystddev*sqrt(p1);
    if (f <= 0.0) return 0.0;

    double dx = (x - xmean_);
    double dy = (y - ymean_);

    double a = -1.0/(2.0*p1);
    double b = dx*dx/xxstddev + dy*dy/yystddev - 2*corr_*dx*dy/xystddev;

    return (1.0/f)*exp(a*b);
  }

  double xStdDev() const { return xstddev_; }
  double yStdDev() const { return ystddev_; }

 private:
  void init() {
    corr_ = CMathCorrelation::calc(xv_, yv_);

    calcMeanStdDev(xv_, xmean_, xstddev_);
    calcMeanStdDev(yv_, ymean_, ystddev_);
  }

  void calcMeanStdDev(const Values &x, double &mean, double &stddev) const {
    mean   = 0.0;
    stddev = 0.0;

    auto nx = x.size();

    for (size_t i = 0; i < nx; i++) {
      double x1 = x[i];

      mean   += x1;
      stddev += x1*x1;
    }

    if (nx > 0) {
      mean /= double(nx);

      stddev = sqrt(stddev/double(nx) - mean*mean);
    }
  }

 private:
  Values xv_;
  Values yv_;
  double corr_    { 0.0 };
  double xmean_   { 0.0 };
  double xstddev_ { 0.0 };
  double ymean_   { 0.0 };
  double ystddev_ { 0.0 };
};

#endif
