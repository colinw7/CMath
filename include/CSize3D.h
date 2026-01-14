#ifndef CSIZE_3D_H
#define CSIZE_3D_H

#include <CPoint3D.h>

class CSize3D {
 public:
  CSize3D() = default;

  CSize3D(double xs, double ys, double zs) :
   set_(true), xs_(xs), ys_(ys), zs_(zs) {
  }

  bool isSet() const { return set_; }

  void set(double xs, double ys, double zs) {
    set_ = true;
    xs_  = xs;
    ys_  = ys;
    zs_  = zs;
  }

  void get(double *xs, double *ys, double *zs) const {
    assert(set_);

    *xs = xs_;
    *ys = ys_;
    *zs = zs_;
  }

  double getXSize() const { assert(set_); return xs_; }
  double getYSize() const { assert(set_); return ys_; }
  double getZSize() const { assert(set_); return zs_; }

  void setXSize(double xs) { set_ = true; xs_ = xs; }
  void setYSize(double ys) { set_ = true; ys_ = ys; }
  void setZSize(double zs) { set_ = true; zs_ = zs; }

  CPoint3D point() const { return CPoint3D(xs_, ys_, zs_); }

 private:
  bool   set_ { false };
  double xs_  { 0.0 };
  double ys_  { 0.0 };
  double zs_  { 0.0 };
};

#endif
