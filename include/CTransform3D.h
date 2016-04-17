#ifndef CTRANSFORM_3D_H
#define CTRANSFORM_3D_H

#include <CMatrix3D.h>

template<typename T>
class CTransform3DT {
 private:
  T xmin1_, ymin1_, zmin1_;
  T xmax1_, ymax1_, zmax1_;
  T xmin2_, ymin2_, zmin2_;
  T xmax2_, ymax2_, zmax2_;

  mutable CMatrix3DT<T> nmatrix_;
  mutable CMatrix3DT<T> imatrix_;
  mutable bool          inverse_set_;

 public:
  CTransform3DT(T xmin1, T ymin1, T zmin1, T xmax1, T ymax1, T zmax1,
                T xmin2, T ymin2, T zmin2, T xmax2, T ymax2, T zmax2) :
   xmin1_(xmin1), ymin1_(ymin1), zmin1_(zmin1),
   xmax1_(xmax1), ymax1_(ymax1), zmax1_(zmax1),
   xmin2_(xmin2), ymin2_(ymin2), zmin2_(zmin2),
   xmax2_(xmax2), ymax2_(ymax2), zmax2_(zmax2),
   nmatrix_(), imatrix_(), inverse_set_(false) {
    calcMatrix();
  }

  CTransform3DT(const CTransform3DT &t) :
   xmin1_(t.xmin1_), ymin1_(t.ymin1_), zmin1_(t.zmin1_),
   xmax1_(t.xmax1_), ymax1_(t.ymax1_), zmax1_(t.zmax1_),
   xmin2_(t.xmin2_), ymin2_(t.ymin2_), zmin2_(t.zmin2_),
   xmax2_(t.xmax2_), ymax2_(t.ymax2_), zmax2_(t.zmax2_),
   nmatrix_(), imatrix_(), inverse_set_(false) {
    calcMatrix();
  }

 ~CTransform3DT() { }

  CTransform3DT &operator=(const CTransform3DT &transform) {
    xmin1_ = transform.xmin1_;
    ymin1_ = transform.ymin1_;
    zmin1_ = transform.zmin1_;

    xmax1_ = transform.xmax1_;
    ymax1_ = transform.ymax1_;
    zmax1_ = transform.zmax1_;

    xmin2_ = transform.xmin2_;
    ymin2_ = transform.ymin2_;
    zmin2_ = transform.zmin2_;

    xmax2_ = transform.xmax2_;
    ymax2_ = transform.ymax2_;
    zmax2_ = transform.zmax2_;

    inverse_set_ = false;

    calcMatrix();

    return *this;
  }

  void conv(T x1, T y1, T z1, T *x2, T *y2, T *z2) const {
    nmatrix_.multiplyPoint(x1, y1, z1, x2, y2, z2);
  }

  void conv(const CPoint3DT<T> &point1, CPoint3DT<T> &point2) const {
    nmatrix_.multiplyPoint(point1, point2);
  }

  void iconv(T x1, T y1, T z1, T *x2, T *y2, T *z2) const {
    if (! inverse_set_) {
      nmatrix_.invert(imatrix_);

      inverse_set_ = true;
    }

    imatrix_.multiplyPoint(x1, y1, z1, x2, y2, z2);
  }

  void iconv(const CPoint3DT<T> &point1, CPoint3DT<T> &point2) const {
    if (! inverse_set_) {
      nmatrix_.invert(imatrix_);

      inverse_set_ = true;
    }

    imatrix_.multiplyPoint(point1, point2);
  }

  T getXMin1() const { return xmin1_; }
  T getYMin1() const { return ymin1_; }
  T getZMin1() const { return zmin1_; }

  T getXMax1() const { return xmax1_; }
  T getYMax1() const { return ymax1_; }
  T getZMax1() const { return zmax1_; }

  T getXMin2() const { return xmin2_; }
  T getYMin2() const { return ymin2_; }
  T getZMin2() const { return zmin2_; }

  T getXMax2() const { return xmax2_; }
  T getYMax2() const { return ymax2_; }
  T getZMax2() const { return zmax2_; }

  CMatrix3DT<T> *getMatrix () { return &nmatrix_; }
  CMatrix3DT<T> *getIMatrix() { return &imatrix_; }

 private:
  void calcMatrix() const {
    nmatrix_.setTransform(xmin1_, ymin1_, zmin1_,
                          xmax1_, ymax1_, zmax1_,
                          xmin2_, ymin2_, zmin2_,
                          xmax2_, ymax2_, zmax2_);

    inverse_set_ = false;
  }
};

typedef CTransform3DT<double> CTransform3D;

#endif
