#ifndef CTRANSFORM_3D_H
#define CTRANSFORM_3D_H

#include <CMatrix3D.h>

class CTransform3D {
 public:
  CTransform3D(double xmin1, double ymin1, double zmin1,
               double xmax1, double ymax1, double zmax1,
               double xmin2, double ymin2, double zmin2,
               double xmax2, double ymax2, double zmax2) :
   xmin1_(xmin1), ymin1_(ymin1), zmin1_(zmin1),
   xmax1_(xmax1), ymax1_(ymax1), zmax1_(zmax1),
   xmin2_(xmin2), ymin2_(ymin2), zmin2_(zmin2),
   xmax2_(xmax2), ymax2_(ymax2), zmax2_(zmax2),
   nmatrix_(), imatrix_(), inverse_set_(false) {
    calcMatrix();
  }

  CTransform3D(const CTransform3D &t) :
   xmin1_(t.xmin1_), ymin1_(t.ymin1_), zmin1_(t.zmin1_),
   xmax1_(t.xmax1_), ymax1_(t.ymax1_), zmax1_(t.zmax1_),
   xmin2_(t.xmin2_), ymin2_(t.ymin2_), zmin2_(t.zmin2_),
   xmax2_(t.xmax2_), ymax2_(t.ymax2_), zmax2_(t.zmax2_),
   nmatrix_(), imatrix_(), inverse_set_(false) {
    calcMatrix();
  }

 ~CTransform3D() { }

  CTransform3D &operator=(const CTransform3D &transform) {
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

  void conv(double x1, double y1, double z1, double *x2, double *y2, double *z2) const {
    nmatrix_.multiplyPoint(x1, y1, z1, x2, y2, z2);
  }

  void conv(const CPoint3D &point1, CPoint3D &point2) const {
    nmatrix_.multiplyPoint(point1, point2);
  }

  void iconv(double x1, double y1, double z1, double *x2, double *y2, double *z2) const {
    if (! inverse_set_) {
      nmatrix_.invert(imatrix_);

      inverse_set_ = true;
    }

    imatrix_.multiplyPoint(x1, y1, z1, x2, y2, z2);
  }

  void iconv(const CPoint3D &point1, CPoint3D &point2) const {
    if (! inverse_set_) {
      nmatrix_.invert(imatrix_);

      inverse_set_ = true;
    }

    imatrix_.multiplyPoint(point1, point2);
  }

  double getXMin1() const { return xmin1_; }
  double getYMin1() const { return ymin1_; }
  double getZMin1() const { return zmin1_; }

  double getXMax1() const { return xmax1_; }
  double getYMax1() const { return ymax1_; }
  double getZMax1() const { return zmax1_; }

  double getXMin2() const { return xmin2_; }
  double getYMin2() const { return ymin2_; }
  double getZMin2() const { return zmin2_; }

  double getXMax2() const { return xmax2_; }
  double getYMax2() const { return ymax2_; }
  double getZMax2() const { return zmax2_; }

  CMatrix3D *getMatrix () { return &nmatrix_; }
  CMatrix3D *getIMatrix() { return &imatrix_; }

 private:
  void calcMatrix() const {
    nmatrix_.setTransform(xmin1_, ymin1_, zmin1_,
                          xmax1_, ymax1_, zmax1_,
                          xmin2_, ymin2_, zmin2_,
                          xmax2_, ymax2_, zmax2_);

    inverse_set_ = false;
  }

 private:
  double xmin1_, ymin1_, zmin1_;
  double xmax1_, ymax1_, zmax1_;
  double xmin2_, ymin2_, zmin2_;
  double xmax2_, ymax2_, zmax2_;

  mutable CMatrix3D nmatrix_;
  mutable CMatrix3D imatrix_;
  mutable bool      inverse_set_;
};

#endif
