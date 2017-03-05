#ifndef CTRANSFORM_2D_H
#define CTRANSFORM_2D_H

#include <CMatrix2D.h>

class CTransform2D {
 public:
  CTransform2D() {
    calcMatrix();
  }

  CTransform2D(double xmin1, double ymin1, double xmax1, double ymax1,
               double xmin2, double ymin2, double xmax2, double ymax2) :
   xmin1_(xmin1), ymin1_(ymin1), xmax1_(xmax1), ymax1_(ymax1),
   xmin2_(xmin2), ymin2_(ymin2), xmax2_(xmax2), ymax2_(ymax2) {
    calcMatrix();
  }

  CTransform2D(const CTransform2D &transform) :
   xmin1_(transform.xmin1_), ymin1_(transform.ymin1_),
   xmax1_(transform.xmax1_), ymax1_(transform.ymax1_),
   xmin2_(transform.xmin2_), ymin2_(transform.ymin2_),
   xmax2_(transform.xmax2_), ymax2_(transform.ymax2_) {
    calcMatrix();
  }

 ~CTransform2D() { }

  CTransform2D &operator=(const CTransform2D &transform) {
    xmin1_ = transform.xmin1_; ymin1_ = transform.ymin1_;
    xmax1_ = transform.xmax1_; ymax1_ = transform.ymax1_;

    xmin2_ = transform.xmin2_; ymin2_ = transform.ymin2_;
    xmax2_ = transform.xmax2_; ymax2_ = transform.ymax2_;

    inverse_set_ = false;

    calcMatrix();

    return *this;
  }

  void setValues(double xmin1, double ymin1, double xmax1, double ymax1,
                 double xmin2, double ymin2, double xmax2, double ymax2) {
    xmin1_ = xmin1; ymin1_ = ymin1; xmax1_ = xmax1; ymax1_ = ymax1;
    xmin2_ = xmin2; ymin2_ = ymin2; xmax2_ = xmax2; ymax2_ = ymax2;

    inverse_set_ = false;

    calcMatrix();
  }

  void setFrom(double xmin1, double ymin1, double xmax1, double ymax1) {
    xmin1_ = xmin1; ymin1_ = ymin1; xmax1_ = xmax1; ymax1_ = ymax1;

    inverse_set_ = false;

    calcMatrix();
  }

  void setTo(double xmin2, double ymin2, double xmax2, double ymax2) {
    xmin2_ = xmin2; ymin2_ = ymin2; xmax2_ = xmax2; ymax2_ = ymax2;

    inverse_set_ = false;

    calcMatrix();
  }

  void conv(double x1, double y1, double *x2, double *y2) const {
    nmatrix_.multiplyPoint(x1, y1, x2, y2);
  }

  void conv(const CPoint2D &point1, CPoint2D &point2) const {
    nmatrix_.multiplyPoint(point1, point2);
  }

  void iconv(double x1, double y1, double *x2, double *y2) const {
    getIMatrix().multiplyPoint(x1, y1, x2, y2);
  }

  void iconv(const CPoint2D &point1, CPoint2D &point2) const {
    getIMatrix().multiplyPoint(point1, point2);
  }

  double getXMin1() const { return xmin1_; }
  double getYMin1() const { return ymin1_; }

  double getXMax1() const { return xmax1_; }
  double getYMax1() const { return ymax1_; }

  double getXMin2() const { return xmin2_; }
  double getYMin2() const { return ymin2_; }

  double getXMax2() const { return xmax2_; }
  double getYMax2() const { return ymax2_; }

  const CMatrix2D &getMatrix () const { return nmatrix_; }

  const CMatrix2D &getIMatrix() const {
    if (! inverse_set_) {
      CTransform2D *th = const_cast<CTransform2D *>(this);

      nmatrix_.invert(th->imatrix_);

      inverse_set_ = true;
    }

    return imatrix_;
  }

 private:
  void calcMatrix() const {
    nmatrix_.setTransform(xmin1_, ymin1_, xmax1_, ymax1_, xmin2_, ymin2_, xmax2_, ymax2_);

    inverse_set_ = false;
  }

 private:
  double xmin1_ { 0 }, ymin1_ { 0 };
  double xmax1_ { 1 }, ymax1_ { 1 };
  double xmin2_ { 0 }, ymin2_ { 0 };
  double xmax2_ { 1 }, ymax2_ { 1 };

  mutable CMatrix2D nmatrix_;
  mutable CMatrix2D imatrix_;
  mutable bool      inverse_set_ { false };
};

#endif
