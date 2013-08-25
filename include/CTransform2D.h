#ifndef CTRANSFORM_2D_H
#define CTRANSFORM_2D_H

#include <CMatrix2D.h>

template<typename T>
class CTransform2DT {
 private:
  T xmin1_, ymin1_;
  T xmax1_, ymax1_;
  T xmin2_, ymin2_;
  T xmax2_, ymax2_;

  mutable CMatrix2DT<T> nmatrix_;
  mutable CMatrix2DT<T> imatrix_;
  mutable bool          inverse_set_;

 public:
  CTransform2DT() :
   xmin1_(0.0), ymin1_(0.0), xmax1_(1.0), ymax1_(1.0),
   xmin2_(0.0), ymin2_(0.0), xmax2_(1.0), ymax2_(1.0),
   nmatrix_(), imatrix_(), inverse_set_(false)
  {
    calcMatrix();
  }

  CTransform2DT(T xmin1, T ymin1, T xmax1, T ymax1, T xmin2, T ymin2, T xmax2, T ymax2) :
   xmin1_(xmin1), ymin1_(ymin1), xmax1_(xmax1), ymax1_(ymax1),
   xmin2_(xmin2), ymin2_(ymin2), xmax2_(xmax2), ymax2_(ymax2),
   nmatrix_(), imatrix_(), inverse_set_(false)
  {
    calcMatrix();
  }

  CTransform2DT(const CTransform2DT &transform) :
   xmin1_(transform.xmin1_), ymin1_(transform.ymin1_),
   xmax1_(transform.xmax1_), ymax1_(transform.ymax1_),
   xmin2_(transform.xmin2_), ymin2_(transform.ymin2_),
   xmax2_(transform.xmax2_), ymax2_(transform.ymax2_),
   nmatrix_(), imatrix_(), inverse_set_(false)
  {
    calcMatrix();
  }

 ~CTransform2DT() { }

  CTransform2DT &operator=(const CTransform2DT &transform) {
    xmin1_ = transform.xmin1_; ymin1_ = transform.ymin1_;
    xmax1_ = transform.xmax1_; ymax1_ = transform.ymax1_;

    xmin2_ = transform.xmin2_; ymin2_ = transform.ymin2_;
    xmax2_ = transform.xmax2_; ymax2_ = transform.ymax2_;

    inverse_set_ = false;

    calcMatrix();

    return *this;
  }

  void setValues(T xmin1, T ymin1, T xmax1, T ymax1, T xmin2, T ymin2, T xmax2, T ymax2) {
    xmin1_ = xmin1; ymin1_ = ymin1; xmax1_ = xmax1; ymax1_ = ymax1;
    xmin2_ = xmin2; ymin2_ = ymin2; xmax2_ = xmax2; ymax2_ = ymax2;

    inverse_set_ = false;

    calcMatrix();
  }

  void setFrom(T xmin1, T ymin1, T xmax1, T ymax1) {
    xmin1_ = xmin1; ymin1_ = ymin1; xmax1_ = xmax1; ymax1_ = ymax1;

    inverse_set_ = false;

    calcMatrix();
  }

  void setTo(T xmin2, T ymin2, T xmax2, T ymax2) {
    xmin2_ = xmin2; ymin2_ = ymin2; xmax2_ = xmax2; ymax2_ = ymax2;

    inverse_set_ = false;

    calcMatrix();
  }

  void conv(T x1, T y1, T *x2, T *y2) const {
    nmatrix_.multiplyPoint(x1, y1, x2, y2);
  }

  void conv(const CPoint2DT<T> &point1, CPoint2DT<T> &point2) const {
    nmatrix_.multiplyPoint(point1, point2);
  }

  void iconv(T x1, T y1, T *x2, T *y2) const {
    getIMatrix().multiplyPoint(x1, y1, x2, y2);
  }

  void iconv(const CPoint2DT<T> &point1, CPoint2DT<T> &point2) const {
    getIMatrix().multiplyPoint(point1, point2);
  }

  T getXMin1() const { return xmin1_; }
  T getYMin1() const { return ymin1_; }

  T getXMax1() const { return xmax1_; }
  T getYMax1() const { return ymax1_; }

  T getXMin2() const { return xmin2_; }
  T getYMin2() const { return ymin2_; }

  T getXMax2() const { return xmax2_; }
  T getYMax2() const { return ymax2_; }

  const CMatrix2DT<T> &getMatrix () const { return nmatrix_; }

  const CMatrix2DT<T> &getIMatrix() const {
    if (! inverse_set_) {
      CTransform2DT *th = const_cast<CTransform2DT *>(this);

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
};

typedef CTransform2DT<double> CTransform2D;

#endif
