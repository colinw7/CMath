#ifndef CMATRIX_4X4_H
#define CMATRIX_4X4_H

#include <CPoint4D.h>
#include <cstring>
#include <iostream>

/* / m00 m01 m02 m03 \ */
/* | m10 m11 m12 m13 | */
/* | m20 m21 m22 m23 | */
/* \ m30 m31 m32 m33 / */

class CMatrix4x4 {
 public:
  CMatrix4x4() { }

  CMatrix4x4(double m00, double m01, double m02, double m03,
             double m10, double m11, double m12, double m13,
             double m20, double m21, double m22, double m23,
             double m30, double m31, double m32, double m33) :
   m00_(m00), m01_(m01), m02_(m02), m03_(m03),
   m10_(m10), m11_(m11), m12_(m12), m13_(m13),
   m20_(m20), m21_(m21), m22_(m22), m23_(m23),
   m30_(m30), m31_(m31), m32_(m32), m33_(m33) {
  }

  CMatrix4x4(const CMatrix4x4 &a) {
    memcpy(&m00_, &a.m00_, 16*sizeof(double));
  }

  ~CMatrix4x4() { }

  void setIdentity() {
    m00_ = 1.0, m01_ = 0.0, m02_ = 0.0; m03_ = 0.0;
    m10_ = 0.0, m11_ = 1.0, m12_ = 0.0; m13_ = 0.0;
    m20_ = 0.0, m21_ = 0.0, m22_ = 1.0; m23_ = 0.0;
    m30_ = 0.0, m31_ = 0.0, m32_ = 0.0; m33_ = 1.0;
  }

  void setTranslation(double tx, double ty, double tz) {
    m00_ = 1.0, m01_ = 0.0, m02_ = 0.0; m03_ = tx ;
    m10_ = 0.0, m11_ = 1.0, m12_ = 0.0; m13_ = ty ;
    m20_ = 0.0, m21_ = 0.0, m22_ = 1.0; m23_ = tz ;
    m30_ = 0.0, m31_ = 0.0, m32_ = 0.0; m33_ = 1.0;
  }

  void setScale(double sx, double sy, double sz) {
    m00_ = sx , m01_ = 0.0, m02_ = 0.0; m03_ = 0.0;
    m10_ = 0.0, m11_ = sy , m12_ = 0.0; m03_ = 0.0;
    m20_ = 0.0, m21_ = 0.0, m22_ = sz ; m03_ = 0.0;
    m30_ = 0.0, m31_ = 0.0, m32_ = 0.0; m03_ = 0.0;
  }

  void setValues(double m00, double m01, double m02, double m03,
                 double m10, double m11, double m12, double m13,
                 double m20, double m21, double m22, double m23,
                 double m30, double m31, double m32, double m33) {
    m00_ = m00; m01_ = m01; m02_ = m02; m03_ = m03;
    m10_ = m10; m11_ = m11; m12_ = m12; m13_ = m13;
    m20_ = m20; m21_ = m21; m22_ = m22; m23_ = m23;
    m30_ = m30; m31_ = m31; m32_ = m32; m33_ = m33;
  }

  void getValues(double *m00, double *m01, double *m02, double *m03,
                 double *m10, double *m11, double *m12, double *m13,
                 double *m20, double *m21, double *m22, double *m23,
                 double *m30, double *m31, double *m32, double *m33) {
    if (m00) *m00 = m00_; if (m01) *m01 = m01_;
    if (m01) *m02 = m02_; if (m03) *m03 = m03_;
    if (m10) *m10 = m10_; if (m11) *m11 = m11_;
    if (m01) *m12 = m12_; if (m13) *m13 = m13_;
    if (m20) *m20 = m20_; if (m21) *m21 = m21_;
    if (m01) *m22 = m22_; if (m23) *m23 = m23_;
    if (m30) *m30 = m30_; if (m31) *m31 = m31_;
    if (m01) *m32 = m32_; if (m33) *m33 = m33_;
  }

  void setValue(unsigned int i, double value) {
    (&m00_)[i] = value;
  }

  void setValue(unsigned int i, unsigned int j, double value) {
    assert(i < 4 && j < 4);

    double &m = (&m00_)[4*j + i];

    m = value;
  }

  double getValue(unsigned int i) const {
    return (&m00_)[i];
  }

  double getValue(unsigned int i, unsigned int j) const {
    assert(i < 4 && j < 4);

    const double &m = (&m00_)[4*j + i];

    return m;
  }

  void setRow(int r, double x, double y, double z, double w) {
    assert(r < 4);

    double *m = &(&m00_)[r*4];

    m[0] = x, m[1] = y; m[2] = z; m[3] = w;
  }

  void setColumn(int c, double x, double y, double z, double w) {
    assert(c < 4);

    double *m = &(&m00_)[c];

    m[0] = x, m[4] = y; m[8] = z; m[12] = w;
  }

  void getRow(int r, double *x, double *y, double *z, double *w) {
    assert(r < 4);

    double *m = &(&m00_)[r*4];

    if (x) *x = m[0];
    if (y) *y = m[1];
    if (z) *z = m[2];
    if (w) *w = m[3];
  }

  void getColumn(int c, double *x, double *y, double *z, double *w) {
    assert(c < 4);

    double *m = &(&m00_)[c];

    if (x) *x = m[0];
    if (y) *y = m[4];
    if (z) *z = m[8];
    if (w) *w = m[12];
  }

  bool invert(CMatrix4x4 &imatrix) const {
    double det = determinant();

    if (::fabs(det) == 0.0)
      return false;

    double idet = 1.0/det;

    imatrix.m00_ =  idet*det3x3(m11_, m12_, m13_,
                                m21_, m22_, m23_,
                                m31_, m32_, m33_);
    imatrix.m10_ = -idet*det3x3(m10_, m12_, m13_,
                                m20_, m22_, m23_,
                                m30_, m32_, m33_);
    imatrix.m20_ =  idet*det3x3(m10_, m11_, m13_,
                                m20_, m21_, m23_,
                                m30_, m31_, m33_);
    imatrix.m30_ = -idet*det3x3(m10_, m11_, m12_,
                                m20_, m21_, m22_,
                                m30_, m31_, m32_);

    imatrix.m01_ = -idet*det3x3(m01_, m02_, m03_,
                                m21_, m22_, m23_,
                                m31_, m32_, m33_);
    imatrix.m11_ =  idet*det3x3(m00_, m02_, m03_,
                                m20_, m22_, m23_,
                                m30_, m32_, m33_);
    imatrix.m21_ = -idet*det3x3(m00_, m01_, m03_,
                                m20_, m21_, m23_,
                                m30_, m31_, m33_);
    imatrix.m31_ =  idet*det3x3(m00_, m01_, m02_,
                                m20_, m21_, m22_,
                                m30_, m31_, m32_);

    imatrix.m02_ =  idet*det3x3(m01_, m02_, m03_,
                                m11_, m12_, m13_,
                                m31_, m32_, m33_);
    imatrix.m12_ = -idet*det3x3(m00_, m02_, m03_,
                                m10_, m12_, m13_,
                                m30_, m32_, m33_);
    imatrix.m22_ =  idet*det3x3(m00_, m01_, m03_,
                                m10_, m11_, m13_,
                                m30_, m31_, m33_);
    imatrix.m32_ = -idet*det3x3(m00_, m01_, m02_,
                                m10_, m11_, m12_,
                                m30_, m31_, m32_);

    imatrix.m03_ = -idet*det3x3(m01_, m02_, m03_,
                                m11_, m12_, m13_,
                                m21_, m22_, m23_);
    imatrix.m13_ =  idet*det3x3(m00_, m02_, m03_,
                                m10_, m12_, m13_,
                                m20_, m22_, m23_);
    imatrix.m23_ = -idet*det3x3(m00_, m01_, m03_,
                                m10_, m11_, m13_,
                                m20_, m21_, m23_);
    imatrix.m33_ =  idet*det3x3(m00_, m01_, m02_,
                                m10_, m11_, m12_,
                                m20_, m21_, m22_);

    return true;
  }

  double determinant() const {
    return
      (m00_*det3x3(m11_, m12_, m13_, m21_, m22_, m23_, m31_, m32_, m33_) -
       m01_*det3x3(m10_, m12_, m13_, m20_, m22_, m23_, m30_, m32_, m33_) +
       m02_*det3x3(m10_, m11_, m13_, m20_, m21_, m23_, m30_, m31_, m33_) -
       m03_*det3x3(m10_, m11_, m12_, m20_, m21_, m22_, m30_, m31_, m32_));
  }

  void transpose() {
    std::swap(m10_, m01_);
    std::swap(m20_, m02_);
    std::swap(m21_, m12_);
    std::swap(m30_, m03_);
    std::swap(m31_, m13_);
    std::swap(m32_, m23_);
  }

  void multiplyPoint(double  xi, double  yi, double  zi, double wi,
                     double *xo, double *yo, double *zo, double *wo) const {
    *xo = m00_*xi + m01_*yi + m02_*zi + m03_*wi;
    *yo = m10_*xi + m11_*yi + m12_*zi + m13_*wi;
    *zo = m20_*xi + m21_*yi + m22_*zi + m23_*wi;
    *wo = m30_*xi + m31_*yi + m32_*zi + m33_*wi;
  }

  void multiplyPoint(const CPoint4D &p1, CPoint4D &p2) const {
    p2.x = m00_*p1.x + m01_*p1.y + m02_*p1.z + m03_*p1.w;
    p2.y = m10_*p1.x + m11_*p1.y + m12_*p1.z + m13_*p1.w;
    p2.z = m20_*p1.x + m21_*p1.y + m22_*p1.z + m23_*p1.w;
    p2.w = m30_*p1.x + m31_*p1.y + m32_*p1.z + m33_*p1.w;
  }

  void zero() { memset(&m00_, 0, 16*sizeof(double)); }

  void print(std::ostream &os) const {
    os << "(" << m00_ << "," << m01_ << "," <<
                 m02_ << "," << m03_ << ")" << std::endl;
    os << "(" << m10_ << "," << m11_ << "," <<
                 m12_ << "," << m13_ << ")" << std::endl;
    os << "(" << m20_ << "," << m21_ << "," <<
                 m22_ << "," << m23_ << ")" << std::endl;
    os << "(" << m30_ << "," << m31_ << "," <<
                 m32_ << "," << m33_ << ")" << std::endl;
  }

  CMatrix4x4 &operator=(const CMatrix4x4 &a) {
    memcpy(&m00_, &a.m00_, 16*sizeof(double));

    return *this;
  }

  CMatrix4x4 &operator+=(const CMatrix4x4 &b) {
    m00_ += b.m00_; m01_ += b.m01_; m02_ += b.m02_; m03_ += b.m03_;
    m10_ += b.m10_; m11_ += b.m11_; m12_ += b.m12_; m13_ += b.m13_;
    m20_ += b.m20_; m21_ += b.m21_; m22_ += b.m22_; m23_ += b.m23_;
    m30_ += b.m30_; m31_ += b.m31_; m32_ += b.m32_; m33_ += b.m33_;

    return *this;
  }

  CMatrix4x4 operator+(const CMatrix4x4 &b) {
    return CMatrix4x4(
      m00_ + b.m00_, m01_ + b.m01_, m02_ + b.m02_, m03_ + b.m03_,
      m10_ + b.m10_, m11_ + b.m11_, m12_ + b.m12_, m13_ + b.m13_,
      m20_ + b.m20_, m21_ + b.m21_, m22_ + b.m22_, m23_ + b.m23_,
      m30_ + b.m30_, m31_ + b.m31_, m32_ + b.m32_, m33_ + b.m33_);
  }

  CMatrix4x4 &operator-=(const CMatrix4x4 &b) {
    m00_ -= b.m00_; m01_ -= b.m01_; m02_ -= b.m02_; m03_ -= b.m03_;
    m10_ -= b.m10_; m11_ -= b.m11_; m12_ -= b.m12_; m13_ -= b.m13_;
    m20_ -= b.m20_; m21_ -= b.m21_; m22_ -= b.m22_; m23_ -= b.m23_;
    m30_ -= b.m30_; m31_ -= b.m31_; m32_ -= b.m32_; m33_ -= b.m33_;

    return *this;
  }

  CMatrix4x4 operator-(const CMatrix4x4 &b) {
    return CMatrix4x4(
      m00_ - b.m00_, m01_ - b.m01_, m02_ - b.m02_, m03_ - b.m03_,
      m10_ - b.m10_, m11_ - b.m11_, m12_ - b.m12_, m13_ - b.m13_,
      m20_ - b.m20_, m21_ - b.m21_, m22_ - b.m22_, m23_ - b.m23_,
      m30_ - b.m30_, m31_ - b.m31_, m32_ - b.m32_, m33_ - b.m33_);
  }

  CMatrix4x4 &operator*=(const CMatrix4x4 &b) {
    CMatrix4x4 a;

    memcpy(&a.m00_, &m00_, 16*sizeof(double));

    m00_ = a.m00_*b.m00_ + a.m01_*b.m10_ + a.m02_*b.m20_ + a.m03_*b.m30_;
    m01_ = a.m00_*b.m01_ + a.m01_*b.m11_ + a.m02_*b.m21_ + a.m03_*b.m31_;
    m02_ = a.m00_*b.m02_ + a.m01_*b.m12_ + a.m02_*b.m22_ + a.m03_*b.m32_;
    m03_ = a.m00_*b.m03_ + a.m01_*b.m13_ + a.m02_*b.m23_ + a.m03_*b.m33_;

    m10_ = a.m10_*b.m00_ + a.m11_*b.m10_ + a.m12_*b.m20_ + a.m13_*b.m30_;
    m11_ = a.m10_*b.m01_ + a.m11_*b.m11_ + a.m12_*b.m21_ + a.m13_*b.m31_;
    m12_ = a.m10_*b.m02_ + a.m11_*b.m12_ + a.m12_*b.m22_ + a.m13_*b.m32_;
    m13_ = a.m10_*b.m03_ + a.m11_*b.m13_ + a.m12_*b.m23_ + a.m13_*b.m33_;

    m20_ = a.m20_*b.m00_ + a.m21_*b.m10_ + a.m22_*b.m20_ + a.m23_*b.m30_;
    m21_ = a.m20_*b.m01_ + a.m21_*b.m11_ + a.m22_*b.m21_ + a.m23_*b.m31_;
    m22_ = a.m20_*b.m02_ + a.m21_*b.m12_ + a.m22_*b.m22_ + a.m23_*b.m32_;
    m23_ = a.m20_*b.m03_ + a.m21_*b.m13_ + a.m22_*b.m23_ + a.m23_*b.m33_;

    m30_ = a.m30_*b.m00_ + a.m31_*b.m10_ + a.m32_*b.m20_ + a.m33_*b.m30_;
    m31_ = a.m30_*b.m01_ + a.m31_*b.m11_ + a.m32_*b.m21_ + a.m33_*b.m31_;
    m32_ = a.m30_*b.m02_ + a.m31_*b.m12_ + a.m32_*b.m22_ + a.m33_*b.m32_;
    m33_ = a.m30_*b.m03_ + a.m31_*b.m13_ + a.m32_*b.m23_ + a.m33_*b.m33_;

    return *this;
  }

  CMatrix4x4 operator*(const CMatrix4x4 &b) const {
    return CMatrix4x4(
      m00_*b.m00_ + m01_*b.m10_ + m02_*b.m20_ + m03_*b.m30_,
      m00_*b.m01_ + m01_*b.m11_ + m02_*b.m21_ + m03_*b.m31_,
      m00_*b.m02_ + m01_*b.m12_ + m02_*b.m22_ + m03_*b.m32_,
      m00_*b.m03_ + m01_*b.m13_ + m02_*b.m23_ + m03_*b.m33_,
      m10_*b.m00_ + m11_*b.m10_ + m12_*b.m20_ + m13_*b.m30_,
      m10_*b.m01_ + m11_*b.m11_ + m12_*b.m21_ + m13_*b.m31_,
      m10_*b.m02_ + m11_*b.m12_ + m12_*b.m22_ + m13_*b.m32_,
      m10_*b.m03_ + m11_*b.m13_ + m12_*b.m23_ + m13_*b.m33_,
      m20_*b.m00_ + m21_*b.m10_ + m22_*b.m20_ + m23_*b.m30_,
      m20_*b.m01_ + m21_*b.m11_ + m22_*b.m21_ + m23_*b.m31_,
      m20_*b.m02_ + m21_*b.m12_ + m22_*b.m22_ + m23_*b.m32_,
      m20_*b.m03_ + m21_*b.m13_ + m22_*b.m23_ + m23_*b.m33_,
      m30_*b.m00_ + m31_*b.m10_ + m32_*b.m20_ + m33_*b.m30_,
      m30_*b.m01_ + m31_*b.m11_ + m32_*b.m21_ + m33_*b.m31_,
      m30_*b.m02_ + m31_*b.m12_ + m32_*b.m22_ + m33_*b.m32_,
      m30_*b.m03_ + m31_*b.m13_ + m32_*b.m23_ + m33_*b.m33_);
  }

  CMatrix4x4 &operator*=(double s) {
    m00_ *= s; m01_ *= s; m02_ *= s; m03_ *= s;
    m10_ *= s; m11_ *= s; m12_ *= s; m13_ *= s;
    m20_ *= s; m21_ *= s; m22_ *= s; m23_ *= s;
    m30_ *= s; m31_ *= s; m32_ *= s; m33_ *= s;

    return *this;
  }

  CMatrix4x4 operator*(double s) {
    return CMatrix4x4(m00_*s, m01_*s, m02_*s, m03_*s,
                      m10_*s, m11_*s, m12_*s, m13_*s,
                      m20_*s, m21_*s, m22_*s, m23_*s,
                      m30_*s, m31_*s, m32_*s, m33_*s);
  }

  CMatrix4x4 &operator/=(const CMatrix4x4 &b) {
    CMatrix4x4 bi;

    if (! b.invert(bi))
      throw "Divide by zero";

    return (*this) *= bi;
  }

  CMatrix4x4 operator/(const CMatrix4x4 &b) {
    CMatrix4x4 bi;

    if (! b.invert(bi))
      throw "Divide by zero";

    return (*this) * bi;
  }

  double        operator[](unsigned int i) { return (&m00_)[i]; }
  const double &operator[](unsigned int i) const { return (&m00_)[i]; }

  friend std::ostream &operator<<(std::ostream &os, const CMatrix4x4 &matrix) {
    matrix.print(os);

    return os;
  }

 private:
  double det3x3(double a, double b, double c, double d, double e, double f,
                double g, double h, double i) const {
    return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
  }

 private:
  double m00_ { 0 }, m01_ { 0 }, m02_ { 0 }, m03_ { 0 };
  double m10_ { 0 }, m11_ { 0 }, m12_ { 0 }, m13_ { 0 };
  double m20_ { 0 }, m21_ { 0 }, m22_ { 0 }, m23_ { 0 };
  double m30_ { 0 }, m31_ { 0 }, m32_ { 0 }, m33_ { 0 };
};

#endif
