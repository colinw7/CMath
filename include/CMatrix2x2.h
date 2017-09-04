#ifndef CMATRIX_2X2_H
#define CMATRIX_2X2_H

#include <CPoint2D.h>
#include <CVector2D.h>
#include <iostream>
#include <cstring>

/* / m00 m01 \ */
/* \ m10 m11 / */

class CMatrix2x2 {
 private:
  enum Type {
    CMATRIX_2x2_IDENTITY
  };

 private:
  typedef CPoint2D  Point;
  typedef CVector2D Vector;

 public:
  CMatrix2x2() { }

  explicit CMatrix2x2(Type type) {
    if (type == CMATRIX_2x2_IDENTITY)
      setIdentity();
    else
      assert(false && "Bad Matrix Type");
  }

  CMatrix2x2(double a, double b, double c, double d) :
    m00_(a), m01_(b), m10_(c), m11_(d) {
  }

  CMatrix2x2(const double *m, int n) {
    if (n == 4)
      setValues(m[0], m[1], m[2], m[3]);
    else
     assert(false && "Invalid size");
  }

  CMatrix2x2(const CMatrix2x2 &a) {
    memcpy(&m00_, &a.m00_, 4*sizeof(double));
  }

 ~CMatrix2x2() { }

  void setIdentity() {
    m00_ = 1.0, m01_ = 0.0;
    m10_ = 0.0, m11_ = 1.0;
  }

  void setScale(double sx, double sy) {
    m00_ = sx , m01_ = 0.0;
    m10_ = 0.0, m11_ = sy ;
  }

  void setRotation(double a) {
    double c = ::cos(a);
    double s = ::sin(s);

    m00_ =  c , m01_ =  s;
    m10_ = -s , m11_ =  c;
  }

  void setValues(double a, double b, double c, double d) {
    m00_ = a, m01_ = b;
    m10_ = c, m11_ = d;
  }

  void getValues(double *a, double *b, double *c, double *d) const {
    if (a ) *a = m00_;
    if (b ) *b = m01_;
    if (c ) *c = m10_;
    if (d ) *d = m11_;
  }

  void setColumn(int c, double x, double y) {
    assert(c < 2);

    double *m = &(&m00_)[c];

    m[0] = x, m[2] = y;
  }

  void setColumn(int c, const Point &point) {
    setColumn(c, point.x, point.y);
  }

  void setColumn(int c, const Vector &vector) {
    setColumn(c, vector.getX(), vector.getY());
  }

  void getColumn(int c, double *x, double *y) {
    assert(c < 2);

    double *m = &(&m00_)[c];

    if (x) *x = m[0];
    if (y) *y = m[2];
  }

  void setRow(int r, double x, double y) {
    assert(r < 2);

    double *m = &(&m00_)[r*2];

    m[0] = x, m[1] = y;
  }

  void setRow(int r, const Point &point) {
    setRow(r, point.x, point.y);
  }

  void setRow(int r, const Vector &vector) {
    setRow(r, vector.getX(), vector.getY());
  }

  void getRow(int r, double *x, double *y) {
    assert(r < 2);

    double *m = &(&m00_)[r*2];

    if (x) *x = m[0];
    if (y) *y = m[1];
  }

  void multiplyPoint(double xi, double yi, double *xo, double *yo) const {
    *xo = m00_*xi + m01_*yi;
    *yo = m10_*xi + m11_*yi;
  }

  void multiplyPoint(const CPoint2D &point1, CPoint2D &point2) const {
    point2.x = m00_*point1.x + m01_*point1.y;
    point2.y = m10_*point1.x + m11_*point1.y;
  }

  void multiplyVector(const Vector &ivector, Vector &ovector) const {
    double ix, iy;

    ivector.getXY(&ix, &iy);

    double ox = m00_*ix + m01_*iy;
    double oy = m10_*ix + m11_*iy;

    ovector.setXY(ox, oy);
  }

  void preMultiplyPoint(double xi, double yi, double *xo, double *yo) const {
    *xo = m00_*xi + m10_*yi;
    *yo = m01_*xi + m11_*yi;
  }

  void preMultiplyPoint(const Point &ipoint, Point &opoint) const {
    opoint.x = m00_*ipoint.x + m10_*ipoint.y;
    opoint.y = m01_*ipoint.x + m11_*ipoint.y;
  }

  void preMultiplyVector(const Vector &ivector, Vector &ovector) const {
    double ix, iy;

    ivector.getXY(&ix, &iy);

    double ox = m00_*ix + m10_*iy;
    double oy = m01_*ix + m11_*iy;

    ovector.setXY(ox, oy);
  }

  void transpose() {
    std::swap(m10_, m01_);
  }

  CMatrix2x2 transposed() const {
    return CMatrix2x2(m00_, m10_, m01_, m11_);
  }

  bool invert(CMatrix2x2 &imatrix) const {
    double d = determinant();

    if (::fabs(d) < 1E-6)
      return false;

    double id = 1.0/d;

    imatrix.m00_ =  m11_*id;
    imatrix.m01_ = -m01_*id;
    imatrix.m10_ = -m10_*id;
    imatrix.m11_ =  m00_*id;

    return true;
  }

  CMatrix2x2 inverse() const {
    CMatrix2x2 imatrix;

    if (! invert(imatrix))
      assert(false && "Divide by zero");

    return imatrix;
  }

  double determinant() const {
    return (m00_*m11_ - m01_*m10_);
  }

  void normalize() {
    double d = determinant();

    double id = 1.0/d;

    for (int i = 0; i < 4; ++i)
      (&m00_)[i] *= id;
  }

  static CMatrix2x2 *newIdentityMatrix() {
    CMatrix2x2 *m = new CMatrix2x2();

    m->setIdentity();

    return m;
  }

  static bool solveAXeqB(const CMatrix2x2 &a, CPoint2D &x, const CPoint2D &b) {
    double det_a = a.determinant();

    if (::fabs(det_a) < 1E-6)
      return false;

    double idet_a = 1.0/det_a;

    CMatrix2x2 t(a);

    t.setColumn(0, b.x, b.y);

    double det_t = t.determinant();

    x.x = det_t*idet_a;

    t = a;

    t.setColumn(1, b.x, b.y);

    det_t = t.determinant();

    x.y = det_t*idet_a;

    return true;
  }

  void zero() { memset(&m00_, 0, 4*sizeof(double)); }

  CMatrix2x2 &operator=(const CMatrix2x2 &a) {
    memcpy(&m00_, &a.m00_, 4*sizeof(double));

    return *this;
  }

  CMatrix2x2 &operator+=(const CMatrix2x2 &b) {
    m00_ += b.m00_; m01_ += b.m01_,
    m10_ += b.m10_; m11_ += b.m11_;

    return *this;
  }

  CMatrix2x2 operator+(const CMatrix2x2 &b) {
    CMatrix2x2 c = *this;

    c += b;

    return c;
  }

  CMatrix2x2 &operator-=(const CMatrix2x2 &b) {
    m00_ -= b.m00_; m01_ -= b.m01_,
    m10_ -= b.m10_; m11_ -= b.m11_;

    return *this;
  }

  CMatrix2x2 operator-(const CMatrix2x2 &b) {
    CMatrix2x2 c = *this;

    c -= b;

    return c;
  }

  CMatrix2x2 &operator*=(const CMatrix2x2 &b) {
    CMatrix2x2 a;

    memcpy(&a.m00_, &m00_, 4*sizeof(double));

    m00_ = a.m00_*b.m00_ + a.m01_*b.m10_;
    m01_ = a.m00_*b.m01_ + a.m01_*b.m11_;

    m10_ = a.m10_*b.m00_ + a.m11_*b.m10_;
    m11_ = a.m10_*b.m01_ + a.m11_*b.m11_;

    return *this;
  }

  CMatrix2x2 operator*(const CMatrix2x2 &b) {
    CMatrix2x2 c = *this;

    c *= b;

    return c;
  }

  CMatrix2x2 &operator*=(double s) {
    CMatrix2x2 a = *this;

    m00_ = a.m00_*s; m01_ = a.m01_*s;
    m10_ = a.m10_*s; m11_ = a.m11_*s;

    return *this;
  }

  CMatrix2x2 operator*(double s) {
    CMatrix2x2 c = *this;

    c *= s;

    return c;
  }

  friend Point operator*(const CMatrix2x2 &m, const Point &p) {
    Point p1;

    m.multiplyPoint(p, p1);

    return p1;
  }

  friend Point operator*(const Point &p, const CMatrix2x2 &m) {
    Point p1;

    m.preMultiplyPoint(p, p1);

    return p1;
  }

  friend Vector operator*(const CMatrix2x2 &m, const Vector &v) {
    Vector v1;

    m.multiplyVector(v, v1);

    return v1;
  }

  friend Vector operator*(const Vector &v, const CMatrix2x2 &m) {
    Vector v1;

    m.preMultiplyVector(v, v1);

    return v1;
  }

  CMatrix2x2 &operator/=(const CMatrix2x2 &b) {
    CMatrix2x2 bi;

    if (! b.invert(bi)) {
      assert(false && "Divide by zero");
      return *this;
    }

    return (*this) *= bi;
  }

  CMatrix2x2 operator/(const CMatrix2x2 &b) {
    CMatrix2x2 c = *this;

    c /= b;

    return c;
  }

  void setValue(unsigned int i, double value) {
    (&m00_)[i] = value;
  }

  void setValue(unsigned int i, unsigned int j, double value) {
    assert(i < 2 && j < 2);

    double &m = (&m00_)[2*j + i];

    m = value;
  }

  double getValue(unsigned int i) { return (&m00_)[i]; }

  double getValue(unsigned int i, unsigned int j) const {
    assert(i < 2 && j < 2);

    const double &m = (&m00_)[2*j + i];

    return m;
  }

  double operator[](unsigned int i) { return (&m00_)[i]; }

  const double &operator[](unsigned int i) const { return (&m00_)[i]; }

  void print(std::ostream &os) const {
    os << "(" << m00_ << "," << m01_ << ")" << std::endl;
    os << "(" << m10_ << "," << m11_ << ")" << std::endl;
  }

  friend std::ostream &operator<<(std::ostream &os, const CMatrix2x2 &matrix) {
    matrix.print(os);

    return os;
  }

 private:
  static double calcDeterminant(double m00, double m01, double m10, double m11) {
    return m00*m11 - m01*m10;
  }

 private:
  double m00_ { 0 }, m01_ { 0 };
  double m10_ { 0 }, m11_ { 0 };
};

#endif
