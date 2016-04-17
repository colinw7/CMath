#ifndef CMATRIX_2X2_H
#define CMATRIX_2X2_H

#include <CMathGen.h>
#include <CPoint2D.h>
#include <CVector2D.h>
#include <CThrow.h>

/* / a b \ */
/* \ c d / */

template<typename T>
class CMatrix2x2T {
 private:
  enum Type {
    CMATRIX_2x2_IDENTITY
  };

 private:
  typedef CPoint2DT<T>  Point;
  typedef CVector2DT<T> Vector;

 private:
  union {
    T m1_[4];
    T m2_[2][2];

    struct {
      T m00_, m01_;
      T m10_, m11_;
    };

    struct {
      T a_, b_;
      T c_, d_;
    };
  };

 public:
  CMatrix2x2T() { }

  explicit CMatrix2x2T(typename CMatrix2x2T::Type type) {
    if (type == CMATRIX_2x2_IDENTITY)
      setIdentity();
    else
      CTHROW("Bad Matrix Type");
  }

  CMatrix2x2T(T a, T b, T c, T d) :
    m00_(a), m01_(b), m10_(c), m11_(d) {
  }

  CMatrix2x2T(const T *m, int n) {
    if (n == 4)
      setValues(m[0], m[1],
                m[2], m[3]);
    else
     CTHROW("Invalid size");
  }

  CMatrix2x2T(const CMatrix2x2T &a) {
    memcpy(m1_, a.m1_, sizeof(m1_));
  }

 ~CMatrix2x2T() { }

  void setIdentity() {
    m00_ = 1.0, m01_ = 0.0;
    m10_ = 0.0, m11_ = 1.0;
  }

  void setScale(T sx, T sy) {
    m00_ = sx , m01_ = 0.0;
    m10_ = 0.0, m11_ = sy ;
  }

  void setRotation(T a) {
    T c = ::cos(a);
    T s = ::sin(s);

    m00_ =  c , m01_ =  s;
    m10_ = -s , m11_ =  c;
  }

  void setValues(T a, T b, T c, T d) {
    m00_ = a, m01_ = b;
    m10_ = c, m11_ = d;
  }

  void getValues(T *a, T *b, T *c, T *d) const {
    if (a ) *a = m00_;
    if (b ) *b = m01_;
    if (c ) *c = m10_;
    if (d ) *d = m11_;
  }

  void setColumn(int c, T x, T y) {
    m2_[0][c] = x, m2_[1][c] = y;
  }

  void setColumn(int c, const Point &point) {
    m2_[0][c] = point.x, m2_[1][c] = point.y;
  }

  void setColumn(int c, const Vector &vector) {
    vector.getXY(&m2_[0][c], &m2_[1][c]);
  }

  void getColumn(int c, T *x, T *y) {
    if (x) *x = m2_[0][c];
    if (y) *y = m2_[1][c];
  }

  void setRow(int r, T x, T y) {
    m2_[r][0] = x, m2_[r][1] = y;
  }

  void setRow(int r, const Point &point) {
    m2_[r][0] = point.x, m2_[r][1] = point.y;
  }

  void setRow(int r, const Vector &vector) {
    vector.getXY(&m2_[r][0], &m2_[r][1]);
  }

  void getRow(int r, T *x, T *y) {
    if (x) *x = m2_[r][0];
    if (y) *y = m2_[r][1];
  }

  void multiplyPoint(T xi, T yi, T *xo, T *yo) const {
    *xo = m00_*xi + m01_*yi;
    *yo = m10_*xi + m11_*yi;
  }

  void multiplyPoint(const CPoint2DT<T> &point1, CPoint2DT<T> &point2) const {
    point2.x = m00_*point1.x + m01_*point1.y;
    point2.y = m10_*point1.x + m11_*point1.y;
  }

  void multiplyVector(const Vector &ivector, Vector &ovector) const {
    T ix, iy;

    ivector.getXY(&ix, &iy);

    T ox = m00_*ix + m01_*iy;
    T oy = m10_*ix + m11_*iy;

    ovector.setXYZ(ox, oy);
  }

  void preMultiplyPoint(T xi, T yi, T *xo, T *yo) const {
    *xo = m00_*xi + m10_*yi;
    *yo = m01_*xi + m11_*yi;
  }

  void preMultiplyPoint(const Point &ipoint, Point &opoint) const {
    opoint.x = m00_*ipoint.x + m10_*ipoint.y;
    opoint.y = m01_*ipoint.x + m11_*ipoint.y;
  }

  void preMultiplyVector(const Vector &ivector, Vector &ovector) const {
    T ix, iy;

    ivector.getXY(&ix, &iy);

    T ox = m00_*ix + m10_*iy;
    T oy = m01_*ix + m11_*iy;

    ovector.setXYZ(ox, oy);
  }

  void transpose() {
    swap(m10_, m01_);
  }

  CMatrix2x2T transposed() const {
    return CMatrix2x2T(m00_, m10_,
                       m01_, m11_);
  }

  bool invert(CMatrix2x2T &imatrix) const {
    T d = determinant();

    if (::fabs(d) < 1E-6)
      return false;

    T id = 1.0/d;

    imatrix.m00_ =  m11_*id;
    imatrix.m01_ = -m01_*id;
    imatrix.m10_ = -m10_*id;
    imatrix.m11_ =  m00_*id;

    return true;
  }

  CMatrix2x2T inverse() const {
    CMatrix2x2T imatrix;

    if (! invert(imatrix))
      CTHROW("Divide by zero");

    return imatrix;
  }

  T determinant() const {
    return (m00_*m11_ - m01_*m10_);
  }

  void normalize() {
    T d = determinant();

    T id = 1.0/d;

    for (int i = 0; i < 4; ++i)
      m1_[i] *= id;
  }

  static CMatrix2x2T *newIdentityMatrix() {
    CMatrix2x2T *m = new CMatrix2x2T();

    m->setIdentity();

    return m;
  }

  static bool solveAXeqB(const CMatrix2x2T &a, CPoint2DT<T> &x,
                         const CPoint2DT<T> &b) {
    T det_a = a.determinant();

    if (::fabs(det_a) < 1E-6)
      return false;

    T idet_a = 1.0/det_a;

    CMatrix2x2T t(a);

    t.setColumn(0, b.x, b.y);

    T det_t = t.determinant();

    x.x = det_t*idet_a;

    t = a;

    t.setColumn(1, b.x, b.y);

    det_t = t.determinant();

    x.y = det_t*idet_a;

    return true;
  }

  void zero() { memset(m1_, 0, sizeof(m1_)); }

  CMatrix2x2T &operator=(const CMatrix2x2T &a) {
    memcpy(m1_, a.m1_, sizeof(m1_));

    return *this;
  }

  CMatrix2x2T &operator+=(const CMatrix2x2T &b) {
    m00_ += b.m00_; m01_ += b.m01_,
    m10_ += b.m10_; m11_ += b.m11_;

    return *this;
  }

  CMatrix2x2T operator+(const CMatrix2x2T &b) {
    CMatrix2x2T c = *this;

    c += b;

    return c;
  }

  CMatrix2x2T &operator-=(const CMatrix2x2T &b) {
    m00_ -= b.m00_; m01_ -= b.m01_,
    m10_ -= b.m10_; m11_ -= b.m11_;

    return *this;
  }

  CMatrix2x2T operator-(const CMatrix2x2T &b) {
    CMatrix2x2T c = *this;

    c -= b;

    return c;
  }

  CMatrix2x2T &operator*=(const CMatrix2x2T &b) {
    CMatrix2x2T a;

    memcpy(a.m1_, m1_, sizeof(m1_));

    m00_ = a.m00_*b.m00_ + a.m01_*b.m10_;
    m01_ = a.m00_*b.m01_ + a.m01_*b.m11_;

    m10_ = a.m10_*b.m00_ + a.m11_*b.m10_;
    m11_ = a.m10_*b.m01_ + a.m11_*b.m11_;

    return *this;
  }

  CMatrix2x2T operator*(const CMatrix2x2T &b) {
    CMatrix2x2T c = *this;

    c *= b;

    return c;
  }

  CMatrix2x2T &operator*=(T s) {
    CMatrix2x2T a = *this;

    m00_ = a.m00_*s; m01_ = a.m01_*s;
    m10_ = a.m10_*s; m11_ = a.m11_*s;

    return *this;
  }

  CMatrix2x2T operator*(T s) {
    CMatrix2x2T c = *this;

    c *= s;

    return c;
  }

  friend Point operator*(const CMatrix2x2T &m, const Point &p) {
    Point p1;

    m.multiplyPoint(p, p1);

    return p1;
  }

  friend Point operator*(const Point &p, const CMatrix2x2T &m) {
    Point p1;

    m.preMultiplyPoint(p, p1);

    return p1;
  }

  friend Vector operator*(const CMatrix2x2T &m, const Vector &v) {
    Vector v1;

    m.multiplyVector(v, v1);

    return v1;
  }

  friend Vector operator*(const Vector &v, const CMatrix2x2T &m) {
    Vector v1;

    m.preMultiplyVector(v, v1);

    return v1;
  }

  CMatrix2x2T &operator/=(const CMatrix2x2T &b) {
    CMatrix2x2T bi;

    if (! b.invert(bi)) {
      CTHROW("Divide by zero");
      return *this;
    }

    return (*this) *= bi;
  }

  CMatrix2x2T operator/(const CMatrix2x2T &b) {
    CMatrix2x2T c = *this;

    c /= b;

    return c;
  }

  void setValue(unsigned int i, T value) {
    m1_[i] = value;
  }

  void setValue(unsigned int i, unsigned int j, T value) {
    m2_[i][j] = value;
  }

  T getValue(unsigned int i) { return m1_[i]; }

  T getValue(unsigned int i, unsigned int j) const {
    return m2_[i][j];
  }

  T operator[](unsigned int i) { return m1_[i]; }

  const T &operator[](unsigned int i) const { return m1_[i]; }

  void print(ostream &os) const {
    os << "(" << m00_ << "," << m01_ << ")" << endl;
    os << "(" << m10_ << "," << m11_ << ")" << endl;
  }

  friend ostream &operator<<(ostream &os, const CMatrix2x2T &matrix) {
    matrix.print(os);

    return os;
  }

 private:
  static T calcDeterminant(T m00, T m01, T m10, T m11) {
    return m00*m11 - m01*m10;
  }
};

typedef CMatrix2x2T<double> CMatrix2x2;

#endif
