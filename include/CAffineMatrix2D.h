#ifndef CMATRIX_2D_H
#define CMATRIX_2D_H

#include <CMathGen.h>
#include <CPoint2D.h>
#include <CVector2D.h>
#include <CMatrixType.h>
#include <CThrow.h>

/* / a b tx \ */
/* | c d ty | */

template<typename T>
class CAffineMatrix2DT {
 private:
  typedef CPoint2DT<T>  Point;
  typedef CVector2DT<T> Vector;

 private:
  T a_, b_, c_, d_;
  T tx_, ty_;

 public:
  // constructor/destructor
  CAffineMatrix2DT() :
   a_(1.0), b_(0.0), c_(0.0), d_(1.0), tx_(0.0), ty_(0.0) {
  }

 ~CAffineMatrix2DT() { }

  CAffineMatrix2DT(T a, T b, T c, T d) :
   a_(a), b_(b), c_(c), d_(d), tx_(0.0), ty_(0.0) {
  }

  CAffineMatrix2DT(T a, T b, T c, T d, T tx, T ty) :
   a_(a), b_(b), c_(c), d_(d), tx_(tx), ty_(ty) {
  }

  CAffineMatrix2DT(const T *m, int n) :
   a_(0.0), b_(0.0), c_(0.0), d_(0.0), tx_(0.0), ty_(0.0) {
    if      (n == 4)
      setValues(m[0], m[1], m[2], m[3]);
    else if (n == 6)
      setValues(m[0], m[1], m[2], m[3], m[4], m[5]);
    else
     CTHROW("Invalid size");
  }

  CAffineMatrix2DT(const Vector &v0, const Vector &v1) :
   a_(v0.getX()), b_(v1.getX()), c_(v0.getY()), d_(v1.getY()), tx_(0.0), ty_(0.0) {
  }

  CAffineMatrix2DT *dup() const {
    return new CAffineMatrix2DT(*this);
  }

  //------

  // copy operations
  CAffineMatrix2DT(const CAffineMatrix2DT &a) :
   a_(a.a_), b_(a.b_), c_(a.c_), d_(a.d_), tx_(a.tx_), ty_(a.ty_) {
  }

  const CAffineMatrix2DT &operator=(const CAffineMatrix2DT &a) {
    memcpy(&a_, &a.a_, 6*sizeof(T));

    return *this;
  }

  //------

  // output
  void print(std::ostream &os) const {
    os << "(" << a_ << "," << b_ << "," << c_ << "," << d_ << "," << tx_ << "," << ty_ << "))";
  }

  friend std::ostream &operator<<(std::ostream &os, const CAffineMatrix2DT &matrix) {
    matrix.print(os);

    return os;
  }

  //------

  static CAffineMatrix2DT identity() {
    return CAffineMatrix2DT();
  }

  static CAffineMatrix2DT translation(const Point &point) {
    return translation(point.x, point.y);
  }

  static CAffineMatrix2DT translation(T tx, T ty) {
    CAffineMatrix2DT matrix;

    matrix.setTranslation(tx, ty);

    return matrix;
  }

  static CAffineMatrix2DT scale(const Point &point) {
    return scale(point.x, point.y);
  }

  static CAffineMatrix2DT scale(T s) {
    CAffineMatrix2DT matrix;

    matrix.setScale(s, s);

    return matrix;
  }

  static CAffineMatrix2DT scale(T sx, T sy) {
    CAffineMatrix2DT matrix;

    matrix.setScale(sx, sy);

    return matrix;
  }

  static CAffineMatrix2DT rotation(T a) {
    CAffineMatrix2DT matrix;

    matrix.setRotation(a);

    return matrix;
  }

  static CAffineMatrix2DT skew(T sx, T sy) {
    CAffineMatrix2DT matrix;

    matrix.setSkew(sx, sy);

    return matrix;
  }

  static CAffineMatrix2DT reflection(T a) {
    CAffineMatrix2DT matrix;

    matrix.setReflection(a);

    return matrix;
  }

  static CAffineMatrix2DT reflection(T dx, T dy) {
    CAffineMatrix2DT matrix;

    matrix.setReflection(dx, dy);

    return matrix;
  }

  void setIdentity() {
    setInnerIdentity ();
    setOuterTranslate(0.0, 0.0);
  }

  void setTranslation(T tx, T ty) {
    setInnerIdentity ();
    setOuterTranslate(tx, ty);
  }

  void setTranslation(const Vector &vector) {
    setInnerIdentity ();
    setOuterTranslate(vector.getX(), vector.getY());
  }

  void setScale(T s) {
    setInnerScale    (s, s);
    setOuterTranslate(0.0, 0.0);
  }

  void setScale(T sx, T sy) {
    setInnerScale    (sx, sy);
    setOuterTranslate(0.0, 0.0);
  }

  void setScaleTranslation(T sx, T sy, T tx, T ty) {
    setInnerScale    (sx, sy);
    setOuterTranslate(tx, ty);
  }

  void setScaleTranslation(T s, T tx, T ty) {
    setInnerScale    (s, s);
    setOuterTranslate(tx, ty);
  }

  void setRotation(T a) {
    setInnerRotation (a);
    setOuterTranslate(0.0, 0.0);
  }

  void setRotationTranslation(T a, T tx, T ty) {
    setInnerRotation (a);
    setOuterTranslate(tx, ty);
  }

  void setReflection(T a) {
    setUnitInnerReflection(cos(a), sin(a));
    setOuterTranslate     (0.0, 0.0);
  }

  void setReflection(T dx, T dy) {
    setInnerReflection(dx, dy);
    setOuterTranslate (0.0, 0.0);
  }

  void setSkew(T x, T y) {
    setInnerSkew     (x, y);
    setOuterTranslate(0.0, 0.0);
  }

  void setSkewX(T a) {
    setInnerSkewX    (a);
    setOuterTranslate(0.0, 0.0);
  }

  void setSkewY(T a) {
    setInnerSkewY    (a);
    setOuterTranslate(0.0, 0.0);
  }

  void setValues(T a, T b, T c, T d) {
    a_ = a, b_ = b; tx_ = 0.0;
    c_ = c, d_ = d; ty_ = 0.0;
  }

  void setValues(T a, T b, T c, T d, T tx, T ty) {
    a_ = a, b_ = b, tx_ = tx;
    c_ = c, d_ = d, ty_ = ty;
  }

  void getValues(T *a, T *b, T *c, T *d) const {
    if (a) *a = a_;
    if (b) *b = b_;
    if (c) *c = c_;
    if (d) *d = d_;
  }

  void getValues(T *a, T *b, T *c, T *d, T *tx, T *ty) const {
    if (a ) *a  = a_;
    if (b ) *b  = b_;
    if (c ) *c  = c_;
    if (d ) *d  = d_;
    if (tx) *tx = tx_;
    if (ty) *ty = ty_;
  }

  void getValues(T *v, int n) const {
    if      (n == 4) {
      v[0] = a_; v[1] = b_;
      v[2] = c_; v[3] = d_;
    }
    else if (n == 6) {
      v[0] = a_ ; v[1] = b_ ;
      v[2] = c_ ; v[3] = d_ ;
      v[4] = tx_; v[5] = ty_;
    }
    else
      CTHROW("Invalid Size");
  }

  //---------

  void setColumn(int c, T x, T y) {
    switch (c) {
      case 0: a_ = x; c_ = y; break;
      case 1: b_ = x; d_ = y; break;
    }
  }

  void setColumn(int c, const Point &point) {
    switch (c) {
      case 0: a_ = point.x; c_ = point.y; break;
      case 1: b_ = point.x; d_ = point.y; break;
    }
  }

  void setColumn(int c, const Vector &vector) {
    switch (c) {
      case 0: vector.getXY(&a_, &c_); break;
      case 1: vector.getXY(&b_, &d_); break;
    }
  }

  void setColumns(Vector &u, Vector &v) {
    setColumn(0, u);
    setColumn(1, v);
  }

  void getColumn(int c, T *x, T *y) const {
    switch (c) {
      case 0: if (x) *x = a_; if (y) *y = c_; break;
      case 1: if (x) *x = b_; if (y) *y = d_; break;
    }
  }

  void getColumn(int c, Vector &vector) {
    switch (c) {
      case 0: vector = Vector(a_, c_); break;
      case 1: vector = Vector(b_, d_); break;
    }
  }

  void getColumns(Vector &u, Vector &v) {
    getColumn(0, u);
    getColumn(1, v);
  }

  //------

  void setRow(int r, T x, T y) {
    switch (r) {
      case 0: a_ = x; b_ = y; break;
      case 1: c_ = x; d_ = y; break;
    }
  }

  void setRow(int r, const Point &point) {
    switch (r) {
      case 0: a_ = point.x; b_ = point.y; break;
      case 1: c_ = point.x; d_ = point.y; break;
    }
  }

  void setRow(int r, const Vector &vector) {
    switch (r) {
      case 0: vector.getXY(&a_, &b_); break;
      case 1: vector.getXY(&c_, &d_); break;
    }
  }

  void getRow(int r, T *x, T *y) const {
    switch (r) {
      case 0: if (x) *x = a_; if (y) *y = b_; break;
      case 1: if (x) *x = c_; if (y) *y = d_; break;
    }
  }

  //------

  void multiplyPoint(T xi, T yi, T *xo, T *yo) const {
    *xo = a_*xi + b_*yi + tx_;
    *yo = c_*xi + d_*yi + ty_;
  }

  void multiplyPoint(const Point &point1, Point &point2) const {
    point2.x = a_*point1.x + b_*point1.y + tx_;
    point2.y = c_*point1.x + d_*point1.y + ty_;
  }

  void multiplyVector(const Vector &ivector, Vector &ovector) const {
    T ix, iy;

    ivector.getXY(&ix, &iy);

    T ox = a_*ix + b_*iy;
    T oy = c_*ix + d_*iy;

    ovector.setXY(ox, oy);
  }

  void preMultiplyPoint(T xi, T yi, T *xo, T *yo) const {
    *xo = a_*xi + c_*yi;
    *yo = b_*xi + d_*yi;
  }

  void preMultiplyPoint(const Point &ipoint, Point &opoint) const {
    opoint.x = a_*ipoint.x + c_*ipoint.y;
    opoint.y = b_*ipoint.x + d_*ipoint.y;
  }

  void preMultiplyVector(const Vector &ivector, Vector &ovector) const {
    T ix, iy;

    ivector.getXY(&ix, &iy);

    T ox = a_*ix + c_*iy;
    T oy = b_*ix + d_*iy;

    ovector.setXY(ox, oy);
  }

  const CAffineMatrix2DT &translate(T x, T y) {
    tx_ += x;
    ty_ += y;

    return *this;
  }

  bool invert(CAffineMatrix2DT &imatrix) const {
    T d = determinant();

    if (::fabs(d) == 0.0)
      return false;

    T id = 1.0/d;

    imatrix.a_ =  id*d_;
    imatrix.d_ =  id*a_;

    imatrix.c_ = -id*c_;
    imatrix.b_ = -id*b_;

    imatrix.tx_ =  id*calcDeterminant(b_, tx_, d_, ty_);
    imatrix.ty_ = -id*calcDeterminant(a_, tx_, c_, ty_);

    return true;
  }

  CAffineMatrix2DT inverse() const {
    CAffineMatrix2DT imatrix;

    if (! invert(imatrix))
      CTHROW("Divide by zero");

    return imatrix;
  }

  T determinant() const {
    return (a_*d_ - b_*c_);
  }

  const CAffineMatrix2DT &normalize() {
    T d = determinant();

    T id = 1.0/d;

    a_ *= id;
    b_ *= id;
    c_ *= id;
    d_ *= id;

    return *this;
  }

  void setTransform(T xmin1, T ymin1, T xmax1, T ymax1, T xmin2, T ymin2, T xmax2, T ymax2) {
    T sx = (xmax2 - xmin2)/(xmax1 - xmin1);
    T sy = (ymax2 - ymin2)/(ymax1 - ymin1);

    T tx = -xmin1*sx + xmin2;
    T ty = -ymin1*sy + ymin2;

    setInnerScale    (sx, sy);
    setOuterTranslate(tx, ty);
  }

  static CAffineMatrix2DT *newIdentityMatrix() {
    CAffineMatrix2DT *m = new CAffineMatrix2DT();

    return m;
  }

  static bool solveAXeqB(const CAffineMatrix2DT &a, CPoint2DT<T> &x, const CPoint2DT<T> &b) {
    T det_a = a.determinant();

    if (::fabs(det_a) <= 0.0)
      return false;

    T idet_a = 1.0/det_a;

    CAffineMatrix2DT t(a);

    t.setColumn(0, b.x, b.y);

    T det_t = t.determinant();

    x.x = det_t*idet_a;

    t = a;

    t.setColumn(1, b.x, b.y);

    det_t = t.determinant();

    x.y = det_t*idet_a;

    return true;
  }

  //------

  bool isIdentity() const {
    return REAL_EQ(a_ , 1.0) && REAL_EQ(d_ , 1.0) &&
           REAL_EQ(c_ , 0.0) && REAL_EQ(b_ , 0.0) &&
           REAL_EQ(tx_, 0.0) && REAL_EQ(ty_, 0.0);
  }

  double getAngle() const { return asin(c_); }

  //------

  void zero() { memset(&a_, 0, 6*sizeof(T)); }

  const CAffineMatrix2DT &operator+=(const CAffineMatrix2DT &b) {
    a_ += b.a_; b_ += b.b_; tx_ += b.tx_;
    c_ += b.c_; d_ += b.d_; ty_ += b.ty_;

    return *this;
  }

  CAffineMatrix2DT operator+(const CAffineMatrix2DT &b) const {
    CAffineMatrix2DT c = *this;

    c += b;

    return c;
  }

  const CAffineMatrix2DT &operator-=(const CAffineMatrix2DT &b) {
    a_ -= b.a_; b_ -= b.b_, tx_ -= b.tx_;
    c_ -= b.c_; d_ -= b.d_; ty_ -= b.ty_;

    return *this;
  }

  CAffineMatrix2DT operator-(const CAffineMatrix2DT &b) const {
    CAffineMatrix2DT c = *this;

    c -= b;

    return c;
  }

  const CAffineMatrix2DT &operator*=(const CAffineMatrix2DT &b) {
    CAffineMatrix2DT a = *this;

    a_  = a.a_*b.a_ + a.b_*b.c_;
    b_  = a.a_*b.b_ + a.b_*b.d_;
    c_  = a.c_*b.a_ + a.d_*b.c_;
    d_  = a.c_*b.b_ + a.d_*b.d_;

    tx_ = a.a_*b.tx_ + a.b_*b.ty_ + a.tx_;
    ty_ = a.c_*b.tx_ + a.d_*b.ty_ + a.ty_;

    return *this;
  }

  CAffineMatrix2DT operator*(const CAffineMatrix2DT &b) const {
    CAffineMatrix2DT c = *this;

    c *= b;

    return c;
  }

  const CAffineMatrix2DT &operator*=(T s) {
    CAffineMatrix2DT a = *this;

    a_ = a.a_*s; b_ = a.b_*s; tx_ = a.tx_*s;
    c_ = a.c_*s; d_ = a.d_*s; ty_ = a.ty_*s;

    return *this;
  }

  CAffineMatrix2DT operator*(T s) const {
    CAffineMatrix2DT c = *this;

    c *= s;

    return c;
  }

  friend Point operator*(const CAffineMatrix2DT &m, const Point &p) {
    Point p1;

    m.multiplyPoint(p, p1);

    return p1;
  }

  friend Point operator*(const Point &p, const CAffineMatrix2DT &m) {
    Point p1;

    m.preMultiplyPoint(p, p1);

    return p1;
  }

  friend Vector operator*(const CAffineMatrix2DT &m, const Vector &v) {
    Vector v1;

    m.multiplyVector(v, v1);

    return v1;
  }

  friend Vector operator*(const Vector &v, const CAffineMatrix2DT &m) {
    Vector v1;

    m.preMultiplyVector(v, v1);

    return v1;
  }

  const CAffineMatrix2DT &operator/=(const CAffineMatrix2DT &b) {
    CAffineMatrix2DT bi;

    if (! b.invert(bi)) {
      CTHROW("Divide by zero");
      return *this;
    }

    return (*this) *= bi;
  }

  CAffineMatrix2DT operator/(const CAffineMatrix2DT &b) const {
    CAffineMatrix2DT c = *this;

    c /= b;

    return c;
  }

  //------

  void setValue(uint i, T value) {
    assert(i < 6);

    (&a_)[i] = value;
  }

  void setValue(uint i, uint j, T value) {
    assert(i < 2 && j < 2);

    T &m = (&a_)[2*j + i];

    m = value;
  }

  const T &getValue(uint i) const {
    assert(i < 6);

    return (&a_)[i];
  }

  const T &getValue(uint i, uint j) const {
    assert(i < 2 && j < 2);

    const T &m = (&a_)[2*j + i];

    return m;
  }

  const T &operator[](uint i) const {
    assert(i < 6);

    return (&a_)[i];
  }

  T &operator[](uint i) {
    assert(i < 6);

    return (&a_)[i];
  }

  //------

  void setInnerIdentity() {
    a_ = 1.0, b_ = 0.0;
    c_ = 0.0, d_ = 1.0;
  }

  void setInnerScale(T sx, T sy) {
    a_ = sx , b_ = 0.0;
    c_ = 0.0, d_ = sy ;
  }

  void setInnerRotation(T a) {
    T c = ::cos(a);
    T s = ::sin(a);

    a_ =  c, b_ = -s;
    c_ =  s, d_ =  c;
  }

  void setInnerSkew(T x, T y) {
    T tx = ::tan(x);
    T ty = ::tan(y);

    a_ = 1 , b_ = tx;
    c_ = ty, d_ = 1 ;
  }

  void setInnerSkewX(T x) {
    T tx = ::tan(x);

    a_ = 1, b_ = tx;
    c_ = 0, d_ = 1 ;
  }

  void setInnerSkewY(T y) {
    T ty = ::tan(y);

    a_ = 1 , b_ = 0;
    c_ = ty, d_ = 1;
  }

  void setInnerReflection(T dx, T dy) {
    double l = sqrt(dx*dx + dy*dy);

    setUnitInnerReflection(dx/l, dy/l);
  }

  void setUnitInnerReflection(T dx, T dy) {
  //a_ = (dx*dx - dy*dy)/l; b_ = 2*dx*dy/l;
  //c_ = b_               ; d_ = -a_;

    a_ = 2.0*dx*dx - 1.0; b_ = 2.0*dx*dy;
    c_ = b_             ; d_ = 2.0*dy*dy - 1.0;
  }

  void setOuterTranslate(T tx, T ty) {
    tx_ = tx; ty_ = ty;
  }

 private:
  static T calcDeterminant(T m00, T m01, T m10, T m11) {
    return m00*m11 - m01*m10;
  }
};

typedef CAffineMatrix2DT<double> CAffineMatrix2D;

#endif
