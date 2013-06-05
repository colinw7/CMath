#ifndef CMATRIX_2D_H
#define CMATRIX_2D_H

#include <CMathGen.h>
#include <CPoint2D.h>
#include <CVector2D.h>
#include <CMatrixType.h>
#include <CThrow.h>

/* / a b tx \ */
/* | c d ty | */
/* \ 0 0 1  / */

/* / m00 m01 m02 \ */
/* | m10 m11 m12 | */
/* \ m20 m21 m22 / */

template<typename T>
class CMatrix2DT {
 private:
  typedef CPoint2DT<T>  Point;
  typedef CVector2DT<T> Vector;

 private:
  T m00_, m01_, m02_;
  T m10_, m11_, m12_;
  T m20_, m21_, m22_;

 public:
  // constructor/destructor
  CMatrix2DT() :
   m00_(0.0), m01_(0.0), m02_(0.0),
   m10_(0.0), m11_(0.0), m12_(0.0),
   m20_(0.0), m21_(0.0), m22_(0.0) {
  }

 ~CMatrix2DT() { }

  explicit
  CMatrix2DT(CMatrixType type) :
   m00_(0.0), m01_(0.0), m02_(0.0),
   m10_(0.0), m11_(0.0), m12_(0.0),
   m20_(0.0), m21_(0.0), m22_(0.0) {
    if (type == CMATRIX_TYPE_IDENTITY)
      setIdentity();
    else
      CTHROW("Bad Matrix Type");
  }

  CMatrix2DT(T a, T b, T c, T d) :
   m00_(a  ), m01_(b  ), m02_(0.0),
   m10_(c  ), m11_(d  ), m12_(0.0),
   m20_(0.0), m21_(0.0), m22_(0.0) {
    setBottomIdentity();
  }

  CMatrix2DT(T a, T b, T c, T d, T tx, T ty) :
   m00_(a), m01_(b), m02_(tx), m10_(c), m11_(d), m12_(ty) {
    setBottomIdentity();
  }

  CMatrix2DT(const T *m, int n) :
   m00_(0.0), m01_(0.0), m02_(0.0),
   m10_(0.0), m11_(0.0), m12_(0.0),
   m20_(0.0), m21_(0.0), m22_(0.0) {
    if      (n == 4)
      setValues(m[0], m[1], m[2], m[3]);
    else if (n == 6)
      setValues(m[0], m[1], m[2], m[3], m[4], m[5]);
    else
     CTHROW("Invalid size");
  }

  CMatrix2DT(const Vector &v0, const Vector &v1) :
   m00_(v0.getX()), m01_(v1.getX()), m02_(0.0),
   m10_(v0.getY()), m11_(v1.getY()), m12_(0.0),
   m20_(0.0      ), m21_(0.0      ), m22_(0.0) {
    setOuterIdentity ();
    setBottomIdentity();
  }

  CMatrix2DT *dup() const {
    return new CMatrix2DT(*this);
  }

  //------

  // copy operations
  CMatrix2DT(const CMatrix2DT &a) :
   m00_(a.m00_), m01_(a.m01_), m02_(a.m02_),
   m10_(a.m10_), m11_(a.m11_), m12_(a.m12_),
   m20_(a.m20_), m21_(a.m21_), m22_(a.m22_) {
  }

  const CMatrix2DT &operator=(const CMatrix2DT &a) {
    memcpy(&m00_, &a.m00_, 9*sizeof(T));

    return *this;
  }

  //------

  // output
  void print(std::ostream &os) const {
    os << "((" << m00_ << "," << m01_ << "," << m02_ << ")" <<
          " (" << m10_ << "," << m11_ << "," << m12_ << ")" <<
          " (" << m20_ << "," << m21_ << "," << m22_ << "))";
  }

  friend std::ostream &operator<<(std::ostream &os, const CMatrix2DT &matrix) {
    matrix.print(os);

    return os;
  }

  //------

  static CMatrix2DT identity() {
    return CMatrix2DT(CMATRIX_TYPE_IDENTITY);
  }

  static CMatrix2DT translation(const Point &point) {
    return translation(point.x, point.y);
  }

  static CMatrix2DT translation(T tx, T ty) {
    CMatrix2DT matrix;

    matrix.setTranslation(tx, ty);

    return matrix;
  }

  static CMatrix2DT scale(const Point &point) {
    return scale(point.x, point.y);
  }

  static CMatrix2DT scale(T s) {
    CMatrix2DT matrix;

    matrix.setScale(s, s);

    return matrix;
  }

  static CMatrix2DT scale(T sx, T sy) {
    CMatrix2DT matrix;

    matrix.setScale(sx, sy);

    return matrix;
  }

  static CMatrix2DT rotation(T a) {
    CMatrix2DT matrix;

    matrix.setRotation(a);

    return matrix;
  }

  static CMatrix2DT skew(T sx, T sy) {
    CMatrix2DT matrix;

    matrix.setSkew(sx, sy);

    return matrix;
  }

  static CMatrix2DT reflection(T a) {
    CMatrix2DT matrix;

    matrix.setReflection(a);

    return matrix;
  }

  static CMatrix2DT reflection(T dx, T dy) {
    CMatrix2DT matrix;

    matrix.setReflection(dx, dy);

    return matrix;
  }

  void setIdentity() {
    setInnerIdentity ();
    setOuterIdentity ();
    setBottomIdentity();
  }

  void setTranslation(T tx, T ty) {
    setInnerIdentity ();
    setOuterTranslate(tx, ty);
    setBottomIdentity();
  }

  void setTranslation(const Vector &vector) {
    setInnerIdentity ();
    setOuterTranslate(vector.getX(), vector.getY());
    setBottomIdentity();
  }

  void setScale(T s) {
    setInnerScale    (s, s);
    setOuterIdentity ();
    setBottomIdentity();
  }

  void setScale(T sx, T sy) {
    setInnerScale    (sx, sy);
    setOuterIdentity ();
    setBottomIdentity();
  }

  void setScaleTranslation(T sx, T sy, T tx, T ty) {
    setInnerScale    (sx, sy);
    setOuterTranslate(tx, ty);
    setBottomIdentity();
  }

  void setScaleTranslation(T s, T tx, T ty) {
    setInnerScale    (s, s);
    setOuterTranslate(tx, ty);
    setBottomIdentity();
  }

  void setRotation(T a) {
    setInnerRotation (a);
    setOuterIdentity ();
    setBottomIdentity();
  }

  void setRotationTranslation(T a, T tx, T ty) {
    setInnerRotation (a);
    setOuterTranslate(tx, ty);
    setBottomIdentity();
  }

  void setReflection(T a) {
    setUnitInnerReflection(cos(a), sin(a));
    setOuterIdentity      ();
    setBottomIdentity     ();
  }

  void setReflection(T dx, T dy) {
    setInnerReflection(dx, dy);
    setOuterIdentity  ();
    setBottomIdentity ();
  }

  void setSkew(T x, T y) {
    setInnerSkew     (x, y);
    setOuterIdentity ();
    setBottomIdentity();
  }

  void setSkewX(T a) {
    setInnerSkewX    (a);
    setOuterIdentity ();
    setBottomIdentity();
  }

  void setSkewY(T a) {
    setInnerSkewY    (a);
    setOuterIdentity ();
    setBottomIdentity();
  }

  void setValues(T a, T b, T c, T d) {
    m00_ = a, m01_ = b;
    m10_ = c, m11_ = d;

    setOuterIdentity ();
    setBottomIdentity();
  }

  void setValues(T a, T b, T c, T d, T tx, T ty) {
    m00_ = a, m01_ = b, m02_ = tx;
    m10_ = c, m11_ = d, m12_ = ty;

    setBottomIdentity();
  }

  void getValues(T *a, T *b, T *c, T *d) const {
    if (a) *a = m00_;
    if (b) *b = m01_;
    if (c) *c = m10_;
    if (d) *d = m11_;
  }

  void getValues(T *a, T *b, T *c, T *d, T *tx, T *ty) const {
    if (a ) *a  = m00_;
    if (b ) *b  = m01_;
    if (c ) *c  = m10_;
    if (d ) *d  = m11_;
    if (tx) *tx = m02_;
    if (ty) *ty = m12_;
  }

  void getValues(T *v, int n) const {
    if      (n == 4) {
      v[0] = m00_; v[1] = m01_;
      v[2] = m10_; v[3] = m11_;
    }
    else if (n == 6) {
      v[0] = m00_; v[1] = m01_;
      v[2] = m10_; v[3] = m11_;
      v[4] = m02_; v[5] = m12_;
    }
    else if (n == 9) {
      v[0] = m00_; v[1] = m01_; v[2] = m02_;
      v[3] = m10_; v[4] = m11_; v[5] = m12_;
      v[6] = m20_; v[7] = m21_; v[8] = m22_;
    }
    else
      CTHROW("Invalid Size");
  }

  //---------

  void setColumn(int c, T x, T y) {
    switch (c) {
      case 0: m00_ = x; m10_ = y; break;
      case 1: m01_ = x; m11_ = y; break;
    }
  }

  void setColumn(int c, const Point &point) {
    switch (c) {
      case 0: m00_ = point.x; m10_ = point.y; break;
      case 1: m01_ = point.x; m11_ = point.y; break;
    }
  }

  void setColumn(int c, const Vector &vector) {
    switch (c) {
      case 0: vector.getXY(&m00_, &m10_); break;
      case 1: vector.getXY(&m01_, &m11_); break;
    }
  }

  void setColumns(Vector &u, Vector &v) {
    setColumn(0, u);
    setColumn(1, v);
  }

  void getColumn(int c, T *x, T *y) const {
    switch (c) {
      case 0: if (x) *x = m00_; if (y) *y = m10_; break;
      case 1: if (x) *x = m01_; if (y) *y = m11_; break;
    }
  }

  void getColumn(int c, Vector &vector) {
    switch (c) {
      case 0: vector = Vector(m00_, m10_); break;
      case 1: vector = Vector(m01_, m11_); break;
    }
  }

  void getColumns(Vector &u, Vector &v) {
    getColumn(0, u);
    getColumn(1, v);
  }

  //------

  void setRow(int r, T x, T y) {
    switch (r) {
      case 0: m00_ = x; m01_ = y; break;
      case 1: m10_ = x; m11_ = y; break;
    }
  }

  void setRow(int r, const Point &point) {
    switch (r) {
      case 0: m00_ = point.x; m01_ = point.y; break;
      case 1: m10_ = point.x; m11_ = point.y; break;
    }
  }

  void setRow(int r, const Vector &vector) {
    switch (r) {
      case 0: vector.getXY(&m00_, &m01_); break;
      case 1: vector.getXY(&m10_, &m11_); break;
    }
  }

  void getRow(int r, T *x, T *y) const {
    switch (r) {
      case 0: if (x) *x = m00_; if (y) *y = m01_; break;
      case 1: if (x) *x = m10_; if (y) *y = m11_; break;
    }
  }

  //------

  void multiplyPoint(T xi, T yi, T *xo, T *yo) const {
    *xo = m00_*xi + m01_*yi + m02_;
    *yo = m10_*xi + m11_*yi + m12_;
  }

  void multiplyPoint(const Point &point1, Point &point2) const {
    point2.x = m00_*point1.x + m01_*point1.y + m02_;
    point2.y = m10_*point1.x + m11_*point1.y + m12_;
  }

  void multiplyVector(const Vector &ivector, Vector &ovector) const {
    T ix, iy;

    ivector.getXY(&ix, &iy);

    T ox = m00_*ix + m01_*iy;
    T oy = m10_*ix + m11_*iy;

    ovector.setXY(ox, oy);
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

    ovector.setXY(ox, oy);
  }

  const CMatrix2DT &translate(T x, T y) {
    m02_ += x;
    m12_ += y;

    return *this;
  }

  const CMatrix2DT &transpose() {
    swap(m10_, m01_);
    swap(m20_, m02_);
    swap(m21_, m12_);

    return *this;
  }

  CMatrix2DT transposed() const {
    return CMatrix2DT(m00_, m10_, m20_, m01_, m11_, m21_, m02_, m12_, m22_);
  }

  bool invert(CMatrix2DT &imatrix) const {
    T d = determinant();

    if (::fabs(d) == 0.0)
      return false;

    T id = 1.0/d;

    imatrix.m00_ =  id*calcDeterminant(m11_, m12_, m21_, m22_);
    imatrix.m10_ = -id*calcDeterminant(m10_, m12_, m20_, m22_);
    imatrix.m20_ =  id*calcDeterminant(m10_, m11_, m20_, m21_);

    imatrix.m01_ = -id*calcDeterminant(m01_, m02_, m21_, m22_);
    imatrix.m11_ =  id*calcDeterminant(m00_, m02_, m20_, m22_);
    imatrix.m21_ = -id*calcDeterminant(m00_, m01_, m20_, m21_);

    imatrix.m02_ =  id*calcDeterminant(m01_, m02_, m11_, m12_);
    imatrix.m12_ = -id*calcDeterminant(m00_, m02_, m10_, m12_);
    imatrix.m22_ =  id*calcDeterminant(m00_, m01_, m10_, m11_);

    return true;
  }

  CMatrix2DT inverse() const {
    CMatrix2DT imatrix;

    if (! invert(imatrix))
      CTHROW("Divide by zero");

    return imatrix;
  }

  T determinant() const {
    return (m00_*calcDeterminant(m11_, m12_, m21_, m22_) -
            m01_*calcDeterminant(m10_, m12_, m20_, m22_) +
            m02_*calcDeterminant(m10_, m11_, m20_, m21_));
  }

  const CMatrix2DT &normalize() {
    T d = determinant();

    T id = 1.0/d;

    m00_ *= id;
    m01_ *= id;
    m10_ *= id;
    m11_ *= id;

    return *this;
  }

  void setTransform(T xmin1, T ymin1, T xmax1, T ymax1, T xmin2, T ymin2, T xmax2, T ymax2) {
    T sx = (xmax2 - xmin2)/(xmax1 - xmin1);
    T sy = (ymax2 - ymin2)/(ymax1 - ymin1);

    T tx = -xmin1*sx + xmin2;
    T ty = -ymin1*sy + ymin2;

    setInnerScale    (sx, sy);
    setOuterTranslate(tx, ty);
    setBottomIdentity();
  }

  static CMatrix2DT *newIdentityMatrix() {
    CMatrix2DT *m = new CMatrix2DT();

    m->setIdentity();

    return m;
  }

  static bool solveAXeqB(const CMatrix2DT &a, CPoint2DT<T> &x, const CPoint2DT<T> &b) {
    T det_a = a.determinant();

    if (::fabs(det_a) <= 0.0)
      return false;

    T idet_a = 1.0/det_a;

    CMatrix2DT t(a);

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
    return REAL_EQ(m00_, 1.0) && REAL_EQ(m11_, 1.0) && REAL_EQ(m22_, 1.0) &&
           REAL_EQ(m10_, 0.0) && REAL_EQ(m01_, 0.0) &&
           REAL_EQ(m20_, 0.0) && REAL_EQ(m02_, 0.0) &&
           REAL_EQ(m12_, 0.0) && REAL_EQ(m21_, 0.0);
  }

  double getAngle() const {
    return asin(m10_);
  }

  //------

  void getSize(double *sx, double *sy) const {
    *sx = fabs(m00_ + m01_);
    *sy = fabs(m10_ + m11_);
  }

  //------

  void zero() { memset(&m00_, 0, 9*sizeof(T)); }

  const CMatrix2DT &operator+=(const CMatrix2DT &b) {
    m00_ += b.m00_; m01_ += b.m01_; m02_ += b.m02_;
    m10_ += b.m10_; m11_ += b.m11_; m12_ += b.m12_;
    m20_ += b.m20_; m21_ += b.m21_; m22_ += b.m22_;

    return *this;
  }

  CMatrix2DT operator+(const CMatrix2DT &b) const {
    CMatrix2DT c = *this;

    c += b;

    return c;
  }

  const CMatrix2DT &operator-=(const CMatrix2DT &b) {
    m00_ -= b.m00_; m01_ -= b.m01_, m02_ -= b.m02_;
    m10_ -= b.m10_; m11_ -= b.m11_; m12_ -= b.m12_;
    m20_ -= b.m20_; m21_ -= b.m21_; m22_ -= b.m22_;

    return *this;
  }

  CMatrix2DT operator-(const CMatrix2DT &b) const {
    CMatrix2DT c = *this;

    c -= b;

    return c;
  }

  const CMatrix2DT &operator*=(const CMatrix2DT &b) {
    CMatrix2DT a = *this;

    m00_ = a.m00_*b.m00_ + a.m01_*b.m10_ + a.m02_*b.m20_;
    m01_ = a.m00_*b.m01_ + a.m01_*b.m11_ + a.m02_*b.m21_;
    m02_ = a.m00_*b.m02_ + a.m01_*b.m12_ + a.m02_*b.m22_;

    m10_ = a.m10_*b.m00_ + a.m11_*b.m10_ + a.m12_*b.m20_;
    m11_ = a.m10_*b.m01_ + a.m11_*b.m11_ + a.m12_*b.m21_;
    m12_ = a.m10_*b.m02_ + a.m11_*b.m12_ + a.m12_*b.m22_;

    m20_ = a.m20_*b.m00_ + a.m21_*b.m10_ + a.m22_*b.m20_;
    m21_ = a.m20_*b.m01_ + a.m21_*b.m11_ + a.m22_*b.m21_;
    m22_ = a.m20_*b.m02_ + a.m21_*b.m12_ + a.m22_*b.m22_;

    return *this;
  }

  CMatrix2DT operator*(const CMatrix2DT &b) const {
    CMatrix2DT c = *this;

    c *= b;

    return c;
  }

  const CMatrix2DT &operator*=(T s) {
    CMatrix2DT a = *this;

    m00_ = a.m00_*s; m01_ = a.m01_*s; m02_ = a.m02_*s;
    m10_ = a.m10_*s; m11_ = a.m11_*s; m12_ = a.m12_*s;
    m20_ = a.m20_*s; m21_ = a.m21_*s; m22_ = a.m22_*s;

    return *this;
  }

  CMatrix2DT operator*(T s) const {
    CMatrix2DT c = *this;

    c *= s;

    return c;
  }

  friend Point operator*(const CMatrix2DT &m, const Point &p) {
    Point p1;

    m.multiplyPoint(p, p1);

    return p1;
  }

  friend Point operator*(const Point &p, const CMatrix2DT &m) {
    Point p1;

    m.preMultiplyPoint(p, p1);

    return p1;
  }

  friend Vector operator*(const CMatrix2DT &m, const Vector &v) {
    Vector v1;

    m.multiplyVector(v, v1);

    return v1;
  }

  friend Vector operator*(const Vector &v, const CMatrix2DT &m) {
    Vector v1;

    m.preMultiplyVector(v, v1);

    return v1;
  }

  const CMatrix2DT &operator/=(const CMatrix2DT &b) {
    CMatrix2DT bi;

    if (! b.invert(bi)) {
      CTHROW("Divide by zero");
      return *this;
    }

    return (*this) *= bi;
  }

  CMatrix2DT operator/(const CMatrix2DT &b) const {
    CMatrix2DT c = *this;

    c /= b;

    return c;
  }

  //------

  void setValue(uint i, T value) {
    (&m00_)[i] = value;
  }

  void setValue(uint i, uint j, T value) {
    assert(i < 3 && j < 3);

    T &m = (&m00_)[3*j + i];

    m = value;
  }

  const T &getValue(uint i) const {
    assert(i < 9);

    return (&m00_)[i];
  }

  const T &getValue(uint i, uint j) const {
    assert(i < 3 && j < 3);

    const T &m = (&m00_)[3*j + i];

    return m;
  }

  const T &operator[](uint i) const {
    assert(i < 9);

    return (&m00_)[i];
  }

  T &operator[](uint i) {
    assert(i < 9);

    return (&m00_)[i];
  }

  //------

  void setInnerIdentity() {
    m00_ = 1.0, m01_ = 0.0;
    m10_ = 0.0, m11_ = 1.0;
  }

  void setInnerScale(T sx, T sy) {
    m00_ = sx , m01_ = 0.0;
    m10_ = 0.0, m11_ = sy ;
  }

  void setInnerRotation(T a) {
    T c = ::cos(a);
    T s = ::sin(a);

    m00_ =  c, m01_ = -s;
    m10_ =  s, m11_ =  c;
  }

  void setInnerSkew(T x, T y) {
    T tx = ::tan(x);
    T ty = ::tan(y);

    m00_ = 1 , m01_ = tx;
    m10_ = ty, m11_ = 1 ;
  }

  void setInnerSkewX(T x) {
    T tx = ::tan(x);

    m00_ = 1, m01_ = tx;
    m10_ = 0, m11_ = 1 ;
  }

  void setInnerSkewY(T y) {
    T ty = ::tan(y);

    m00_ = 1 , m01_ = 0;
    m10_ = ty, m11_ = 1;
  }

  void setInnerReflection(T dx, T dy) {
    double l = sqrt(dx*dx + dy*dy);

    setUnitInnerReflection(dx/l, dy/l);
  }

  void setUnitInnerReflection(T dx, T dy) {
    //m00_ = (dx*dx - dy*dy)/l; m01_ = 2*dx*dy/l;
    //m10_ = m01_             ; m11_ = -m00_;

    m00_ = 2.0*dx*dx - 1.0; m01_ = 2.0*dx*dy;
    m10_ = m01_           ; m11_ = 2.0*dy*dy - 1.0;
  }

  void setOuterIdentity() {
    m02_ = 0.0; m12_ = 0.0;
  }

  void setOuterTranslate(T tx, T ty) {
    m02_ = tx; m12_ = ty;
  }

  void setBottomIdentity() {
    m20_ = 0.0, m21_ = 0.0, m22_ = 1.0;
  }

 private:
  static T calcDeterminant(T m00, T m01, T m10, T m11) {
    return m00*m11 - m01*m10;
  }
};

typedef CMatrix2DT<double> CMatrix2D;

#endif
