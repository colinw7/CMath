#ifndef CMATRIX_2D_H
#define CMATRIX_2D_H

#include <CPoint2D.h>
#include <CVector2D.h>

/* / a b tx \ */
/* | c d ty | */

class CAffineMatrix2D {
 private:
  typedef CPoint2D  Point;
  typedef CVector2D Vector;

 public:
  // constructor/destructor
  CAffineMatrix2D() :
   a_(1.0), b_(0.0), c_(0.0), d_(1.0), tx_(0.0), ty_(0.0) {
  }

 ~CAffineMatrix2D() { }

  CAffineMatrix2D(double a, double b, double c, double d) :
   a_(a), b_(b), c_(c), d_(d), tx_(0.0), ty_(0.0) {
  }

  CAffineMatrix2D(double a, double b, double c, double d, double tx, double ty) :
   a_(a), b_(b), c_(c), d_(d), tx_(tx), ty_(ty) {
  }

  CAffineMatrix2D(const double *m, int n) :
   a_(0.0), b_(0.0), c_(0.0), d_(0.0), tx_(0.0), ty_(0.0) {
    if      (n == 4)
      setValues(m[0], m[1], m[2], m[3]);
    else if (n == 6)
      setValues(m[0], m[1], m[2], m[3], m[4], m[5]);
    else
     assert(false && "Invalid size");
  }

  CAffineMatrix2D(const Vector &v0, const Vector &v1) :
   a_(v0.getX()), b_(v1.getX()), c_(v0.getY()), d_(v1.getY()), tx_(0.0), ty_(0.0) {
  }

  CAffineMatrix2D *dup() const {
    return new CAffineMatrix2D(*this);
  }

  //------

  // copy operations
  CAffineMatrix2D(const CAffineMatrix2D &a) :
   a_(a.a_), b_(a.b_), c_(a.c_), d_(a.d_), tx_(a.tx_), ty_(a.ty_) {
  }

  const CAffineMatrix2D &operator=(const CAffineMatrix2D &a) {
    memcpy(&a_, &a.a_, 6*sizeof(double));

    return *this;
  }

  //------

  // output
  void print(std::ostream &os) const {
    os << "(" << a_ << "," << b_ << "," << c_ << "," << d_ << "," << tx_ << "," << ty_ << "))";
  }

  friend std::ostream &operator<<(std::ostream &os, const CAffineMatrix2D &matrix) {
    matrix.print(os);

    return os;
  }

  //------

  static CAffineMatrix2D identity() {
    return CAffineMatrix2D();
  }

  static CAffineMatrix2D translation(const Point &point) {
    return translation(point.x, point.y);
  }

  static CAffineMatrix2D translation(double tx, double ty) {
    CAffineMatrix2D matrix;

    matrix.setTranslation(tx, ty);

    return matrix;
  }

  static CAffineMatrix2D scale(const Point &point) {
    return scale(point.x, point.y);
  }

  static CAffineMatrix2D scale(double s) {
    CAffineMatrix2D matrix;

    matrix.setScale(s, s);

    return matrix;
  }

  static CAffineMatrix2D scale(double sx, double sy) {
    CAffineMatrix2D matrix;

    matrix.setScale(sx, sy);

    return matrix;
  }

  static CAffineMatrix2D rotation(double a) {
    CAffineMatrix2D matrix;

    matrix.setRotation(a);

    return matrix;
  }

  static CAffineMatrix2D skew(double sx, double sy) {
    CAffineMatrix2D matrix;

    matrix.setSkew(sx, sy);

    return matrix;
  }

  static CAffineMatrix2D reflection(double a) {
    CAffineMatrix2D matrix;

    matrix.setReflection(a);

    return matrix;
  }

  static CAffineMatrix2D reflection(double dx, double dy) {
    CAffineMatrix2D matrix;

    matrix.setReflection(dx, dy);

    return matrix;
  }

  void setIdentity() {
    setInnerIdentity ();
    setOuterTranslate(0.0, 0.0);
  }

  void setTranslation(double tx, double ty) {
    setInnerIdentity ();
    setOuterTranslate(tx, ty);
  }

  void setTranslation(const Vector &vector) {
    setInnerIdentity ();
    setOuterTranslate(vector.getX(), vector.getY());
  }

  void setScale(double s) {
    setInnerScale    (s, s);
    setOuterTranslate(0.0, 0.0);
  }

  void setScale(double sx, double sy) {
    setInnerScale    (sx, sy);
    setOuterTranslate(0.0, 0.0);
  }

  void setScaleTranslation(double sx, double sy, double tx, double ty) {
    setInnerScale    (sx, sy);
    setOuterTranslate(tx, ty);
  }

  void setScaleTranslation(double s, double tx, double ty) {
    setInnerScale    (s, s);
    setOuterTranslate(tx, ty);
  }

  void setRotation(double a) {
    setInnerRotation (a);
    setOuterTranslate(0.0, 0.0);
  }

  void setRotationTranslation(double a, double tx, double ty) {
    setInnerRotation (a);
    setOuterTranslate(tx, ty);
  }

  void setReflection(double a) {
    setUnitInnerReflection(cos(a), sin(a));
    setOuterTranslate     (0.0, 0.0);
  }

  void setReflection(double dx, double dy) {
    setInnerReflection(dx, dy);
    setOuterTranslate (0.0, 0.0);
  }

  void setSkew(double x, double y) {
    setInnerSkew     (x, y);
    setOuterTranslate(0.0, 0.0);
  }

  void setSkewX(double a) {
    setInnerSkewX    (a);
    setOuterTranslate(0.0, 0.0);
  }

  void setSkewY(double a) {
    setInnerSkewY    (a);
    setOuterTranslate(0.0, 0.0);
  }

  void setValues(double a, double b, double c, double d) {
    a_ = a, b_ = b; tx_ = 0.0;
    c_ = c, d_ = d; ty_ = 0.0;
  }

  void setValues(double a, double b, double c, double d, double tx, double ty) {
    a_ = a, b_ = b, tx_ = tx;
    c_ = c, d_ = d, ty_ = ty;
  }

  void getValues(double *a, double *b, double *c, double *d) const {
    if (a) *a = a_;
    if (b) *b = b_;
    if (c) *c = c_;
    if (d) *d = d_;
  }

  void getValues(double *a, double *b, double *c, double *d, double *tx, double *ty) const {
    if (a ) *a  = a_;
    if (b ) *b  = b_;
    if (c ) *c  = c_;
    if (d ) *d  = d_;
    if (tx) *tx = tx_;
    if (ty) *ty = ty_;
  }

  void getValues(double *v, int n) const {
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
      assert(false && "Invalid Size");
  }

  //---------

  void setColumn(int c, double x, double y) {
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

  void getColumn(int c, double *x, double *y) const {
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

  void setRow(int r, double x, double y) {
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

  void getRow(int r, double *x, double *y) const {
    switch (r) {
      case 0: if (x) *x = a_; if (y) *y = b_; break;
      case 1: if (x) *x = c_; if (y) *y = d_; break;
    }
  }

  //------

  void multiplyPoint(double xi, double yi, double *xo, double *yo) const {
    *xo = a_*xi + b_*yi + tx_;
    *yo = c_*xi + d_*yi + ty_;
  }

  void multiplyPoint(const Point &point1, Point &point2) const {
    point2.x = a_*point1.x + b_*point1.y + tx_;
    point2.y = c_*point1.x + d_*point1.y + ty_;
  }

  void multiplyVector(const Vector &ivector, Vector &ovector) const {
    double ix, iy;

    ivector.getXY(&ix, &iy);

    double ox = a_*ix + b_*iy;
    double oy = c_*ix + d_*iy;

    ovector.setXY(ox, oy);
  }

  void preMultiplyPoint(double xi, double yi, double *xo, double *yo) const {
    *xo = a_*xi + c_*yi;
    *yo = b_*xi + d_*yi;
  }

  void preMultiplyPoint(const Point &ipoint, Point &opoint) const {
    opoint.x = a_*ipoint.x + c_*ipoint.y;
    opoint.y = b_*ipoint.x + d_*ipoint.y;
  }

  void preMultiplyVector(const Vector &ivector, Vector &ovector) const {
    double ix, iy;

    ivector.getXY(&ix, &iy);

    double ox = a_*ix + c_*iy;
    double oy = b_*ix + d_*iy;

    ovector.setXY(ox, oy);
  }

  const CAffineMatrix2D &translate(double x, double y) {
    tx_ += x;
    ty_ += y;

    return *this;
  }

  bool invert(CAffineMatrix2D &imatrix) const {
    double d = determinant();

    if (::fabs(d) == 0.0)
      return false;

    double id = 1.0/d;

    imatrix.a_ =  id*d_;
    imatrix.d_ =  id*a_;

    imatrix.c_ = -id*c_;
    imatrix.b_ = -id*b_;

    imatrix.tx_ =  id*calcDeterminant(b_, tx_, d_, ty_);
    imatrix.ty_ = -id*calcDeterminant(a_, tx_, c_, ty_);

    return true;
  }

  CAffineMatrix2D inverse() const {
    CAffineMatrix2D imatrix;

    if (! invert(imatrix))
      assert(false && "Divide by zero");

    return imatrix;
  }

  double determinant() const {
    return (a_*d_ - b_*c_);
  }

  const CAffineMatrix2D &normalize() {
    double d = determinant();

    double id = 1.0/d;

    a_ *= id;
    b_ *= id;
    c_ *= id;
    d_ *= id;

    return *this;
  }

  void setTransform(double xmin1, double ymin1, double xmax1, double ymax1,
                    double xmin2, double ymin2, double xmax2, double ymax2) {
    double sx = (xmax2 - xmin2)/(xmax1 - xmin1);
    double sy = (ymax2 - ymin2)/(ymax1 - ymin1);

    double tx = -xmin1*sx + xmin2;
    double ty = -ymin1*sy + ymin2;

    setInnerScale    (sx, sy);
    setOuterTranslate(tx, ty);
  }

  static CAffineMatrix2D *newIdentityMatrix() {
    CAffineMatrix2D *m = new CAffineMatrix2D();

    return m;
  }

  static bool solveAXeqB(const CAffineMatrix2D &a, CPoint2D &x, const CPoint2D &b) {
    double det_a = a.determinant();

    if (::fabs(det_a) <= 0.0)
      return false;

    double idet_a = 1.0/det_a;

    CAffineMatrix2D t(a);

    t.setColumn(0, b.x, b.y);

    double det_t = t.determinant();

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

  void zero() { memset(&a_, 0, 6*sizeof(double)); }

  const CAffineMatrix2D &operator+=(const CAffineMatrix2D &b) {
    a_ += b.a_; b_ += b.b_; tx_ += b.tx_;
    c_ += b.c_; d_ += b.d_; ty_ += b.ty_;

    return *this;
  }

  CAffineMatrix2D operator+(const CAffineMatrix2D &b) const {
    CAffineMatrix2D c = *this;

    c += b;

    return c;
  }

  const CAffineMatrix2D &operator-=(const CAffineMatrix2D &b) {
    a_ -= b.a_; b_ -= b.b_, tx_ -= b.tx_;
    c_ -= b.c_; d_ -= b.d_; ty_ -= b.ty_;

    return *this;
  }

  CAffineMatrix2D operator-(const CAffineMatrix2D &b) const {
    CAffineMatrix2D c = *this;

    c -= b;

    return c;
  }

  const CAffineMatrix2D &operator*=(const CAffineMatrix2D &b) {
    CAffineMatrix2D a = *this;

    a_  = a.a_*b.a_ + a.b_*b.c_;
    b_  = a.a_*b.b_ + a.b_*b.d_;
    c_  = a.c_*b.a_ + a.d_*b.c_;
    d_  = a.c_*b.b_ + a.d_*b.d_;

    tx_ = a.a_*b.tx_ + a.b_*b.ty_ + a.tx_;
    ty_ = a.c_*b.tx_ + a.d_*b.ty_ + a.ty_;

    return *this;
  }

  CAffineMatrix2D operator*(const CAffineMatrix2D &b) const {
    CAffineMatrix2D c = *this;

    c *= b;

    return c;
  }

  const CAffineMatrix2D &operator*=(double s) {
    CAffineMatrix2D a = *this;

    a_ = a.a_*s; b_ = a.b_*s; tx_ = a.tx_*s;
    c_ = a.c_*s; d_ = a.d_*s; ty_ = a.ty_*s;

    return *this;
  }

  CAffineMatrix2D operator*(double s) const {
    CAffineMatrix2D c = *this;

    c *= s;

    return c;
  }

  friend Point operator*(const CAffineMatrix2D &m, const Point &p) {
    Point p1;

    m.multiplyPoint(p, p1);

    return p1;
  }

  friend Point operator*(const Point &p, const CAffineMatrix2D &m) {
    Point p1;

    m.preMultiplyPoint(p, p1);

    return p1;
  }

  friend Vector operator*(const CAffineMatrix2D &m, const Vector &v) {
    Vector v1;

    m.multiplyVector(v, v1);

    return v1;
  }

  friend Vector operator*(const Vector &v, const CAffineMatrix2D &m) {
    Vector v1;

    m.preMultiplyVector(v, v1);

    return v1;
  }

  const CAffineMatrix2D &operator/=(const CAffineMatrix2D &b) {
    CAffineMatrix2D bi;

    if (! b.invert(bi)) {
      assert(false && "Divide by zero");
      return *this;
    }

    return (*this) *= bi;
  }

  CAffineMatrix2D operator/(const CAffineMatrix2D &b) const {
    CAffineMatrix2D c = *this;

    c /= b;

    return c;
  }

  //------

  void setValue(uint i, double value) {
    assert(i < 6);

    (&a_)[i] = value;
  }

  void setValue(uint i, uint j, double value) {
    assert(i < 2 && j < 2);

    double &m = (&a_)[2*j + i];

    m = value;
  }

  const double &getValue(uint i) const {
    assert(i < 6);

    return (&a_)[i];
  }

  const double &getValue(uint i, uint j) const {
    assert(i < 2 && j < 2);

    const double &m = (&a_)[2*j + i];

    return m;
  }

  const double &operator[](uint i) const {
    assert(i < 6);

    return (&a_)[i];
  }

  double &operator[](uint i) {
    assert(i < 6);

    return (&a_)[i];
  }

  //------

  void setInnerIdentity() {
    a_ = 1.0, b_ = 0.0;
    c_ = 0.0, d_ = 1.0;
  }

  void setInnerScale(double sx, double sy) {
    a_ = sx , b_ = 0.0;
    c_ = 0.0, d_ = sy ;
  }

  void setInnerRotation(double a) {
    double c = ::cos(a);
    double s = ::sin(a);

    a_ =  c, b_ = -s;
    c_ =  s, d_ =  c;
  }

  void setInnerSkew(double x, double y) {
    double tx = ::tan(x);
    double ty = ::tan(y);

    a_ = 1 , b_ = tx;
    c_ = ty, d_ = 1 ;
  }

  void setInnerSkewX(double x) {
    double tx = ::tan(x);

    a_ = 1, b_ = tx;
    c_ = 0, d_ = 1 ;
  }

  void setInnerSkewY(double y) {
    double ty = ::tan(y);

    a_ = 1 , b_ = 0;
    c_ = ty, d_ = 1;
  }

  void setInnerReflection(double dx, double dy) {
    double l = sqrt(dx*dx + dy*dy);

    setUnitInnerReflection(dx/l, dy/l);
  }

  void setUnitInnerReflection(double dx, double dy) {
  //a_ = (dx*dx - dy*dy)/l; b_ = 2*dx*dy/l;
  //c_ = b_               ; d_ = -a_;

    a_ = 2.0*dx*dx - 1.0; b_ = 2.0*dx*dy;
    c_ = b_             ; d_ = 2.0*dy*dy - 1.0;
  }

  void setOuterTranslate(double tx, double ty) {
    tx_ = tx; ty_ = ty;
  }

 private:
  static double calcDeterminant(double m00, double m01, double m10, double m11) {
    return m00*m11 - m01*m10;
  }

 private:
  double a_ { 0 }, b_ { 0 }, c_ { 0 }, d_ { 0 };
  double tx_ { 0 }, ty_ { 0 };
};

#endif
