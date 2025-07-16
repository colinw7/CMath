#ifndef CVECTOR_4D_H
#define CVECTOR_4D_H

// TODO: use CMath::sqrt() and CMath::isqrt() to allow replacement

#include <CMathGen.h>
#include <CMathUtil.h>
#include <CMathMacros.h>

class CVector4D {
 public:
  // constructor/destructor
  CVector4D() { }

  CVector4D(double x, double y, double z, double w) :
   x_(x), y_(y), z_(z), w_(w) {
  }

 ~CVector4D() { }

  //------

  // copy operations
  CVector4D(const CVector4D &vector) :
   x_(vector.x_), y_(vector.y_), z_(vector.z_), w_(vector.w_), normalized_(vector.normalized_) {
  }

  CVector4D &operator=(const CVector4D &v) {
    x_ = v.x_; y_ = v.y_; z_ = v.z_; w_ = v.w_;

    normalized_ = v.normalized_;

    return *this;
  }

  //------

  // output
  void print(std::ostream &os) const {
    os << "(" << x_ << "," << y_ << "," << z_ << "," << w_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CVector4D &vector) {
    vector.print(os);

    return os;
  }

  //------

  // accessors

  // get
  double x() const { return x_; }
  double y() const { return y_; }
  double z() const { return z_; }
  double w() const { return w_; }

  double getX() const { return x(); }
  double getY() const { return y(); }
  double getZ() const { return z(); }
  double getW() const { return w(); }

  void getXYZ(double *x, double *y, double *z, double *w) const {
    *x = x_; *y = y_; *z = z_; *w = w_;
  }

  const double *getValues() const { return &x_; }

  // Reference routine would break encapsulation

  double operator[](uint i) const { assert(i < 3); return (&x_)[i]; }

  // set
  CVector4D &setX(double x) { x_ = x; normalized_ = false; return *this; }
  CVector4D &setY(double y) { y_ = y; normalized_ = false; return *this; }
  CVector4D &setZ(double z) { z_ = z; normalized_ = false; return *this; }
  CVector4D &setW(double w) { w_ = w; normalized_ = false; return *this; }

  CVector4D &setXYZW(double x, double y, double z, double w) {
    x_ = x; y_ = y; z_ = z; w_ = w;

    normalized_ = false;

    return *this;
  }

  void iset(uint i, double v) {
    assert(i < 4);

    (&x_)[i] = v;

    normalized_ = false;
  }

  // more get accessors
  bool getNormalized() const {
    return normalized_;
  }

  //------

  // comparison
  // TODO: tolerance ? use eq()
  int cmp(const CVector4D &v) const {
    if      (x_ < v.x_) return -1;
    else if (x_ > v.x_) return  1;
    else if (y_ < v.y_) return -1;
    else if (y_ > v.y_) return  1;
    else if (z_ < v.z_) return -1;
    else if (z_ > v.z_) return  1;
    else if (w_ < v.w_) return -1;
    else if (w_ > v.w_) return  1;
    else                return  0;
  }

  friend bool operator==(const CVector4D &lhs, const CVector4D &rhs) {
    return lhs.cmp(rhs) == 0;
  }

  friend bool operator!=(const CVector4D &lhs, const CVector4D &rhs) {
    return lhs.cmp(rhs) != 0;
  }

  friend bool operator< (const CVector4D &lhs, const CVector4D &rhs) {
    return lhs.cmp(rhs) <  0;
  }

  friend bool operator<=(const CVector4D &lhs, const CVector4D &rhs) {
    return lhs.cmp(rhs) <= 0;
  }

  friend bool operator> (const CVector4D &lhs, const CVector4D &rhs) {
    return lhs.cmp(rhs) >  0;
  }

  friend bool operator>=(const CVector4D &lhs, const CVector4D &rhs) {
    return lhs.cmp(rhs) >= 0;
  }

  bool eq(const CVector4D &rhs) const {
    return fabs(x_ - rhs.x_) < 1E-6 &&
           fabs(y_ - rhs.y_) < 1E-6 &&
           fabs(z_ - rhs.z_) < 1E-6 &&
           fabs(w_ - rhs.w_) < 1E-6;
  }

  //------

#if 0
  explicit CVector4D(const CPoint4D &point) :
    x_(point.x), y_(point.y), z_(point.z), w_(point.w) {
  }
#endif

  //------

#if 0
  // to point
  CPoint4D point() const {
    return CPoint4D(x_, y_, z_, w_);
  }
#endif

  //------

  CVector4D &zero() {
    x_ = 0.0; y_ = 0.0; z_ = 0.0; w_ = 0.0;

    normalized_ = false;

    return *this;
  }

  //------

  void incX(double x = 1.0) { x_ += x; normalized_ = false; }
  void incY(double y = 1.0) { y_ += y; normalized_ = false; }
  void incZ(double z = 1.0) { z_ += z; normalized_ = false; }
  void incW(double w = 1.0) { w_ += w; normalized_ = false; }

  void decX(double x = 1.0) { x_ -= x; normalized_ = false; }
  void decY(double y = 1.0) { y_ -= y; normalized_ = false; }
  void decZ(double z = 1.0) { z_ -= z; normalized_ = false; }
  void decW(double w = 1.0) { w_ -= w; normalized_ = false; }

  void scaleX(double x = 1.0) { x_ *= x; normalized_ = false; }
  void scaleY(double y = 1.0) { y_ *= y; normalized_ = false; }
  void scaleZ(double z = 1.0) { z_ *= z; normalized_ = false; }
  void scaleW(double w = 1.0) { w_ *= w; normalized_ = false; }

  //------

  double minComponent() { return std::min(std::min(std::min(x_, y_), z_), w_); }
  double maxComponent() { return std::max(std::max(std::max(x_, y_), z_), w_); }

  //------

#if 0
  CVector4D &operator=(const CPoint4D &point) {
    x_ = point.x; y_ = point.y; z_ = point.z;

    normalized_ = false;

    return *this;
  }
#endif

  //------

  // operators

  // unary +/-
  CVector4D operator+() const {
    return CVector4D(x_, y_, z_, w_, normalized_);
  }

  CVector4D operator-() const {
    return CVector4D(-x_, -y_, -z_, -w_, normalized_);
  }

  // addition
  CVector4D &operator+=(const CVector4D &rhs) {
    x_ += rhs.x_; y_ += rhs.y_; z_ += rhs.z_; w_ += rhs.w_;

    normalized_ = false;

    return *this;
  }

  CVector4D operator+(const CVector4D &rhs) const {
    return CVector4D(x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_, w_ + rhs.w_);
  }

  // subtraction
  CVector4D &operator-=(const CVector4D &rhs) {
    x_ -= rhs.x_; y_ -= rhs.y_; z_ -= rhs.z_; w_ -= rhs.w_;

    normalized_ = false;

    return *this;
  }

  CVector4D operator-(const CVector4D &rhs) const {
    return CVector4D(x_ - rhs.x_, y_ - rhs.y_, z_ - rhs.z_, w_ - rhs.w_);
  }

  // scalar multiplication/division
  CVector4D &operator*=(const CVector4D &rhs) {
    x_ *= rhs.x_; y_ *= rhs.y_; z_ *= rhs.z_; w_ *= rhs.w_;

    normalized_ = false;

    return *this;
  }

  CVector4D &operator*=(double rhs) {
    x_ *= rhs; y_ *= rhs; z_ *= rhs; w_ *= rhs;

    normalized_ = false;

    return *this;
  }

  CVector4D operator*(const CVector4D &rhs) {
    CVector4D t(*this);

    t *= rhs;

    return t;
  }

  friend CVector4D operator*(const CVector4D &lhs, double rhs) {
    return CVector4D(lhs.x_*rhs, lhs.y_*rhs, lhs.z_*rhs, lhs.w_*rhs);
  }

  friend CVector4D operator*(double lhs, const CVector4D &rhs) {
    return CVector4D(lhs*rhs.x_, lhs*rhs.y_, lhs*rhs.z_, lhs*rhs.w_);
  }

  CVector4D &operator/=(double rhs) {
    double irhs = 1.0/rhs;

    x_ *= irhs; y_ *= irhs; z_ *= irhs; w_ *= irhs;

    normalized_ = false;

    return *this;
  }

  CVector4D operator/(double rhs) const {
    double irhs = 1.0/rhs;

    return CVector4D(x_*irhs, y_*irhs, z_*irhs, w_*irhs);
  }

  //------

#if 0
  friend CPoint4D operator+(const CPoint4D &lhs, const CVector4D &rhs) {
    return CPoint4D(lhs.x + rhs.x_, lhs.y + rhs.y_, lhs.z + rhs.z_, lhs.w + rhs.w_);
  }

  friend CPoint4D operator+=(CPoint4D &lhs, const CVector4D &rhs) {
    lhs.x += rhs.x_; lhs.y += rhs.y_; lhs.z += rhs.z_; lhs.w += rhs.w_;

    return lhs;
  }

  friend CPoint4D operator-(const CPoint4D &lhs, const CVector4D &rhs) {
    return CPoint4D(lhs.x - rhs.x_, lhs.y - rhs.y_, lhs.z - rhs.z_, lhs.w - rhs.w_);
  }

  friend CPoint4D operator-=(CPoint4D &lhs, const CVector4D &rhs) {
    lhs.x -= rhs.x_; lhs.y -= rhs.y_; lhs.z -= rhs.z_; lhs.w -= rhs.w_;

    return lhs;
  }
#endif

#if 0
  friend CVector4D operator-(const CPoint4D &lhs, const CPoint4D &rhs) {
    return CVector4D(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z, lhs.w - rhs.w);
  }
#endif

  //------

 private:
  CVector4D(double x, double y, double z, double w, bool normalized) :
   x_(x), y_(y), z_(z), w_(w), normalized_(normalized) {
  }

 private:
  double x_ { 0 }, y_ { 0 }, z_ { 0 }, w_ { 1.0 };
  bool   normalized_ { false };
};

#endif
