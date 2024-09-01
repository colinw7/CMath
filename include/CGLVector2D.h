#ifndef CGLVector2D_H
#define CGLVector2D_H

// TODO: use CMath::sqrt() and CMath::isqrt() to allow replacement

#include <CMathGen.h>
#include <CMathUtil.h>
#include <CPoint2D.h>
#include <CMathMacros.h>

class CGLVector2D {
 public:
  enum Type {
    ZERO,
    UNIT,
    UNIT_X,
    UNIT_Y
  };

 public:
  // constructor/destructor
  CGLVector2D() { }

  CGLVector2D(float x, float y) :
   x_(x), y_(y) {
  }

 ~CGLVector2D() { }

  //------

  // copy operations
  CGLVector2D(const CGLVector2D &vector) :
   x_(vector.x_), y_(vector.y_) {
  }

  CGLVector2D &operator=(const CGLVector2D &v) {
    x_ = v.x_; y_ = v.y_;

    return *this;
  }

  //------

  // output
  void print(std::ostream &os) const {
    os << "(" << x_ << "," << y_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CGLVector2D &vector) {
    vector.print(os);

    return os;
  }

  //------

  // accessors

  // get
  float x() const { return x_; }
  float y() const { return y_; }

  float getX() const { return x_; }
  float getY() const { return y_; }

  void getXY(float *x, float *y) const {
    *x = x_; *y = y_;
  }

  const float *getValues() const { return &x_; }

  // Reference routine would break encapsulation

  float operator[](uint i) const { assert(i < 3); return (&x_)[i]; }

  // set
  CGLVector2D &setX(float x) {
    x_ = x;

    return *this;
  }

  CGLVector2D &setY(float y) {
    y_ = y;

    return *this;
  }

  CGLVector2D &setXY(float x, float y) {
    x_ = x; y_ = y;

    return *this;
  }

  void iset(uint i, float v) {
    assert(i < 3);

    (&x_)[i] = v;
  }

  // more get accessors
  bool getNormalized() const {
    return fabs(length() - 1.0) < 1E-6;
  }

  float length() const {
    return float(std::sqrt(lengthSqr()));
  }

  float fastLength() const {
    return float(CMathGen::fastDistance(double(x_), double(y_)));
  }

  float modulus() const {
    return length();
  }

  float lengthSqr() const {
    return (x_*x_ + y_*y_);
  }

  //------

  // comparison
  // TODO: tolerance ? use eq()
  int cmp(const CGLVector2D &v) const {
    if      (x_ < v.x_) return -1;
    else if (x_ > v.x_) return  1;
    else if (y_ < v.y_) return -1;
    else if (y_ > v.y_) return  1;
    else                return  0;
  }

  friend bool operator==(const CGLVector2D &lhs, const CGLVector2D &rhs) {
    return lhs.cmp(rhs) == 0;
  }

  friend bool operator!=(const CGLVector2D &lhs, const CGLVector2D &rhs) {
    return lhs.cmp(rhs) != 0;
  }

  friend bool operator< (const CGLVector2D &lhs, const CGLVector2D &rhs) {
    return lhs.cmp(rhs) <  0;
  }

  friend bool operator<=(const CGLVector2D &lhs, const CGLVector2D &rhs) {
    return lhs.cmp(rhs) <= 0;
  }

  friend bool operator> (const CGLVector2D &lhs, const CGLVector2D &rhs) {
    return lhs.cmp(rhs) >  0;
  }

  friend bool operator>=(const CGLVector2D &lhs, const CGLVector2D &rhs) {
    return lhs.cmp(rhs) >= 0;
  }

  bool eq(const CGLVector2D &rhs) const {
    return fabs(x_ - rhs.x_) < 1E-6 && fabs(y_ - rhs.y_) < 1E-6;
  }

  //------

  explicit CGLVector2D(const CPoint2D &point) :
    x_(float(point.x)), y_(float(point.y)) {
  }

  // v = point2 - point1
  CGLVector2D(const CPoint2D &point1, const CPoint2D &point2) :
   x_(float(point2.x - point1.x)), y_(float(point2.y - point1.y)) {
  }

  // v = vector2 - vector1
  CGLVector2D(const CGLVector2D &vector1, const CGLVector2D &vector2) :
   x_(float(vector2.x_ - vector1.x_)), y_(float(vector2.y_ - vector1.y_)) {
  }

  //------

  // parametric creation
  CGLVector2D(Type type, float s=1.0) {
    if      (type == ZERO  ) { x_ = 0.0; y_ = 0.0; }
    else if (type == UNIT  ) { x_ =   s; y_ =   s; }
    else if (type == UNIT_X) { x_ =   s; y_ = 0.0; }
    else if (type == UNIT_Y) { x_ = 0.0; y_ =   s; }
  }

  //------

  // to point
  CPoint2D point() const {
    return CPoint2D(x_, y_);
  }

  //------

  CGLVector2D &zero() {
    x_ = 0.0; y_ = 0.0;

    return *this;
  }

  bool isZero() const {
    return CMathUtil::realEq(lengthSqr(), 0);
  }

  bool isUnit() const {
    return CMathUtil::realEq(lengthSqr(), 1);
  }

  //------

  CGLVector2D &normalize() {
    float len = length();

    // assert on len == 0 ?

    float factor = 0.0;

    if (len > 0.0)
      factor = 1.0f/len;

    x_ *= factor;
    y_ *= factor;

    return *this;
  }

  CGLVector2D normalized() const {
    float len = length();

    // assert on len == 0 ?

    float factor = 0.0;

    if (len > 0.0)
      factor = 1.0f/len;

    return CGLVector2D(x_*factor, y_*factor);
  }

  CGLVector2D unit() const {
    return normalized();
  }

  //------

  CGLVector2D &setMagnitude(float magnitude) {
    float factor = 0.0;

    float len = length();

    if (len > 0.0)
      factor = magnitude/len;

    x_ *= factor; y_ *= factor;

    return *this;
  }

  //------

  float getDistance(const CGLVector2D &vector) const {
    CGLVector2D diff = *this - vector;

    return diff.length();
  }

  float getDistanceSqr(const CGLVector2D &vector) const {
    CGLVector2D diff = *this - vector;

    return diff.lengthSqr();
  }

  //------

  void incX(float x = 1.0) { x_ += x; }
  void incY(float y = 1.0) { y_ += y; }

  void decX(float x = 1.0) { x_ -= x; }
  void decY(float y = 1.0) { y_ -= y; }

  void scaleX(float x = 1.0) { x_ *= x; }
  void scaleY(float y = 1.0) { y_ *= y; }

  //------

  float minComponent() { return std::min(x_, y_); }
  float maxComponent() { return std::max(x_, y_); }

  float minAbsComponent() { return float(std::min(::fabs(x_), ::fabs(y_))); }
  float maxAbsComponent() { return float(std::max(::fabs(x_), ::fabs(y_))); }

  //------

  static CGLVector2D min(const CGLVector2D &lhs, const CGLVector2D &rhs) {
    return CGLVector2D(std::min(lhs.x_, rhs.x_), std::min(lhs.y_, rhs.y_));
  }

  static CGLVector2D max(const CGLVector2D &lhs, const CGLVector2D &rhs) {
    return CGLVector2D(std::max(lhs.x_, rhs.x_), std::max(lhs.y_, rhs.y_));
  }

  //------

  CGLVector2D &operator=(const CPoint2D &point) {
    x_ = float(point.x); y_ = float(point.y);

    return *this;
  }

  //------

  // operators

  // unary +/-
  CGLVector2D operator+() const {
    return CGLVector2D(x_, y_);
  }

  CGLVector2D operator-() const {
    return CGLVector2D(-x_, -y_);
  }

  // addition
  CGLVector2D &operator+=(const CGLVector2D &rhs) {
    x_ += rhs.x_; y_ += rhs.y_;

    return *this;
  }

  CGLVector2D operator+(const CGLVector2D &rhs) const {
    return CGLVector2D(x_ + rhs.x_, y_ + rhs.y_);
  }

  // subtraction
  CGLVector2D &operator-=(const CGLVector2D &rhs) {
    x_ -= rhs.x_; y_ -= rhs.y_;

    return *this;
  }

  CGLVector2D operator-(const CGLVector2D &rhs) const {
    return CGLVector2D(x_ - rhs.x_, y_ - rhs.y_);
  }

  // scalar multiplication/division
  CGLVector2D &operator*=(const CGLVector2D &rhs) {
    x_ *= rhs.x_; y_ *= rhs.y_;

    return *this;
  }

  CGLVector2D &operator*=(float rhs) {
    x_ *= rhs; y_ *= rhs;

    return *this;
  }

  CGLVector2D operator*(const CGLVector2D &rhs) {
    CGLVector2D t(*this);

    t *= rhs;

    return t;
  }

  friend CGLVector2D operator*(const CGLVector2D &lhs, float rhs) {
    return CGLVector2D(lhs.x_*rhs, lhs.y_*rhs);
  }

  friend CGLVector2D operator*(float lhs, const CGLVector2D &rhs) {
    return CGLVector2D(lhs*rhs.x_, lhs*rhs.y_);
  }

  CGLVector2D &operator/=(float rhs) {
    float irhs = 1.0f/rhs;

    x_ *= irhs; y_ *= irhs;

    return *this;
  }

  CGLVector2D operator/(float rhs) const {
    float irhs = 1.0f/rhs;

    return CGLVector2D(x_*irhs, y_*irhs);
  }

  //------

  // dot product
  float dotProduct(const CGLVector2D &v) const {
    return (x_*v.x_ + y_*v.y_);
  }

//float dotProduct(const Normal &normal) const {
//  return (x_*normal.getX() + y_*normal.getY());
//}

  static float dotProduct(const CGLVector2D &v1, const CGLVector2D &v2) {
    return v1.dotProduct(v2);
  }

  static float dotProduct(const CPoint2D &v1, const CGLVector2D &v2) {
    return CGLVector2D(float(v1.x), float(v1.y)).dotProduct(CGLVector2D(v2));
  }

  float dotProduct(const CPoint2D &point) const {
    return dotProduct(CGLVector2D(float(point.x), float(point.y)));
  }

  float dotProduct(float x, float y) const {
    return dotProduct(CGLVector2D(x, y));
  }

  static float dotProduct(const CGLVector2D &v1, float x2, float y2) {
    return v1.dotProduct(CGLVector2D(x2, y2));
  }

  float dotProductSelf() const {
    return dotProduct(*this);
  }

  static float absDotProduct(const CGLVector2D &v1, const CGLVector2D &v2) {
    return float(std::fabs(dotProduct(v1, v2)));
  }

  //------

  // cross product
  CGLVector2D crossProduct(const CGLVector2D &v) const {
    return CGLVector2D(y_*v.x_ - x_*v.y_, x_*v.y_ - y_*v.x_);
  }

  static CGLVector2D crossProduct(const CGLVector2D &v1, const CGLVector2D &v2) {
    return v1.crossProduct(v2);
  }

  CGLVector2D crossProduct(const CPoint2D &point) const {
    return crossProduct(CGLVector2D(float(point.x), float(point.y)));
  }

  CGLVector2D crossProduct(float x, float y) const {
    return crossProduct(CGLVector2D(x, y));
  }

  CGLVector2D unitCrossProduct(const CGLVector2D &v) const {
    return crossProduct(v).normalize();
  }

  // Note: area of parallelogram is
  // CVector v1(x2 - x1, y2 - y1), v2(x4 - x1, y4 - y1);
  // area = v1.crossProduct(v2).length();

  //------

  friend CPoint2D operator+(const CPoint2D &lhs, const CGLVector2D &rhs) {
    return CPoint2D(float(lhs.x + rhs.x_), float(lhs.y + rhs.y_));
  }

  friend CPoint2D operator+=(CPoint2D &lhs, const CGLVector2D &rhs) {
    lhs.x += rhs.x_; lhs.y += rhs.y_;

    return lhs;
  }

  friend CPoint2D operator-(const CPoint2D &lhs, const CGLVector2D &rhs) {
    return CPoint2D(lhs.x - rhs.x_, lhs.y - rhs.y_);
  }

  friend CPoint2D operator-=(CPoint2D &lhs, const CGLVector2D &rhs) {
    lhs.x -= rhs.x_; lhs.y -= rhs.y_;

    return lhs;
  }

#if 0
  friend CGLVector2D operator-(const CPoint2D &lhs, const CPoint2D &rhs) {
    return CGLVector2D(lhs.x - rhs.x, lhs.y - rhs.y);
  }
#endif

  //------

  CGLVector2D normal(const CGLVector2D &vector2) const {
    return unitCrossProduct(vector2);
  }

  static CGLVector2D normal(const CGLVector2D &vector1, const CGLVector2D &vector2) {
    return vector1.normal(vector2);
  }

  //------

  float cosIncluded(const CGLVector2D &vector1) const {
    float dot = dotProduct(vector1);

    dot /= length();
    dot /= vector1.length();

    return dot;
  }

  static float cosIncluded(const CGLVector2D &vector1, const CGLVector2D &vector2) {
    return vector1.cosIncluded(vector2);
  }

  //------

 private:
  float x_ { 0.0f };
  float y_ { 0.0f };
};

#endif
