#ifndef CIVECTOR_2D_H
#define CIVECTOR_2D_H

#include <iostream>

template<typename T>
class CIVector2DT {
 private:
  typedef CIPoint2DT<T> Point;

  T x_, y_;

 public:
  CIVector2DT(T x=0, T y=0) :
   x_(x), y_(y) {
  }

  CIVector2DT(const CIVector2DT &vector) :
   x_(vector.x_), y_(vector.y_) {
  }

  explicit CIVector2DT(const Point &point) :
   x_(point.x), y_(point.y) {
  }

  CIVector2DT(const Point &point1, const Point &point2) :
   x_(point2.x - point1.x), y_(point2.y - point1.y) {
  }

  CIVector2DT(const CIVector2DT &vector1, const CIVector2DT &vector2) :
   x_(vector2.x_ - vector1.x_), y_(vector2.y_ - vector1.y_) {
  }

  T getX() const { return x_; }
  T getY() const { return y_; }

  void getXY(T *x, T *y) const {
    *x = x_; *y = y_;
  }

  void setX(T x) {
    x_ = x;
  }

  void setY(T y) {
    y_ = y;
  }

  void setXY(T x, T y) {
    x_ = x; y_ = y;
  }

  Point point() const {
    return Point(x_, y_);
  }

  double length() const {
    return sqrt(x_*x_ + y_*y_);
  }

  T lengthSqr() const {
    return (x_*x_ + y_*y_);
  }

  double modulus() const {
    return length();
  }

  CIVector2DT &zero() {
    x_ = 0.0; y_ = 0.0;

    return *this;
  }

  void incX(T x) {
    x_ += x;
  }

  void incY(T y) {
    y_ += y;
  }

  T minComponent() {
    return min(x_, y_);
  }

  T maxComponent() {
    return max(x_, y_);
  }

  T minAbsComponent() {
    return min(abs(x_), abs(y_));
  }

  T maxAbsComponent() {
    return max(abs(x_), abs(y_));
  }

  CIVector2DT operator+() const {
    return CIVector2DT(x_, y_);
  }

  CIVector2DT operator-() const {
    return CIVector2DT(-x_, -y_);
  }

  friend bool operator==(const CIVector2DT &lhs, const CIVector2DT &rhs) {
    return (lhs.x_ == rhs.x_ && lhs.y_ == rhs.y_);
  }

  friend bool operator!=(const CIVector2DT &lhs, const CIVector2DT &rhs) {
    return (lhs.x_ != rhs.x_ || lhs.y_ != rhs.y_);
  }

  CIVector2DT &operator+=(const CIVector2DT &rhs) {
    x_ += rhs.x_; y_ += rhs.y_;

    return *this;
  }

  CIVector2DT operator+(const CIVector2DT &rhs) const {
    return CIVector2DT(x_ + rhs.x_, y_ + rhs.y_);
  }

  CIVector2DT &operator-=(const CIVector2DT &rhs) {
    x_ -= rhs.x_; y_ -= rhs.y_;

    return *this;
  }

  CIVector2DT operator-(const CIVector2DT &rhs) const {
    return CIVector2DT(x_ - rhs.x_, y_ - rhs.y_);
  }

  CIVector2DT &operator*=(T rhs) {
    x_ *= rhs; y_ *= rhs;

    return *this;
  }

  friend CIVector2DT operator*(const CIVector2DT &lhs, T rhs) {
    return CIVector2DT(lhs.x_*rhs, lhs.y_*rhs);
  }

  friend CIVector2DT operator*(T lhs, CIVector2DT &rhs) {
    return CIVector2DT(lhs*rhs.x_, lhs*rhs.y_);
  }

  void print(std::ostream &os) const {
    os << "(" << x_ << "," << y_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CIVector2DT &vector) {
    vector.print(os);

    return os;
  }

  CIVector2DT crossProduct(const CIVector2DT &vector) const {
    return CIVector2DT(y_*vector.x_ - x_*vector.y_,
                       x_*vector.y_ - y_*vector.x_);
  }

  CIVector2DT crossProduct(T x, T y) const {
    return CIVector2DT(y_*x - x_*y, x_*y - y_*x);
  }

  static CIVector2DT crossProduct(const CIVector2DT &vector1,
                                const CIVector2DT &vector2) {
    return CIVector2DT(vector1.y_*vector2.x_ - vector1.x_*vector2.y_,
                       vector1.x_*vector2.y_ - vector1.y_*vector2.x_);
  }

  T dotProduct(const CIVector2DT &vector) const {
    return (x_*vector.x_ + y_*vector.y_);
  }

  T dotProduct(T x, T y) const {
    return (x_*x + y_*y);
  }

  static T dotProduct(const CIVector2DT &vector1,
                      const CIVector2DT &vector2) {
    return (vector1.x_*vector2.x_ + vector1.y_*vector2.y_);
  }

  static T dotProduct(const CIVector2DT &vector1,
                      T x2, T y2) {
    return (vector1.x_*x2 + vector1.y_*y2);
  }

  T dotProductSelf() const {
    return (x_*x_ + y_*y_);
  }
};

typedef class CIVector2DT<int> CIVector2D;

#endif
