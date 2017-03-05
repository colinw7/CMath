#ifndef CIVECTOR_2D_H
#define CIVECTOR_2D_H

#include <iostream>
#include <CMathGen.h>

class CIVector2D {
 public:
  explicit CIVector2D(int x=0, int y=0) :
   x_(x), y_(y) {
  }

  CIVector2D(const CIVector2D &vector) :
   x_(vector.x_), y_(vector.y_) {
  }

  explicit CIVector2D(const CIPoint2D &point) :
   x_(point.x), y_(point.y) {
  }

  CIVector2D(const CIPoint2D &point1, const CIPoint2D &point2) :
   x_(point2.x - point1.x), y_(point2.y - point1.y) {
  }

  CIVector2D(const CIVector2D &vector1, const CIVector2D &vector2) :
   x_(vector2.x_ - vector1.x_), y_(vector2.y_ - vector1.y_) {
  }

  int getX() const { return x_; }
  int getY() const { return y_; }

  void getXY(int *x, int *y) const {
    *x = x_; *y = y_;
  }

  void setX(int x) {
    x_ = x;
  }

  void setY(int y) {
    y_ = y;
  }

  void setXY(int x, int y) {
    x_ = x; y_ = y;
  }

  CIPoint2D point() const {
    return CIPoint2D(x_, y_);
  }

  double length() const {
    return sqrt(x_*x_ + y_*y_);
  }

  int fastLength() const {
    return CMathGen::fastDistance(x_, y_);
  }

  int lengthSqr() const {
    return (x_*x_ + y_*y_);
  }

  double modulus() const {
    return length();
  }

  CIVector2D &zero() {
    x_ = 0.0; y_ = 0.0;

    return *this;
  }

  void incX(int x) {
    x_ += x;
  }

  void incY(int y) {
    y_ += y;
  }

  int minComponent() {
    return std::min(x_, y_);
  }

  int maxComponent() {
    return std::max(x_, y_);
  }

  int minAbsComponent() {
    return std::min(std::abs(x_), std::abs(y_));
  }

  int maxAbsComponent() {
    return std::max(std::abs(x_), std::abs(y_));
  }

  CIVector2D operator+() const {
    return CIVector2D(x_, y_);
  }

  CIVector2D operator-() const {
    return CIVector2D(-x_, -y_);
  }

  friend bool operator==(const CIVector2D &lhs, const CIVector2D &rhs) {
    return (lhs.x_ == rhs.x_ && lhs.y_ == rhs.y_);
  }

  friend bool operator!=(const CIVector2D &lhs, const CIVector2D &rhs) {
    return (lhs.x_ != rhs.x_ || lhs.y_ != rhs.y_);
  }

  CIVector2D &operator+=(const CIVector2D &rhs) {
    x_ += rhs.x_; y_ += rhs.y_;

    return *this;
  }

  CIVector2D operator+(const CIVector2D &rhs) const {
    return CIVector2D(x_ + rhs.x_, y_ + rhs.y_);
  }

  CIVector2D &operator-=(const CIVector2D &rhs) {
    x_ -= rhs.x_; y_ -= rhs.y_;

    return *this;
  }

  CIVector2D operator-(const CIVector2D &rhs) const {
    return CIVector2D(x_ - rhs.x_, y_ - rhs.y_);
  }

  //---

  CIVector2D &operator*=(int rhs) {
    x_ *= rhs; y_ *= rhs;

    return *this;
  }

  friend CIVector2D operator*(const CIVector2D &lhs, int rhs) {
    return CIVector2D(lhs.x_*rhs, lhs.y_*rhs);
  }

  friend CIVector2D operator*(int lhs, CIVector2D &rhs) {
    return CIVector2D(lhs*rhs.x_, lhs*rhs.y_);
  }

  //---

  CIVector2D &operator/=(int rhs) {
    x_ /= rhs; y_ /= rhs;

    return *this;
  }

  friend CIVector2D operator/(const CIVector2D &lhs, int rhs) {
    return CIVector2D(lhs.x_/rhs, lhs.y_/rhs);
  }

  //---

  void print(std::ostream &os) const {
    os << "(" << x_ << "," << y_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CIVector2D &vector) {
    vector.print(os);

    return os;
  }

  CIVector2D crossProduct(const CIVector2D &vector) const {
    return CIVector2D(y_*vector.x_ - x_*vector.y_,
                       x_*vector.y_ - y_*vector.x_);
  }

  CIVector2D crossProduct(int x, int y) const {
    return CIVector2D(y_*x - x_*y, x_*y - y_*x);
  }

  static CIVector2D crossProduct(const CIVector2D &vector1,
                                const CIVector2D &vector2) {
    return CIVector2D(vector1.y_*vector2.x_ - vector1.x_*vector2.y_,
                       vector1.x_*vector2.y_ - vector1.y_*vector2.x_);
  }

  int dotProduct(const CIVector2D &vector) const {
    return (x_*vector.x_ + y_*vector.y_);
  }

  int dotProduct(int x, int y) const {
    return (x_*x + y_*y);
  }

  static int dotProduct(const CIVector2D &vector1,
                      const CIVector2D &vector2) {
    return (vector1.x_*vector2.x_ + vector1.y_*vector2.y_);
  }

  static int dotProduct(const CIVector2D &vector1,
                      int x2, int y2) {
    return (vector1.x_*x2 + vector1.y_*y2);
  }

  int dotProductSelf() const {
    return (x_*x_ + y_*y_);
  }

 private:
  int x_ { 0 }, y_ { 0 };
};

#endif
