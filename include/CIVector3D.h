#ifndef CIVECTOR_3D_H
#define CIVECTOR_3D_H

#include <CMathGen.h>

class CIVector3D {
 private:
  int x_, y_, z_;

 public:
  CIVector3D(int x=0, int y=0, int z=0) :
   x_(x), y_(y), z_(z) {
  }

  CIVector3D(const CIVector3D &vector) :
    x_(vector.x_), y_(vector.y_), z_(vector.z_) {
  }

  CIVector3D(const CIVector3D &vector1, const CIVector3D &vector2) :
    x_(vector2.x_ - vector1.x_), y_(vector2.y_ - vector1.y_),
    z_(vector2.z_ - vector1.z_) {
  }

  int getX() const { return x_; }
  int getY() const { return y_; }
  int getZ() const { return z_; }

  void getXYZ(int *x, int *y, int *z) const {
    *x = x_; *y = y_; *z = z_;
  }

  void setX(int x) {
    x_ = x;
  }

  void setY(int y) {
    y_ = y;
  }

  void setZ(int z) {
    z_ = z;
  }

  void setXYZ(int x, int y, int z) {
    x_ = x; y_ = y; z_ = z;
  }

  double length() const {
    return sqrt(x_*x_ + y_*y_ + z_*z_);
  }

  double fastLength() const {
    return CMathGen::fastDistance(x_, y_, z_);
  }

  int modulus() const {
    return length();
  }

  int lengthSqr() const {
    return (x_*x_ + y_*y_ + z_*z_);
  }

  CIVector3D &zero() {
    x_ = 0; y_ = 0; z_ = 0;

    return *this;
  }

  int getDistance(const CIVector3D &vector) {
    CIVector3D diff = *this - vector;

    return diff.length();
  }

  void incX(int x = 1) {
    x_ += x;
  }

  void incY(int y = 1) {
    y_ += y;
  }

  void incZ(int z = 1) {
    z_ += z;
  }

  void decX(int x = 1) {
    x_ -= x;
  }

  void decY(int y = 1) {
    y_ -= y;
  }

  void decZ(int z = 1) {
    z_ -= z;
  }

  int minComponent() {
    return min(min(x_, y_), z_);
  }

  int maxComponent() {
    return max(max(x_, y_), z_);
  }

  int minAbsComponent() {
    return min(min(fabs(x_), fabs(y_)), fabs(z_));
  }

  int maxAbsComponent() {
    return max(max(fabs(x_), fabs(y_)), fabs(z_));
  }

  CIVector3D operator+() const {
    return CIVector3D(x_, y_, z_);
  }

  CIVector3D operator-() const {
    return CIVector3D(-x_, -y_, -z_);
  }

  friend bool operator==(const CIVector3D &lhs, const CIVector3D &rhs) {
    return (lhs.x_ == rhs.x_ && lhs.y_ == rhs.y_ && lhs.z_ == rhs.z_);
  }

  friend bool operator!=(const CIVector3D &lhs, const CIVector3D &rhs) {
    return (lhs.x_ != rhs.x_ || lhs.y_ != rhs.y_ || lhs.z_ != rhs.z_);
  }

  CIVector3D &operator+=(const CIVector3D &rhs) {
    x_ += rhs.x_; y_ += rhs.y_; z_ += rhs.z_;

    return *this;
  }

  CIVector3D operator+(const CIVector3D &rhs) const {
    return CIVector3D(x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_);
  }

  CIVector3D &operator-=(const CIVector3D &rhs) {
    x_ -= rhs.x_; y_ -= rhs.y_; z_ -= rhs.z_;

    return *this;
  }

  CIVector3D operator-(const CIVector3D &rhs) const {
    return CIVector3D(x_ - rhs.x_, y_ - rhs.y_, z_ - rhs.z_);
  }

  CIVector3D &operator*=(int rhs) {
    x_ *= rhs; y_ *= rhs; z_ *= rhs;

    return *this;
  }

  friend CIVector3D operator*(const CIVector3D &lhs, int rhs) {
    return CIVector3D(lhs.x_*rhs, lhs.y_*rhs, lhs.z_*rhs);
  }

  friend CIVector3D operator*(int lhs, const CIVector3D &rhs) {
    return CIVector3D(lhs*rhs.x_, lhs*rhs.y_, lhs*rhs.z_);
  }

  void print(ostream &os) const {
    os << "(" << x_ << "," << y_ << "," << z_ << ")";
  }

  friend ostream &operator<<(ostream &os, const CIVector3D &vector) {
    vector.print(os);

    return os;
  }

  CIVector3D crossProduct(const CIVector3D &vector) const {
    return CIVector3D(y_*vector.z_ - z_*vector.y_,
                      z_*vector.x_ - x_*vector.z_,
                      x_*vector.y_ - y_*vector.x_);
  }

  CIVector3D crossProduct(int x, int y, int z) const {
    return CIVector3D(y_*z - z_*y, z_*x - x_*z, x_*y - y_*x);
  }

  static CIVector3D crossProduct(const CIVector3D &vector1,
                                const CIVector3D &vector2) {
    return CIVector3D(vector1.y_*vector2.z_ - vector1.z_*vector2.y_,
                      vector1.z_*vector2.x_ - vector1.x_*vector2.z_,
                      vector1.x_*vector2.y_ - vector1.y_*vector2.x_);
  }

  int dotProduct(const CIVector3D &vector) const {
    return (x_*vector.x_ + y_*vector.y_ + z_*vector.z_);
  }

  int dotProduct(int x, int y, int z) const {
    return (x_*x + y_*y + z_*z);
  }

  static int dotProduct(const CIVector3D &vector1,
                        const CIVector3D &vector2) {
    return (vector1.x_*vector2.x_ +
            vector1.y_*vector2.y_ +
            vector1.z_*vector2.z_);
  }

  static int dotProduct(const CIVector3D &vector1,
                        int x2, int y2, int z2) {
    return (vector1.x_*x2 + vector1.y_*y2 + vector1.z_*z2);
  }

  int dotProductSelf() const {
    return (x_*x_ + y_*y_ + z_*z_);
  }

  CIVector3D normal(const CIVector3D &vector2) {
    return crossProduct(vector2).normalize();
  }

  static CIVector3D normal(const CIVector3D &vector1,
                           const CIVector3D &vector2) {
    return crossProduct(vector1, vector2).normalize();
  }
};

#endif
