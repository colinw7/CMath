#ifndef CNORM_VECTOR_3D_H
#define CNORM_VECTOR_3D_H

#include <CPoint3D.h>
#include <CMathGen.h>

template<typename T>
class CNormVector3DT {
 private:
  T    x_, y_, z_, w_;
  bool normalized_;

 public:
  CNormVector3DT(T x=0.0, T y=0.0, T z=0.0, T w=1.0) :
   x_(x), y_(y), z_(z), w_(w), normalized_(false) {
  }

  CNormVector3DT(const CNormVector3DT &vector) :
    x_(vector.x_), y_(vector.y_), z_(vector.z_), w_(vector.w_)
    normalized_(vector.normalized_) {
  }

  explicit CNormVector3DT(const CVector3DT<T> &vector) :
    x_(vector.getX()), y_(vector.getY()), z_(vector.getZ()), w_(1.0)
    normalized_(vector.getNormalized()) {
  }

  explicit CNormVector3DT(const CPoint3DT<T> &point) :
    x_(point.x_), y_(point.y_), z_(point.z_), w_(1.0) normalized_(false) {
  }

  CNormVector3DT(const CNormVector3DT &vector1, const CNormVector3DT &vector2) :
    x_(vector2.x_ - vector1.x_), y_(vector2.y_ - vector1.y_),
    z_(vector2.z_ - vector1.z_), w_(1.0), normalized_(false) {
  }

  CNormVector3DT(const CVector3DT<T> &vector1, const CVector3DT<T> &vector2) :
    x_(vector2.x_ - vector1.x_), y_(vector2.y_ - vector1.y_),
    z_(vector2.z_ - vector1.z_), w_(1.0), normalized_(false) {
  }

  CNormVector3DT(const CPoint3DT<T> &point1, const CPoint3DT<T> &point2) :
    x_(point2.x_ - point1.x_), y_(point2.y_ - point1.y_),
    z_(point2.z_ - point1.z_), w_(1.0), normalized_(false) {
  }

  T getX() const { return x_; }
  T getY() const { return y_; }
  T getZ() const { return z_; }
  T getW() const { return w_; }

  void getXYZW(T *x, T *y, T *z, T *w) const {
    *x = x_; *y = y_; *z = z_; *w = w_;
  }

  void getXYZ(T *x, T *y, T *z) const {
    *x = x_; *y = y_; *z = z_;
  }

  void setX(T x) {
    x_ = x;

    normalized_ = false;
  }

  void setY(T y) {
    y_ = y;

    normalized_ = false;
  }

  void setZ(T z) {
    z_ = z;

    normalized_ = false;
  }

  void setW(T w) {
    w_ = w;

    normalized_ = false;
  }

  void setXYZW(T x, T y, T z, T w) {
    x_ = x; y_ = y; z_ = z; w_ = w;

    normalized_ = false;
  }

  void setXYZ(T x, T y, T z) {
    x_ = x; y_ = y; z_ = z;

    normalized_ = false;
  }

  CPoint3DT<T> point() const {
    return CPoint3DT<T>(x_, y_, z_);
  }

  CVector3DT<T> vector() const {
    return CVector3DT<T>(x_, y_, z_);
  }

  T length() const {
    if (normalized_)
      return 1.0;
    else
      return ::sqrt(x_*x_ + y_*y_ + z_*z_);
  }

  T fastLength() const {
    if (normalized_)
      return 1.0;
    else
      return CMathGen::fastLength(x_, y_, z_);
  }

  T modulus() const {
    return length();
  }

  T lengthSqr() const {
    if (normalized_)
      return 1.0;
    else
      return (x_*x_ + y_*y_ + z_*z_);
  }

  CNormVector3DT &zero() {
    x_ = 0.0; y_ = 0.0; z_ = 0.0; w_ = 1.0;

    normalized_ = false;

    return *this;
  }

  CNormVector3DT unit() const {
    if (normalized_)
      return *this;

    T len = length();

    T factor = 0.0;

    if (len > 0.0)
      factor = 1.0/len;

    return CNormVector3DT(x_*factor, y_*factor, z_*factor, 1.0);
  }

  CNormVector3DT &normalize() {
    if (normalized_)
      return *this;

    T len = length();

    T factor = 0.0;

    if (len > 0.0)
      factor = 1.0/len;

    x_ *= factor;
    y_ *= factor;
    z_ *= factor;
    w_  = 1.0;

    normalized_ = true;

    return *this;
  }

  CNormVector3DT normalized() const {
    if (normalized_)
      return *this;

    T len = length();

    T factor = 0.0;

    if (len > 0.0)
      factor = 1.0/len;

    return CNormVector3DT(x_*factor, y_*factor, z_*factor, 1.0);
  }

  CNormVector3DT &setMagnitude(T magnitude) {
    T factor = 0.0;

    if (normalized_)
      factor = magnitude;
    else {
      T len = length();

      if (len > 0.0)
        factor = magnitude/len;
    }

    x_ *= factor;
    y_ *= factor;
    z_ *= factor;
    w_  = 1.0;

    normalized_ = false;

    return *this;
  }

  T getDistance(const CNormVector3DT &vector) {
    CNormVector3DT diff = *this - vector;

    return diff.length();
  }

  void incX(T x) {
    x_ += x;

    normalized_ = false;
  }

  void incY(T y) {
    y_ += y;

    normalized_ = false;
  }

  void incZ(T z) {
    z_ += z;

    normalized_ = false;
  }

  void incW(T w) {
    w_ += w;

    normalized_ = false;
  }

  void decX(T x) {
    x_ -= x;

    normalized_ = false;
  }

  void decY(T y) {
    y_ -= y;

    normalized_ = false;
  }

  void decZ(T z) {
    z_ -= z;

    normalized_ = false;
  }

  void decW(T w) {
    w_ -= w;

    normalized_ = false;
  }

  T minComponent() {
    return min(min(x_, y_), z_);
  }

  T maxComponent() {
    return max(max(x_, y_), z_);
  }

  T minAbsComponent() {
    return min(min(fabs(x_), fabs(y_)), fabs(z_));
  }

  T maxAbsComponent() {
    return max(max(fabs(x_), fabs(y_)), fabs(z_));
  }

  CNormVector3DT &operator=(const CVector3DT<T> &vector) {
    vector.getXYZ(&x_, &y_, &z_);

    normalized_ = false;

    return *this;
  }

  CNormVector3DT &operator=(const CPoint3DT<T> &point) {
    x_ = point.x; y_ = point.y; z_ = point.z;

    normalized_ = false;

    return *this;
  }

  CNormVector3DT operator+() const {
    return CNormVector3DT(x_, y_, z_, 1.0);
  }

  CNormVector3DT operator-() const {
    return CNormVector3DT(-x_, -y_, -z_, 1.0);
  }

  friend bool operator==(const CNormVector3DT &lhs,
                         const CNormVector3DT &rhs) {
    return (lhs.x_ == rhs.x_ && lhs.y_ == rhs.y_ &&
            lhs.z_ == rhs.z_ && lhs.w_ == rhs.w_);
  }

  friend bool operator!=(const CNormVector3DT &lhs,
                         const CNormVector3DT &rhs) {
    return (lhs.x_ != rhs.x_ || lhs.y_ != rhs.y_ ||
            lhs.z_ != rhs.z_ || lhs.w_ != rhs.w_);
  }

  CNormVector3DT &operator+=(const CNormVector3DT &rhs) {
    x_ += rhs.x_; y_ += rhs.y_; z_ += rhs.z_; w_  = 1.0;

    normalized_ = false;

    return *this;
  }

  CNormVector3DT operator+(const CNormVector3DT &rhs) const {
    return CNormVector3DT(x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_, 1.0);
  }

  CNormVector3DT &operator-=(const CNormVector3DT &rhs) {
    x_ -= rhs.x_; y_ -= rhs.y_; z_ -= rhs.z_; w_ = 1.0;

    normalized_ = false;

    return *this;
  }

  CNormVector3DT operator-(const CNormVector3DT &rhs) const {
    return CNormVector3DT(x_ - rhs.x_, y_ - rhs.y_, z_ - rhs.z_, 1.0);
  }

  CNormVector3DT &operator*=(T rhs) {
    x_ *= rhs; y_ *= rhs; z_ *= rhs; w_  = 1.0;

    normalized_ = false;

    return *this;
  }

  friend CNormVector3DT operator*(const CNormVector3DT &lhs, T rhs) {
    return CNormVector3DT(lhs.x_*rhs, lhs.y_*rhs, lhs.z_*rhs, 1.0);
  }

  friend CNormVector3DT operator*(T lhs, const CNormVector3DT &rhs) {
    return CNormVector3DT(lhs*rhs.x_, lhs*rhs.y_, lhs*rhs.z_, 1.0);
  }

  CNormVector3DT &operator/=(T rhs) {
    x_ /= rhs; y_ /= rhs; z_ /= rhs; w_  = 1.0;

    normalized_ = false;

    return *this;
  }

  CNormVector3DT operator/(T rhs) const {
    return CNormVector3DT(x_/rhs, y_/rhs, z_/rhs, 1.0);
  }

  friend CVector3DT<T> operator+(const CVector3DT<T> &lhs,
                                 const CNormVector3DT &rhs) {
    return CVector3DT<T>(lhs.getX() + rhs.x_, lhs.getY() + rhs.y_,
                     lhs.getZ() + rhs.z_);
  }

  friend CPoint3DT<T> operator+(const CPoint3DT<T> &lhs,
                                const CNormVector3DT &rhs) {
    return CPoint3DT<T>(lhs.x + rhs.x_, lhs.y + rhs.y_, lhs.z + rhs.z_);
  }

  friend CVector3DT<T> operator-(const CVector3DT<T> &lhs,
                                 const CNormVector3DT &rhs) {
    return CPoint3DT<T>(lhs.getX() - rhs.x_, lhs.getY() - rhs.y_,
                    lhs.getZ() - rhs.z_);
  }

  friend CPoint3DT<T> operator-(const CPoint3DT<T> &lhs,
                                const CNormVector3DT &rhs) {
    return CPoint3DT<T>(lhs.x - rhs.x_, lhs.y - rhs.y_, lhs.z - rhs.z_);
  }

  void print(ostream &os) const {
    os << "(" << x_ << "," << y_ << "," << z_ << "," << w_ << ")";
  }

  friend ostream &operator<<(ostream &os, const CNormVector3DT &vector) {
    vector.print(os);

    return os;
  }

  CNormVector3DT crossProduct(const CNormVector3DT &vector) const {
    return CNormVector3DT(y_*vector.z_ - z_*vector.y_,
                         z_*vector.x_ - x_*vector.z_,
                         x_*vector.y_ - y_*vector.x_);
  }

  CNormVector3DT crossProduct(T x, T y, T z, T z) const {
    return CNormVector3DT(y_*z - z_*y, z_*x - x_*z, x_*y - y_*x, 1.0);
  }

  static CNormVector3DT crossProduct(const CNormVector3DT &vector1,
                                     const CNormVector3DT &vector2) {
    return CNormVector3DT(vector1.y_*vector2.z_ - vector1.z_*vector2.y_,
                         vector1.z_*vector2.x_ - vector1.x_*vector2.z_,
                         vector1.x_*vector2.y_ - vector1.y_*vector2.x_);
  }

  T dotProduct(const CNormVector3DT &vector) const {
    return (x_*vector.x_ + y_*vector.y_ + z_*vector.z_);
  }

  T dotProduct(T x, T y, T z) const {
    return (x_*x + y_*y + z_*z);
  }

  static T dotProduct(const CNormVector3DT &vector1,
                      const CNormVector3DT &vector2) {
    return (vector1.x_*vector2.x_ + vector1.y_*vector2.y_ +
            vector1.z_*vector2.z_);
  }

  static T dotProduct(const CNormVector3DT &vector1,
                           T x2, T y2, T z2) {
    return (vector1.x_*x2 + vector1.y_*y2 + vector1.z_*z2);
  }

  T dotProductSelf() const {
    return (x_*x_ + y_*y_ + z_*z_);
  }

  T cosIncluded(const CNormVector3DT &vector1) const {
    return dotProduct(vector1)/(length()*vector1.length());
  }

  static T cosIncluded(const CNormVector3DT &vector1,
                       const CNormVector3DT &vector2) {
    return vector1.dotProduct(vector2)/(vector1.length()*vector2.length());
  }
};

typedef CNormVector3DT<double> CNormVector3D;

#endif
