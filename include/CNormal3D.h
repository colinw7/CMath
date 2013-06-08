#ifndef CNORMAL_3D_H
#define CNORMAL_3D_H

#include <CMathGen.h>
#include <CVector3D.h>

template<typename T>
class CVector3DT;

template<typename T>
class CNormal3DT {
 private:
  union {
    struct {
      T x_, y_, z_;
    };

    T m_[3];
  };

  bool normalized_;

 public:
  explicit CNormal3DT(T x=0.0,
                      T y=0.0,
                      T z=0.0) :
   x_(x), y_(y), z_(z), normalized_(false) {
  }

  CNormal3DT(const CNormal3DT &normal) :
    x_(normal.x_), y_(normal.y_), z_(normal.z_),
    normalized_(normal.normalized_) {
  }

  explicit CNormal3DT(const CVector3DT<T> &vector) :
    x_(vector.getX()), y_(vector.getY()), z_(vector.getZ()),
    normalized_(vector.getNormalized()) {
  }

  //------

  T getX() const { return x_; }
  T getY() const { return y_; }
  T getZ() const { return z_; }

  void getXYZ(T *x, T *y, T *z) const {
    *x = x_; *y = y_; *z = z_;
  }

  T operator[](int i) const { return m_[i]; }

  // Reference routine would break encapsulation

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

  void setXYZ(T x, T y, T z) {
    x_ = x; y_ = y; z_ = z;

    normalized_ = false;
  }

  bool getNormalized() const {
    return normalized_;
  }

  //------

  T length() const {
    if (normalized_)
      return 1.0 ;
    else
      return ::sqrt(x_*x_ + y_*y_ + z_*z_);
  }

  T modulus() const {
    return length();
  }

  T lengthSqr() const {
    if (normalized_)
      return 1.0 ;
    else
      return (x_*x_ + y_*y_ + z_*z_);
  }

  //------

  CNormal3DT &zero() {
    x_ = 0.0; y_ = 0.0; z_ = 0.0;

    normalized_ = false;

    return *this;
  }

  CNormal3DT unit() const {
    if (normalized_)
      return *this;

    T len = length();

    T factor = 0.0;

    if (len > 0.0)
      factor = 1.0 /len;

    return CNormal3DT(x_*factor, y_*factor, z_*factor, true);
  }

  //------

  CNormal3DT &normalize() {
    if (normalized_)
      return *this;

    T len = length();

    T factor = 0.0;

    if (len > 0.0)
      factor = 1.0 /len;

    x_ *= factor;
    y_ *= factor;
    z_ *= factor;

    normalized_ = true;

    return *this;
  }

  CNormal3DT normalized() const {
    if (normalized_)
      return *this;

    T len = length();

    T factor = 0.0;

    if (len > 0.0)
      factor = 1.0 /len;

    return CNormal3DT(x_*factor, y_*factor, z_*factor, true);
  }

  //------

  CNormal3DT &setMagnitude(T magnitude) {
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

    normalized_ = false;

    return *this;
  }

  //------

  T getDistance(const CNormal3DT &normal) {
    CNormal3DT diff = *this - normal;

    return diff.length();
  }

  //------

  void incX(T x = 1.0 ) {
    x_ += x;

    normalized_ = false;
  }

  void incY(T y = 1.0 ) {
    y_ += y;

    normalized_ = false;
  }

  void incZ(T z = 1.0 ) {
    z_ += z;

    normalized_ = false;
  }

  void decX(T x = 1.0 ) {
    x_ -= x;

    normalized_ = false;
  }

  void decY(T y = 1.0 ) {
    y_ -= y;

    normalized_ = false;
  }

  void decZ(T z = 1.0 ) {
    z_ -= z;

    normalized_ = false;
  }

  //------

  T minComponent() {
    return min(min(x_, y_), z_);
  }

  T maxComponent() {
    return max(max(x_, y_), z_);
  }

  T minAbsComponent() {
    return min(min(::fabs(x_), ::fabs(y_)),
               ::fabs(z_));
  }

  T maxAbsComponent() {
    return max(max(::fabs(x_), ::fabs(y_)),
               ::fabs(z_));
  }

  //------

  CNormal3DT &operator=(const CNormal3DT &normal) {
    x_ = normal.x_; y_ = normal.y_; z_ = normal.z_;

    normalized_ = normal.normalized_;

    return *this;
  }

  //------

  CNormal3DT operator+() const {
    return CNormal3DT(x_, y_, z_);
  }

  CNormal3DT operator-() const {
    return CNormal3DT(-x_, -y_, -z_);
  }

  //------

  friend bool operator==(const CNormal3DT &lhs, const CNormal3DT &rhs) {
    return (lhs.x_ == rhs.x_ && lhs.y_ == rhs.y_ && lhs.z_ == rhs.z_);
  }

  friend bool operator!=(const CNormal3DT &lhs, const CNormal3DT &rhs) {
    return (lhs.x_ != rhs.x_ || lhs.y_ != rhs.y_ || lhs.z_ != rhs.z_);
  }

  //------

  CNormal3DT &operator+=(const CNormal3DT &rhs) {
    x_ += rhs.x_; y_ += rhs.y_; z_ += rhs.z_;

    normalized_ = false;

    return *this;
  }

  CNormal3DT operator+(const CNormal3DT &rhs) const {
    return CNormal3DT(x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_);
  }

  //------

  CNormal3DT &operator-=(const CNormal3DT &rhs) {
    x_ -= rhs.x_; y_ -= rhs.y_; z_ -= rhs.z_;

    normalized_ = false;

    return *this;
  }

  CNormal3DT operator-(const CNormal3DT &rhs) const {
    return CNormal3DT(x_ - rhs.x_, y_ - rhs.y_, z_ - rhs.z_);
  }

  //------

  CNormal3DT &operator*=(T rhs) {
    x_ *= rhs; y_ *= rhs; z_ *= rhs;

    normalized_ = false;

    return *this;
  }

  friend CNormal3DT operator*(const CNormal3DT &lhs, T rhs) {
    return CNormal3DT(lhs.x_*rhs, lhs.y_*rhs, lhs.z_*rhs);
  }

  friend CNormal3DT operator*(T lhs, const CNormal3DT &rhs) {
    return CNormal3DT(lhs*rhs.x_, lhs*rhs.y_, lhs*rhs.z_);
  }

  //------

  CNormal3DT &operator/=(T rhs) {
    T irhs = 1.0 /rhs;

    x_ *= irhs; y_ *= irhs; z_ *= irhs;

    normalized_ = false;

    return *this;
  }

  CNormal3DT operator/(T rhs) const {
    T irhs = 1.0 /rhs;

    return CNormal3DT(x_*irhs, y_*irhs, z_*irhs);
  }

  //------

  void print(std::ostream &os=std::cout) const {
    os << "(" << x_ << "," << y_ << "," << z_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CNormal3DT &normal) {
    normal.print(os);

    return os;
  }

  //------

  T dotProduct(const CNormal3DT &normal) const {
    return (x_*normal.x_ + y_*normal.y_ + z_*normal.z_);
  }

  T dotProduct(const CVector3DT<T> &vector) const {
    return (x_*vector.getX() + y_*vector.getY() + z_*vector.getZ());
  }

  T dotProduct(T x, T y, T z) const {
    return (x_*x + y_*y + z_*z);
  }

  static T dotProduct(const CNormal3DT &normal1,
                      const CNormal3DT &normal2) {
    return (normal1.x_*normal2.x_ +
            normal1.y_*normal2.y_ +
            normal1.z_*normal2.z_);
  }

  static T dotProduct(const CNormal3DT &normal1,
                      const CVector3DT<T> &vector2) {
    return (normal1.x_*vector2.getX() +
            normal1.y_*vector2.getY() +
            normal1.z_*vector2.getZ());
  }

  static T dotProduct(const CVector3DT<T> &vector1,
                      const CNormal3DT &normal2) {
    return (vector1.getX()*normal2.x_ +
            vector1.getY()*normal2.y_ +
            vector1.getZ()*normal2.z_);
  }

  static T dotProduct(const CNormal3DT &normal1, T x2, T y2, T z2) {
    return (normal1.x_*x2 + normal1.y_*y2 + normal1.z_*z2);
  }

  T dotProductSelf() const {
    return (x_*x_ + y_*y_ + z_*z_);
  }

  static T absDotProduct(const CNormal3DT &normal1,
                         const CNormal3DT &normal2) {
    return ::fabs(normal1.x_*normal2.x_ +
                           normal1.y_*normal2.y_ +
                           normal1.z_*normal2.z_);
  }

  static T absDotProduct(const CNormal3DT &normal1,
                         const CVector3DT<T> &vector2) {
    return ::fabs(normal1.x_*vector2.getX() +
                           normal1.y_*vector2.getY() +
                           normal1.z_*vector2.getZ());
  }

  static T absDotProduct(const CVector3DT<T> &vector1,
                         const CNormal3DT &normal2) {
    return ::fabs(vector1.getX()*normal2.x_ +
                           vector1.getY()*normal2.y_ +
                           vector1.getZ()*normal2.z_);
  }

  //------

  void flip() {
    x_ = -x_; y_ = -y_; z_ = -z_;
  }

  //------

  T cosIncluded(const CNormal3DT &normal1) const {
    T dot = dotProduct(normal1);

    if (! normalized_)
      dot /= length();

    if (! normal1.normalized_)
      dot /= normal1.length();

    return dot;
  }

  static T cosIncluded(const CNormal3DT &normal1,
                       const CNormal3DT &normal2) {
    T dot = normal1.dotProduct(normal2);

    if (! normal1.normalized_)
      dot /= normal1.length();

    if (! normal2.normalized_)
      dot /= normal2.length();

    return dot;
  }

  //------

 private:
  CNormal3DT(T x, T y, T z, bool normalized) :
   x_(x), y_(y), z_(z), normalized_(normalized) {
  }
};

typedef CNormal3DT<double> CNormal3D;

#endif
