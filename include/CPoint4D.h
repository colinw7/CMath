#ifndef CPOINT_4D_H
#define CPOINT_4D_H

#include <iostream>

template<typename T>
class CPoint4DT {
 public:
  T x, y, z, w;

 public:
  CPoint4DT() :
   x(0), y(0), z(0), w(1) {
  }

  CPoint4DT(T x1, T y1, T z1, T w1=1) :
   x(x1), y(y1), z(z1), w(w1) {
  }

  CPoint4DT(const CPoint4DT &point) :
   x(point.x), y(point.y), z(point.z), w(point.w) {
  }

  T getX() const { return x; }
  T getY() const { return y; }
  T getZ() const { return z; }
  T getW() const { return w; }

  void getXYZW(T *x1, T *y1, T *z1, T *w1) const {
    *x1 = x; *y1 = y; *z1 = z; *w1 = w;
  }

  T operator[](int i) const {
    switch (i) {
      case 0 : return x;
      case 1 : return y;
      case 2 : return z;
      case 3 : return w;
      default: assert(false);
    }
  }

  T &operator[](int i) {
    switch (i) {
      case 0 : return x;
      case 1 : return y;
      case 2 : return z;
      case 3 : return w;
      default: assert(false);
    }
  }

  void setX(T x1) { x = x1; }
  void setY(T y1) { y = y1; }
  void setZ(T z1) { z = z1; }
  void setW(T w1) { w = w1; }

  void setXYZW(T x1, T y1, T z1, T w1) {
    x = x1; y = y1; z = z1; w = w1;
  }

  T modulus() const {
    return ::sqrt(x*x + y*y + z*z + w*w);
  }

  T modulusSquared() const {
    return (x*x + y*y + z*z + w*w);
  }

  T getDistance(const CPoint4DT &point) const {
    CPoint4DT diff = *this - point;

    return diff.modulus();
  }

  T getDistanceSquared(const CPoint4DT &point) const {
    CPoint4DT diff = *this - point;

    return diff.modulusSquared();
  }

  CPoint4DT &zero() {
    x = 0.0; y = 0.0; z = 0.0; w = 0.0;

    return *this;
  }

  //-----

  CPoint4DT operator+() {
    return CPoint4DT(x, y, z);
  }

  CPoint4DT operator-() {
    return CPoint4DT(-x, -y, -z);
  }

  //-----

  friend bool operator==(const CPoint4DT &lhs, const CPoint4DT &rhs) {
    return (lhs.x == rhs.x && lhs.y == rhs.y &&
            lhs.z == rhs.z && lhs.w == rhs.w);
  }

  friend bool operator!=(const CPoint4DT &lhs, const CPoint4DT &rhs) {
    return (lhs.x != rhs.x || lhs.y != rhs.y ||
            lhs.z != rhs.z || lhs.w != rhs.w);
  }

  //------

  CPoint4DT &operator+=(const CPoint4DT &rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z; w += rhs.w;

    return *this;
  }

  CPoint4DT operator+(const CPoint4DT &rhs) const {
    CPoint4DT t = *this;

    t += rhs;

    return t;
  }

  //------

  CPoint4DT &operator-=(const CPoint4DT &rhs) {
    x -= rhs.x; y -= rhs.y; z -= rhs.z; w -= rhs.w;

    return *this;
  }

  CPoint4DT operator-(const CPoint4DT &rhs) const {
    CPoint4DT t = *this;

    t -= rhs;

    return t;
  }

  //------

  CPoint4DT &operator*=(T rhs) {
    x *= rhs; y *= rhs; z *= rhs; w *= rhs;

    return *this;
  }

  CPoint4DT operator*(T rhs) const {
    CPoint4DT t = *this;

    t *= rhs;

    return t;
  }

  //------

  CPoint4DT &operator/=(T rhs) {
    T irhs = 1.0/rhs;

    x *= irhs; y *= irhs; z *= irhs; w *= irhs;

    return *this;
  }

  CPoint4DT operator/(T rhs) const {
    CPoint4DT t = *this;

    t /= rhs;

    return t;
  }

  void print(std::ostream &os) const {
    return os << "(" << x << "," << y << "," << z << "," << w << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CPoint4DT &point) {
    point.print(os);

    return os;
  }
};

typedef CPoint4DT<double> CPoint4D;

#endif
