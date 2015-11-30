#ifndef CPOINT_3D_H
#define CPOINT_3D_H

#include <cassert>
#include <iostream>
#include <cmath>
#include <CPoint2D.h>
#include <CIPoint3D.h>

template<typename T>
class CPoint3DT {
 private:
  typedef CIPoint3DT<int> IPoint;
  typedef CPoint2DT<T>    Point2D;

 public:
  T x, y, z;

 public:
  CPoint3DT() :
   x(0), y(0), z(0) {
  }

  CPoint3DT(T x1, T y1, T z1) :
   x(x1), y(y1), z(z1) {
  }

  CPoint3DT(const CPoint3DT &point) :
   x(point.x), y(point.y), z(point.z) {
  }

  CPoint3DT(const IPoint &point) :
   x(point.x), y(point.y), z(point.z) {
  }

  CPoint3DT(const Point2D &point, double z=0) :
   x(point.x), y(point.y), z(z) {
  }

  IPoint toIPoint() const { return IPoint(int(x), int(y), int(z)); }

  Point2D toPoint2D() const { return Point2D(x, y); }

  T getX() const { return x; }
  T getY() const { return y; }
  T getZ() const { return z; }

  void getXYZ(T *x1, T *y1, T *z1) const {
    *x1 = x; *y1 = y; *z1 = z;
  }

  void getXYZ(T xyz[3]) const {
    xyz[0] = x; xyz[1] = y; xyz[2] = z;
  }

  T operator[](int i) const {
    switch (i) {
      case 0 : return x;
      case 1 : return y;
      case 2 : return z;
      default: assert(false);
    }
  }

  T &operator[](int i) {
    switch (i) {
      case 0 : return x;
      case 1 : return y;
      case 2 : return z;
      default: assert(false);
    }
  }

  void setX(T x1) { x = x1; }
  void setY(T y1) { y = y1; }
  void setZ(T z1) { z = z1; }

  void setXYZ(T x1, T y1, T z1) {
    x = x1; y = y1; z = z1;
  }

  void setXYZ(T xyz[3]) {
    x = xyz[0]; y = xyz[1]; z = xyz[2];
  }

  //------

  CPoint3DT &zero() {
    x = 0.0; y = 0.0; z = 0.0;

    return *this;
  }

  //-----

  CPoint3DT operator+() const {
    return CPoint3DT(x, y, z);
  }

  CPoint3DT operator-() const {
    return CPoint3DT(-x, -y, -z);
  }

  //-----

  friend bool operator==(const CPoint3DT &lhs, const CPoint3DT &rhs) {
    return isEqual(lhs, rhs);
  }

  friend bool operator!=(const CPoint3DT &lhs, const CPoint3DT &rhs) {
    return ! (lhs == rhs);
  }

  //------

  static bool isEqual(const CPoint3DT &lhs, const CPoint3DT &rhs, double tol=1E-6) {
    T dx = fabs(lhs.x - rhs.x);
    T dy = fabs(lhs.y - rhs.y);
    T dz = fabs(lhs.z - rhs.z);

    return (dx < tol && dy < tol && dz < tol);
  }

  //------

  // Addition of points makes no mathematical sense but
  // is useful for weighted sum

  CPoint3DT &operator+=(const CPoint3DT &rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z;

    return *this;
  }

  CPoint3DT &operator+=(T rhs) {
    x += rhs; y += rhs; z += rhs;

    return *this;
  }

  CPoint3DT operator+(const CPoint3DT &rhs) const {
    return CPoint3DT(x + rhs.x, y + rhs.y, z + rhs.z);
  }

  friend CPoint3DT operator+(const CPoint3DT &lhs, T rhs) {
    return CPoint3DT(lhs.x + rhs, lhs.y + rhs, lhs.z + rhs);
  }

  friend CPoint3DT operator+(T lhs, const CPoint3DT &rhs) {
    return CPoint3DT(rhs.x + lhs, rhs.y + lhs, rhs.z + lhs);
  }

  //------

  // Subtraction of points makes no mathematical sense but
  // is useful for weighted sum

  CPoint3DT &operator-=(const CPoint3DT &rhs) {
    x -= rhs.x; y -= rhs.y; z -= rhs.z;

    return *this;
  }

  CPoint3DT operator-(const CPoint3DT &rhs) const {
    return CPoint3DT(x - rhs.x, y - rhs.y, z - rhs.z);
  }

  //-----

  // Multiplication of points makes no mathematical sense but
  // is useful for weighted sum

  CPoint3DT &operator*=(T rhs) {
    x *= rhs; y *= rhs; z *= rhs;

    return *this;
  }

  CPoint3DT &operator*=(const CPoint3DT &rhs) {
    x *= rhs.x; y *= rhs.y; z *= rhs.z;

    return *this;
  }

  CPoint3DT operator*(const CPoint3DT &rhs) const {
    return CPoint3DT(x*rhs.x, y*rhs.y, z*rhs.z);
  }

  friend CPoint3DT operator*(const CPoint3DT &lhs, T rhs) {
    return CPoint3DT(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs);
  }

  friend CPoint3DT operator*(T lhs, const CPoint3DT &rhs) {
    return CPoint3DT(rhs.x*lhs, rhs.y*lhs, rhs.z*lhs);
  }

  //------

  // Division of points makes no mathematical sense but
  // is useful for weighted sum

  CPoint3DT &operator/=(T rhs) {
    T irhs = 1.0/rhs;

    x *= irhs; y *= irhs; z *= irhs;

    return *this;
  }

  CPoint3DT &operator/=(const CPoint3DT &rhs) {
    x /= rhs.x; y /= rhs.y; z /= rhs.z;

    return *this;
  }

  CPoint3DT operator/(const CPoint3DT &rhs) const {
    return CPoint3DT(x/rhs.x, y/rhs.y, z/rhs.z);
  }

  friend CPoint3DT operator/(const CPoint3DT &lhs, T rhs) {
    T irhs = 1.0/rhs;

    return CPoint3DT(lhs.x*irhs, lhs.y*irhs, lhs.z*irhs);
  }

  friend CPoint3DT operator/(T lhs, const CPoint3DT &rhs) {
    return CPoint3DT(lhs/rhs.x, lhs/rhs.y, lhs/rhs.z);
  }

  //-----

  T minComponent() {
    return std::min(x, std::min(y, z));
  }

  T maxComponent() {
    return std::max(x, std::max(y, z));
  }

  //-----

  static CPoint3DT min(const CPoint3DT &lhs, const CPoint3DT &rhs) {
    return CPoint3DT(std::min(lhs.x, rhs.x), std::min(lhs.y, rhs.y), std::min(lhs.z, rhs.z));
  }

  static CPoint3DT max(const CPoint3DT &lhs, const CPoint3DT &rhs) {
    return CPoint3DT(std::max(lhs.x, rhs.x), std::max(lhs.y, rhs.y), std::max(lhs.z, rhs.z));
  }

  static CPoint3DT mid(const CPoint3DT &lhs, const CPoint3DT &rhs) {
    return CPoint3DT(0.5*(lhs.x + rhs.x), 0.5*(lhs.y + rhs.y), 0.5*(lhs.z + rhs.z));
  }

  //-----

  void print(std::ostream &os=std::cout) const {
    os << "(" << x << "," << y << "," << z << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CPoint3DT &point) {
    point.print(os);

    return os;
  }
};

typedef CPoint3DT<double> CPoint3D;

#endif
