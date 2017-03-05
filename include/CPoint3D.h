#ifndef CPOINT_3D_H
#define CPOINT_3D_H

#include <cassert>
#include <iostream>
#include <cmath>
#include <CPoint2D.h>
#include <CIPoint3D.h>

class CPoint3D {
 public:
  CPoint3D() :
   x(0), y(0), z(0) {
  }

  CPoint3D(double x1, double y1, double z1) :
   x(x1), y(y1), z(z1) {
  }

  CPoint3D(const CPoint3D &point) :
   x(point.x), y(point.y), z(point.z) {
  }

  CPoint3D(const CPoint2D &point, double z=0) :
   x(point.x), y(point.y), z(z) {
  }

  CPoint2D toPoint2D() const { return CPoint2D(x, y); }

  double getX() const { return x; }
  double getY() const { return y; }
  double getZ() const { return z; }

  void getXYZ(double *x1, double *y1, double *z1) const {
    *x1 = x; *y1 = y; *z1 = z;
  }

  void getXYZ(double xyz[3]) const {
    xyz[0] = x; xyz[1] = y; xyz[2] = z;
  }

  double operator[](int i) const {
    switch (i) {
      case 0 : return x;
      case 1 : return y;
      case 2 : return z;
      default: assert(false);
    }
  }

  double &operator[](int i) {
    switch (i) {
      case 0 : return x;
      case 1 : return y;
      case 2 : return z;
      default: assert(false);
    }
  }

  CPoint3D &setX(double x1) { x = x1; return *this; }
  CPoint3D &setY(double y1) { y = y1; return *this; }
  CPoint3D &setZ(double z1) { z = z1; return *this; }

  CPoint3D &setXYZ(double x1, double y1, double z1) {
    x = x1; y = y1; z = z1;
    return *this;
  }

  CPoint3D &setXYZ(double xyz[3]) {
    x = xyz[0]; y = xyz[1]; z = xyz[2];
    return *this;
  }

  //------

  CPoint3D &zero() {
    x = 0.0; y = 0.0; z = 0.0;

    return *this;
  }

  //-----

  CPoint3D operator+() const {
    return CPoint3D(x, y, z);
  }

  CPoint3D operator-() const {
    return CPoint3D(-x, -y, -z);
  }

  //-----

  friend bool operator==(const CPoint3D &lhs, const CPoint3D &rhs) {
    return isEqual(lhs, rhs);
  }

  friend bool operator!=(const CPoint3D &lhs, const CPoint3D &rhs) {
    return ! (lhs == rhs);
  }

  //------

  static bool isEqual(const CPoint3D &lhs, const CPoint3D &rhs, double tol=1E-6) {
    double dx = fabs(lhs.x - rhs.x);
    double dy = fabs(lhs.y - rhs.y);
    double dz = fabs(lhs.z - rhs.z);

    return (dx < tol && dy < tol && dz < tol);
  }

  //------

  // Addition of points makes no mathematical sense but
  // is useful for weighted sum

  CPoint3D &operator+=(const CPoint3D &rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z;

    return *this;
  }

  CPoint3D &operator+=(double rhs) {
    x += rhs; y += rhs; z += rhs;

    return *this;
  }

  CPoint3D operator+(const CPoint3D &rhs) const {
    return CPoint3D(x + rhs.x, y + rhs.y, z + rhs.z);
  }

  friend CPoint3D operator+(const CPoint3D &lhs, double rhs) {
    return CPoint3D(lhs.x + rhs, lhs.y + rhs, lhs.z + rhs);
  }

  friend CPoint3D operator+(double lhs, const CPoint3D &rhs) {
    return CPoint3D(rhs.x + lhs, rhs.y + lhs, rhs.z + lhs);
  }

  //------

  // Subtraction of points makes no mathematical sense but
  // is useful for weighted sum

  CPoint3D &operator-=(const CPoint3D &rhs) {
    x -= rhs.x; y -= rhs.y; z -= rhs.z;

    return *this;
  }

  CPoint3D operator-(const CPoint3D &rhs) const {
    return CPoint3D(x - rhs.x, y - rhs.y, z - rhs.z);
  }

  //-----

  // Multiplication of points makes no mathematical sense but
  // is useful for weighted sum

  CPoint3D &operator*=(double rhs) {
    x *= rhs; y *= rhs; z *= rhs;

    return *this;
  }

  CPoint3D &operator*=(const CPoint3D &rhs) {
    x *= rhs.x; y *= rhs.y; z *= rhs.z;

    return *this;
  }

  CPoint3D operator*(const CPoint3D &rhs) const {
    return CPoint3D(x*rhs.x, y*rhs.y, z*rhs.z);
  }

  friend CPoint3D operator*(const CPoint3D &lhs, double rhs) {
    return CPoint3D(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs);
  }

  friend CPoint3D operator*(double lhs, const CPoint3D &rhs) {
    return CPoint3D(rhs.x*lhs, rhs.y*lhs, rhs.z*lhs);
  }

  //------

  // Division of points makes no mathematical sense but
  // is useful for weighted sum

  CPoint3D &operator/=(double rhs) {
    double irhs = 1.0/rhs;

    x *= irhs; y *= irhs; z *= irhs;

    return *this;
  }

  CPoint3D &operator/=(const CPoint3D &rhs) {
    x /= rhs.x; y /= rhs.y; z /= rhs.z;

    return *this;
  }

  CPoint3D operator/(const CPoint3D &rhs) const {
    return CPoint3D(x/rhs.x, y/rhs.y, z/rhs.z);
  }

  friend CPoint3D operator/(const CPoint3D &lhs, double rhs) {
    double irhs = 1.0/rhs;

    return CPoint3D(lhs.x*irhs, lhs.y*irhs, lhs.z*irhs);
  }

  friend CPoint3D operator/(double lhs, const CPoint3D &rhs) {
    return CPoint3D(lhs/rhs.x, lhs/rhs.y, lhs/rhs.z);
  }

  //-----

  double minComponent() {
    return std::min(x, std::min(y, z));
  }

  double maxComponent() {
    return std::max(x, std::max(y, z));
  }

  //-----

  static CPoint3D min(const CPoint3D &lhs, const CPoint3D &rhs) {
    return CPoint3D(std::min(lhs.x, rhs.x), std::min(lhs.y, rhs.y), std::min(lhs.z, rhs.z));
  }

  static CPoint3D max(const CPoint3D &lhs, const CPoint3D &rhs) {
    return CPoint3D(std::max(lhs.x, rhs.x), std::max(lhs.y, rhs.y), std::max(lhs.z, rhs.z));
  }

  static CPoint3D mid(const CPoint3D &lhs, const CPoint3D &rhs) {
    return CPoint3D(0.5*(lhs.x + rhs.x), 0.5*(lhs.y + rhs.y), 0.5*(lhs.z + rhs.z));
  }

  //------

  // comparison
  int cmp(const CPoint3D &v) const {
    if      (x < v.x) return -1;
    else if (x > v.x) return  1;
    else if (y < v.y) return -1;
    else if (y > v.y) return  1;
    else if (z < v.z) return -1;
    else if (z > v.z) return  1;
    else              return  0;
  }

  friend bool operator< (const CPoint3D &lhs, const CPoint3D &rhs) {
    return lhs.cmp(rhs) <  0;
  }

  friend bool operator<=(const CPoint3D &lhs, const CPoint3D &rhs) {
    return lhs.cmp(rhs) <= 0;
  }

  friend bool operator> (const CPoint3D &lhs, const CPoint3D &rhs) {
    return lhs.cmp(rhs) >  0;
  }

  friend bool operator>=(const CPoint3D &lhs, const CPoint3D &rhs) {
    return lhs.cmp(rhs) >= 0;
  }

  //-----

  void print(std::ostream &os=std::cout) const {
    os << "(" << x << "," << y << "," << z << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CPoint3D &point) {
    point.print(os);

    return os;
  }

 public:
  double x { 0 }, y { 0 }, z { 0 };
};

#endif
