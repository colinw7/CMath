#ifndef CPOINT_4D_H
#define CPOINT_4D_H

#include <cmath>
#include <iostream>
#include <cassert>

class CPoint4D {
 public:
  CPoint4D() :
   x(0), y(0), z(0), w(1) {
  }

  CPoint4D(double x1, double y1, double z1, double w1=1) :
   x(x1), y(y1), z(z1), w(w1) {
  }

  CPoint4D(const CPoint4D &point) :
   x(point.x), y(point.y), z(point.z), w(point.w) {
  }

  double getX() const { return x; }
  double getY() const { return y; }
  double getZ() const { return z; }
  double getW() const { return w; }

  void getXYZW(double *x1, double *y1, double *z1, double *w1) const {
    *x1 = x; *y1 = y; *z1 = z; *w1 = w;
  }

  double operator[](int i) const {
    switch (i) {
      case 0 : return x;
      case 1 : return y;
      case 2 : return z;
      case 3 : return w;
      default: assert(false);
    }
  }

  double &operator[](int i) {
    switch (i) {
      case 0 : return x;
      case 1 : return y;
      case 2 : return z;
      case 3 : return w;
      default: assert(false);
    }
  }

  void setX(double x1) { x = x1; }
  void setY(double y1) { y = y1; }
  void setZ(double z1) { z = z1; }
  void setW(double w1) { w = w1; }

  void setXYZW(double x1, double y1, double z1, double w1) {
    x = x1; y = y1; z = z1; w = w1;
  }

  double modulus() const {
    return ::sqrt(x*x + y*y + z*z + w*w);
  }

  double modulusSquared() const {
    return (x*x + y*y + z*z + w*w);
  }

  double getDistance(const CPoint4D &point) const {
    CPoint4D diff = *this - point;

    return diff.modulus();
  }

  double getDistanceSquared(const CPoint4D &point) const {
    CPoint4D diff = *this - point;

    return diff.modulusSquared();
  }

  CPoint4D &zero() {
    x = 0.0; y = 0.0; z = 0.0; w = 0.0;

    return *this;
  }

  //-----

  CPoint4D operator+() {
    return CPoint4D(x, y, z);
  }

  CPoint4D operator-() {
    return CPoint4D(-x, -y, -z);
  }

  //-----

  friend bool operator==(const CPoint4D &lhs, const CPoint4D &rhs) {
    return (lhs.x == rhs.x && lhs.y == rhs.y &&
            lhs.z == rhs.z && lhs.w == rhs.w);
  }

  friend bool operator!=(const CPoint4D &lhs, const CPoint4D &rhs) {
    return (lhs.x != rhs.x || lhs.y != rhs.y ||
            lhs.z != rhs.z || lhs.w != rhs.w);
  }

  //------

  CPoint4D &operator+=(const CPoint4D &rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z; w += rhs.w;

    return *this;
  }

  CPoint4D operator+(const CPoint4D &rhs) const {
    CPoint4D t = *this;

    t += rhs;

    return t;
  }

  //------

  CPoint4D &operator-=(const CPoint4D &rhs) {
    x -= rhs.x; y -= rhs.y; z -= rhs.z; w -= rhs.w;

    return *this;
  }

  CPoint4D operator-(const CPoint4D &rhs) const {
    CPoint4D t = *this;

    t -= rhs;

    return t;
  }

  //------

  CPoint4D &operator*=(double rhs) {
    x *= rhs; y *= rhs; z *= rhs; w *= rhs;

    return *this;
  }

  CPoint4D operator*(double rhs) const {
    CPoint4D t = *this;

    t *= rhs;

    return t;
  }

  //------

  CPoint4D &operator/=(double rhs) {
    double irhs = 1.0/rhs;

    x *= irhs; y *= irhs; z *= irhs; w *= irhs;

    return *this;
  }

  CPoint4D operator/(double rhs) const {
    CPoint4D t = *this;

    t /= rhs;

    return t;
  }

  void print(std::ostream &os) const {
    os << "(" << x << "," << y << "," << z << "," << w << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CPoint4D &point) {
    point.print(os);

    return os;
  }

 public:
  double x { 0 }, y { 0 }, z { 0 }, w { 1 };
};

#endif
