#ifndef CPOINT_2D_H
#define CPOINT_2D_H

#include <CIPoint2D.h>

#include <cmath>
#include <cassert>

template<typename T>
class CPoint2DT {
 private:
  typedef CPoint2DT<T>    Point;
  typedef CIPoint2DT<int> IPoint;

 public:
   T x, y;

  CPoint2DT() :
   x(0), y(0) {
  }

  CPoint2DT(T x1, T y1) :
   x(x1), y(y1) {
  }

  CPoint2DT(const CPoint2DT &point) :
   x(point.x), y(point.y) {
  }

  CPoint2DT(const IPoint &point) :
   x(point.x), y(point.y) {
  }

  IPoint toIPoint() const { return IPoint(int(x), int(y)); }

  T getX() const { return x; }
  T getY() const { return y; }

  void getXY(T *x1, T *y1) const {
    *x1 = x; *y1 = y;
  }

  void getXY(T xy[2]) const {
    xy[0] = x; xy[1] = y;
  }

  T operator[](int i) const {
    switch (i) {
      case 0 : return x;
      case 1 : return y;
      default: assert(false);
    }
  }

  T &operator[](int i) {
    switch (i) {
      case 0 : return x;
      case 1 : return y;
      default: assert(false); }
  }

  void setX(T x1) { x = x1; }
  void setY(T y1) { y = y1; }

  void setXY(T x1, T y1) {
    x = x1; y = y1;
  }

  void setXY(T xy[2]) {
    x = xy[0]; y = xy[1];
  }

  //------

  Point &zero() {
    x = 0.0; y = 0.0;

    return *this;
  }

  //-----

  Point operator+() const {
    return Point(x, y);
  }

  Point operator-() const {
    return Point(-x, -y);
  }

  //-----

  bool equal(const Point &rhs, double tol=1E-6) const {
    T dx = fabs(x - rhs.x);
    T dy = fabs(y - rhs.y);

    return (dx < tol && dy < tol);
  }

  //-----

  friend bool operator==(const Point &lhs, const Point &rhs) {
    return lhs.equal(rhs, 1E-6);
  }

  friend bool operator!=(const Point &lhs, const Point &rhs) {
    return ! (lhs == rhs);
  }

  //------

  // Addition of points makes no mathematical sense but
  // is useful for weighted sum

  Point &operator+=(const Point &rhs) {
    x += rhs.x; y += rhs.y;

    return *this;
  }

  Point &operator+=(T rhs) {
    x += rhs; y += rhs;

    return *this;
  }

  Point operator+(const Point &rhs) const {
    return Point(x + rhs.x, y + rhs.y);
  }

  friend Point operator+(const Point &lhs, T rhs) {
    return Point(lhs.x + rhs, lhs.y + rhs);
  }

  friend Point operator+(T lhs, const Point &rhs) {
    return Point(rhs.x + lhs, rhs.y + lhs);
  }

  //------

  // Subtraction of points makes no mathematical sense but
  // is useful for weighted sum

  Point &operator-=(const Point &rhs) {
    x -= rhs.x; y -= rhs.y;

    return *this;
  }

  Point operator-(const Point &rhs) const {
    return Point(x - rhs.x, y - rhs.y);
  }

  //------

  // Multiplication of points makes no mathematical sense but
  // is useful for weighted sum

  Point &operator*=(T rhs) {
    x *= rhs; y *= rhs;

    return *this;
  }

  Point &operator*=(const Point &rhs) {
    x *= rhs.x; y *= rhs.y;

    return *this;
  }

  Point operator*(const Point &rhs) const {
    return Point(x*rhs.x, y*rhs.y);
  }

  friend Point operator*(const Point &lhs, T rhs) {
    return Point(lhs.x*rhs, lhs.y*rhs);
  }

  friend Point operator*(T lhs, const Point &rhs) {
    return Point(rhs.x*lhs, rhs.y*lhs);
  }

  //------

  // Division of points makes no mathematical sense but
  // is useful for weighted sum

  Point &operator/=(T rhs) {
    T irhs = 1.0/rhs;

    x *= irhs; y *= irhs;

    return *this;
  }

  Point &operator/=(const Point &rhs) {
    x /= rhs.x; y /= rhs.y;

    return *this;
  }

  Point operator/(const Point &rhs) const {
    return Point(x/rhs.x, y/rhs.y);
  }

  friend Point operator/(const Point &lhs, T rhs) {
    T irhs = 1.0/rhs;

    return Point(lhs.x*irhs, lhs.y*irhs);
  }

  friend Point operator/(T lhs, const Point &rhs) {
    return Point(lhs/rhs.x, lhs/rhs.y);
  }

  //------

  T minComponent() const {
    return std::min(x, y);
  }

  T maxComponent() const {
    return std::max(x, y);
  }

  //-----

  T distanceSqrTo(const Point &rhs) const {
    T dx = x - rhs.x;
    T dy = y - rhs.y;

    return (dx*dx + dy*dy);
  }

  T distanceTo(const Point &rhs) const {
    return sqrt(distanceSqrTo(rhs));
  }

  //-----

  Point rotate(const Point &center, T da) const {
    T s = sin(da);
    T c = cos(da);

    T x1 = x - center.x;
    T y1 = y - center.y;

    T x2 = x1*c - y1*s;
    T y2 = x1*s + y1*c;

    return Point(x2 + center.x, y2 + center.y);
  }

  //-----

  Point flip(const Point &c, bool x_axis=true) const {
    Point p = *this;

    if (x_axis)
      p.x = 2*c.x - p.x;
    else
      p.y = 2*c.y - p.y;

    return p;
  }

  //-----

  static Point min(const Point &lhs, const Point &rhs) {
    return Point(std::min(lhs.x, rhs.x), std::min(lhs.y, rhs.y));
  }

  static Point max(const Point &lhs, const Point &rhs) {
    return Point(std::max(lhs.x, rhs.x), std::max(lhs.y, rhs.y));
  }

  //-----

  void print(std::ostream &os) const {
    os << "(" << x << "," << y << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const Point &point) {
    point.print(os);

    return os;
  }
};

typedef CPoint2DT<double> CPoint2D;

#endif
