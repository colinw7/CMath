#ifndef CIPOINT_2D_H
#define CIPOINT_2D_H

#include <iostream>
#include <cmath>

class CIPoint2D {
 public:
  explicit CIPoint2D(int x1=0, int y1=0) :
   x(x1), y(y1) {
  }

  CIPoint2D(const CIPoint2D &point) :
   x(point.x), y(point.y) {
  }

  CIPoint2D &operator=(const CIPoint2D &point) {
    x = point.x;
    y = point.y;

    return *this;
  }

  int getX() const { return x; }
  int getY() const { return y; }

  void setX(int x1) { x = x1; }
  void setY(int y1) { y = y1; }

  void setXY(int x1, int y1) {
    x = x1; y = y1;
  }

  void getXY(int *x1, int *y1) const {
    *x1 = x; *y1 = y;
  }

  CIPoint2D &zero() {
    x = 0; y = 0;

    return *this;
  }

  double modulus() const {
    return std::sqrt(x*x + y*y);
  }

  CIPoint2D &operator+=(const CIPoint2D &rhs) {
    x += rhs.x;
    y += rhs.y;

    return *this;
  }

  CIPoint2D &operator+=(int rhs) {
    x += rhs;
    y += rhs;

    return *this;
  }

  CIPoint2D operator+(const CIPoint2D &rhs) const {
    CIPoint2D t = *this;

    t += rhs;

    return t;
  }

  CIPoint2D &operator-=(const CIPoint2D &rhs) {
    x -= rhs.x;
    y -= rhs.y;

    return *this;
  }

  CIPoint2D operator-(const CIPoint2D &rhs) const {
    CIPoint2D t = *this;

    t -= rhs;

    return t;
  }

  CIPoint2D &operator*=(int rhs) {
    x *= rhs;
    y *= rhs;

    return *this;
  }

  friend CIPoint2D operator*(const CIPoint2D &lhs, int rhs) {
    CIPoint2D t = lhs;

    t *= rhs;

    return t;
  }

  friend CIPoint2D operator*(int lhs, const CIPoint2D &rhs) {
    CIPoint2D t = rhs;

    t *= lhs;

    return t;
  }

  CIPoint2D &operator/=(int rhs) {
    x /= rhs;
    y /= rhs;

    return *this;
  }

  friend CIPoint2D operator/(const CIPoint2D &lhs, int rhs) {
    CIPoint2D t = lhs;

    t /= rhs;

    return t;
  }

  friend CIPoint2D operator/(int lhs, const CIPoint2D &rhs) {
    CIPoint2D t = rhs;

    t /= lhs;

    return t;
  }

  //-----

  int cmp(const CIPoint2D &rhs) const {
    int dx = x - rhs.x;

    if (dx != 0) return dx;

    int dy = y - rhs.y;

    if (dy != 0) return dy;

    return 0;
  }

  //-----

  int minComponent() {
    return std::min(x, y);
  }

  int maxComponent() {
    return std::max(x, y);
  }

  //-----

  static CIPoint2D min(const CIPoint2D &lhs, const CIPoint2D &rhs) {
    return CIPoint2D(std::min(lhs.x, rhs.x),
                      std::min(lhs.y, rhs.y));
  }

  static CIPoint2D max(const CIPoint2D &lhs, const CIPoint2D &rhs) {
    return CIPoint2D(std::max(lhs.x, rhs.x),
                      std::max(lhs.y, rhs.y));
  }

  //-----

  void print(std::ostream &os) const {
    os << "(" << x << "," << y << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CIPoint2D &point) {
    point.print(os);

    return os;
  }

 public:
  int x { 0 }, y { 0 };
};

//------

inline bool operator==(const CIPoint2D &lhs, const CIPoint2D &rhs) {
  return (lhs.cmp(rhs) == 0);
}

inline bool operator!=(const CIPoint2D &lhs, const CIPoint2D &rhs) {
  return ! operator==(lhs, rhs);
}

inline bool operator<(const CIPoint2D &lhs, const CIPoint2D &rhs) {
  return (lhs.cmp(rhs) < 0);
}

inline bool operator>(const CIPoint2D &lhs, const CIPoint2D &rhs) {
  return (lhs.cmp(rhs) > 0);
}

inline bool operator>=(const CIPoint2D &lhs, const CIPoint2D &rhs) {
  return ! operator<(lhs, rhs);
}

inline bool operator<=(const CIPoint2D &lhs, const CIPoint2D &rhs) {
  return ! operator>(lhs, rhs);
}

#endif
