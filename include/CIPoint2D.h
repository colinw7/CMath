#ifndef CIPOINT_2D_H
#define CIPOINT_2D_H

#include <iostream>

template<typename T>
class CIPoint2DT {
 public:
  T x, y;

  CIPoint2DT(T x1=0, T y1=0) :
   x(x1), y(y1) {
  }

  CIPoint2DT(const CIPoint2DT &point) :
   x(point.x), y(point.y) {
  }

  const CIPoint2DT &operator=(const CIPoint2DT &point) {
    x = point.x;
    y = point.y;

    return *this;
  }

  T getX() const { return x; }
  T getY() const { return y; }

  void setX(T x1) { x = x1; }
  void setY(T y1) { y = y1; }

  void setXY(T x1, T y1) {
    x = x1; y = y1;
  }

  void getXY(T *x1, T *y1) const {
    *x1 = x; *y1 = y;
  }

  CIPoint2DT &zero() {
    x = 0; y = 0;

    return *this;
  }

  double modulus() const {
    return sqrt(x*x + y*y);
  }

  CIPoint2DT &operator+=(const CIPoint2DT &rhs) {
    x += rhs.x;
    y += rhs.y;

    return *this;
  }

  CIPoint2DT operator+(const CIPoint2DT &rhs) const {
    CIPoint2DT t = *this;

    t += rhs;

    return t;
  }

  CIPoint2DT &operator-=(const CIPoint2DT &rhs) {
    x -= rhs.x;
    y -= rhs.y;

    return *this;
  }

  CIPoint2DT operator-(const CIPoint2DT &rhs) const {
    CIPoint2DT t = *this;

    t -= rhs;

    return t;
  }

  CIPoint2DT &operator*=(T rhs) {
    x *= rhs;
    y *= rhs;

    return *this;
  }

  friend CIPoint2DT operator*(const CIPoint2DT &lhs, T rhs) {
    CIPoint2DT t = lhs;

    t *= rhs;

    return t;
  }

  friend CIPoint2DT operator*(T lhs, const CIPoint2DT &rhs) {
    CIPoint2DT t = rhs;

    t *= lhs;

    return t;
  }

  CIPoint2DT &operator/=(T rhs) {
    x /= rhs;
    y /= rhs;

    return *this;
  }

  friend CIPoint2DT operator/(const CIPoint2DT &lhs, T rhs) {
    CIPoint2DT t = lhs;

    t /= rhs;

    return t;
  }

  friend CIPoint2DT operator/(T lhs, const CIPoint2DT &rhs) {
    CIPoint2DT t = rhs;

    t /= lhs;

    return t;
  }

  //-----

  int cmp(const CIPoint2DT &rhs) const {
    int dx = x - rhs.x;

    if (dx != 0) return dx;

    int dy = y - rhs.y;

    if (dy != 0) return dy;

    return 0;
  }

  //-----

  T minComponent() {
    return std::min(x, y);
  }

  T maxComponent() {
    return std::max(x, y);
  }

  //-----

  static CIPoint2DT min(const CIPoint2DT &lhs, const CIPoint2DT &rhs) {
    return CIPoint2DT(std::min(lhs.x, rhs.x),
                      std::min(lhs.y, rhs.y));
  }

  static CIPoint2DT max(const CIPoint2DT &lhs, const CIPoint2DT &rhs) {
    return CIPoint2DT(std::max(lhs.x, rhs.x),
                      std::max(lhs.y, rhs.y));
  }

  //-----

  void print(std::ostream &os) const {
    os << "(" << x << "," << y << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CIPoint2DT &point) {
    point.print(os);

    return os;
  }
};

typedef CIPoint2DT<int> CIPoint2D;

//------

template<typename T>
inline bool operator==(const CIPoint2DT<T> &lhs, const CIPoint2DT<T> &rhs) {
  return (lhs.cmp(rhs) == 0);
}

template<typename T>
inline bool operator!=(const CIPoint2DT<T> &lhs, const CIPoint2DT<T> &rhs) {
  return ! operator==(lhs, rhs);
}

template<typename T>
inline bool operator<(const CIPoint2DT<T> &lhs, const CIPoint2DT<T> &rhs) {
  return (lhs.cmp(rhs) < 0);
}

template<typename T>
inline bool operator>(const CIPoint2DT<T> &lhs, const CIPoint2DT<T> &rhs) {
  return (lhs.cmp(rhs) > 0);
}

template<typename T>
inline bool operator>=(const CIPoint2DT<T> &lhs, const CIPoint2DT<T> &rhs) {
  return ! operator<(lhs, rhs);
}

template<typename T>
inline bool operator<=(const CIPoint2DT<T> &lhs, const CIPoint2DT<T> &rhs) {
  return ! operator>(lhs, rhs);
}

#endif
