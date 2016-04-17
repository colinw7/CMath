#ifndef CPOINT_3D_H
#define CPOINT_3D_H

// 3D Homogeneous point

#include <CPoint3D.h>

template<typename T>
class CPoint3DHT {
  typedef CPoint3DT<T> Point;

 public:
  union {
    T m[4];

    struct {
      T x, y, z, w;
    };
  };

 public:
  CPoint3DHT() { }

  CPoint3DHT(T x1, T y1, T z1, T w1) :
   x(x1), y(y1), z(z1), w(w1) {
  }

  CPoint3DHT(const CPoint3DHT &point) :
   x(point.x), y(point.y), z(point.z), w(point.w) {
  }

  T getX() const { return x; }
  T getY() const { return y; }
  T getZ() const { return z; }
  T getW() const { return w; }

  void getXYZW(T *x1, T *y1, T *z1, T *w1) const {
    *x1 = x; *y1 = y; *z1 = z; *w1 = w;
  }

  Point getPoint() const {
    return Point(x/w, y/w, z/w);
  }

  T operator[](int i) const { return m[i]; }

  T &operator[](int i) { return m[i]; }

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

  T getDistance(const CPoint3DHT &point) const {
    CPoint3DHT diff = *this - point;

    return diff.modulus();
  }

  T getDistanceSquared(const CPoint3DHT &point) const {
    CPoint3DHT diff = *this - point;

    return diff.modulusSquared();
  }

  CPoint3DHT &zero() {
    x = 0.0; y = 0.0; z = 0.0; w = 0.0;

    return *this;
  }

  //-----

  CPoint3DHT operator+() {
    return CPoint3DHT(x, y, z);
  }

  CPoint3DHT operator-() {
    return CPoint3DHT(-x, -y, -z);
  }

  //-----

  friend bool operator==(const CPoint3DHT &lhs, const CPoint3DHT &rhs) {
    return (lhs.x == rhs.x && lhs.y == rhs.y &&
            lhs.z == rhs.z && lhs.w == rhs.w);
  }

  friend bool operator!=(const CPoint3DHT &lhs, const CPoint3DHT &rhs) {
    return (lhs.x != rhs.x || lhs.y != rhs.y ||
            lhs.z != rhs.z || lhs.w != rhs.w);
  }

  //------

  CPoint3DHT &operator+=(const CPoint3DHT &rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z; w += rhs.w;

    return *this;
  }

  CPoint3DHT operator+(const CPoint3DHT &rhs) const {
    CPoint3DHT t = *this;

    t += rhs;

    return t;
  }

  //------

  CPoint3DHT &operator-=(const CPoint3DHT &rhs) {
    x -= rhs.x; y -= rhs.y; z -= rhs.z; w -= rhs.w;

    return *this;
  }

  CPoint3DHT operator-(const CPoint3DHT &rhs) const {
    CPoint3DHT t = *this;

    t -= rhs;

    return t;
  }

  //------

  CPoint3DHT &operator*=(T rhs) {
    x *= rhs; y *= rhs; z *= rhs; w *= rhs;

    return *this;
  }

  CPoint3DHT operator*(T rhs) const {
    CPoint3DHT t = *this;

    t *= rhs;

    return t;
  }

  //------

  CPoint3DHT &operator/=(T rhs) {
    T irhs = 1.0/rhs;

    x *= irhs; y *= irhs; z *= irhs; w *= irhs;

    return *this;
  }

  CPoint3DHT operator/(T rhs) const {
    CPoint3DHT t = *this;

    t /= rhs;

    return t;
  }

  void print(ostream &os) const {
    return os << "(" << x << "," << y << "," << z << "," << w << ")";
  }

  friend ostream &operator<<(ostream &os, const CPoint3DHT &point) {
    point.print(os);

    return os;
  }
};

typedef CPoint3DHT<double> CPoint3DH;

#endif
