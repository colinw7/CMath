#ifndef CPOINT_3DH_H
#define CPOINT_3DH_H

// 3D Homogeneous point

#include <CPoint3D.h>

class CPoint3DH {
  typedef CPoint3D Point;

 public:
  CPoint3DH() { }

  CPoint3DH(double x1, double y1, double z1, double w1=1) :
   x(x1), y(y1), z(z1), w(w1) {
  }

  CPoint3DH(const CPoint3DH &point) :
   x(point.x), y(point.y), z(point.z), w(point.w) {
  }

  double getX() const { return x; }
  double getY() const { return y; }
  double getZ() const { return z; }
  double getW() const { return w; }

  void getXYZW(double *x1, double *y1, double *z1, double *w1) const {
    *x1 = x; *y1 = y; *z1 = z; *w1 = w;
  }

  Point getPoint() const {
    return Point(x/w, y/w, z/w);
  }

  double operator[](int i) const { return (&x)[i]; }

  double &operator[](int i) { return (&x)[i]; }

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

  double getDistance(const CPoint3DH &point) const {
    CPoint3DH diff = *this - point;

    return diff.modulus();
  }

  double getDistanceSquared(const CPoint3DH &point) const {
    CPoint3DH diff = *this - point;

    return diff.modulusSquared();
  }

  CPoint3DH &zero() {
    x = 0.0; y = 0.0; z = 0.0; w = 0.0;

    return *this;
  }

  //-----

  CPoint3DH operator+() {
    return CPoint3DH(x, y, z);
  }

  CPoint3DH operator-() {
    return CPoint3DH(-x, -y, -z);
  }

  //-----

  friend bool operator==(const CPoint3DH &lhs, const CPoint3DH &rhs) {
    return (lhs.x == rhs.x && lhs.y == rhs.y &&
            lhs.z == rhs.z && lhs.w == rhs.w);
  }

  friend bool operator!=(const CPoint3DH &lhs, const CPoint3DH &rhs) {
    return (lhs.x != rhs.x || lhs.y != rhs.y ||
            lhs.z != rhs.z || lhs.w != rhs.w);
  }

  //------

  CPoint3DH &operator+=(const CPoint3DH &rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z; w += rhs.w;

    return *this;
  }

  CPoint3DH operator+(const CPoint3DH &rhs) const {
    CPoint3DH t = *this;

    t += rhs;

    return t;
  }

  //------

  CPoint3DH &operator-=(const CPoint3DH &rhs) {
    x -= rhs.x; y -= rhs.y; z -= rhs.z; w -= rhs.w;

    return *this;
  }

  CPoint3DH operator-(const CPoint3DH &rhs) const {
    CPoint3DH t = *this;

    t -= rhs;

    return t;
  }

  //------

  CPoint3DH &operator*=(double rhs) {
    x *= rhs; y *= rhs; z *= rhs; w *= rhs;

    return *this;
  }

  CPoint3DH operator*(double rhs) const {
    CPoint3DH t = *this;

    t *= rhs;

    return t;
  }

  //------

  CPoint3DH &operator/=(double rhs) {
    double irhs = 1.0/rhs;

    x *= irhs; y *= irhs; z *= irhs; w *= irhs;

    return *this;
  }

  CPoint3DH operator/(double rhs) const {
    CPoint3DH t = *this;

    t /= rhs;

    return t;
  }

  void print(std::ostream &os) const {
    os << "(" << x << "," << y << "," << z << "," << w << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CPoint3DH &point) {
    point.print(os);

    return os;
  }

 public:
  double x { 0 }, y { 0 }, z { 0 }, w { 1 };
};

#endif
