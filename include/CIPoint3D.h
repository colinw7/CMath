#ifndef CIPOINT_3D_H
#define CIPOINT_3D_H

template<typename T>
class CIPoint3DT {
 public:
  T x, y, z;

  CIPoint3DT(T x1=0, T y1=0, T z1=0) :
   x(x1), y(y1), z(z1) {
  }

  CIPoint3DT(const CIPoint3DT &point) :
   x(point.x), y(point.y), z(point.z) {
  }

  const CIPoint3DT &operator=(const CIPoint3DT &point) {
    x = point.x;
    y = point.y;

    return *this;
  }

  T getX() const { return x; }
  T getY() const { return y; }
  T getZ() const { return z; }

  void setX(T x1) { x = x1; }
  void setY(T y1) { y = y1; }
  void setZ(T z1) { z = z1; }

  void setXYZ(T x1, T y1, T z1) {
    x = x1; y = y1; z = z1;
  }

  void getXYZ(T *x1, T *y1, T *z1) const {
    *x1 = x; *y1 = y; *z1 = z;
  }

  CIPoint3DT &zero() {
    x = 0; y = 0; z = 0;

    return *this;
  }

  double modulus() const {
    return sqrt(x*x + y*y + z*z);
  }

  CIPoint3DT &operator+=(const CIPoint3DT &rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;

    return *this;
  }

  CIPoint3DT operator+(const CIPoint3DT &rhs) const {
    CIPoint3DT t = *this;

    t += rhs;

    return t;
  }

  CIPoint3DT &operator-=(const CIPoint3DT &rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;

    return *this;
  }

  CIPoint3DT operator-(const CIPoint3DT &rhs) const {
    CIPoint3DT t = *this;

    t -= rhs;

    return t;
  }

  CIPoint3DT &operator*=(T rhs) {
    x *= rhs;
    y *= rhs;
    z *= rhs;

    return *this;
  }

  CIPoint3DT operator*(T rhs) const {
    CIPoint3DT t = *this;

    t *= rhs;

    return t;
  }

  CIPoint3DT &operator/=(T rhs) {
    x /= rhs;
    y /= rhs;
    z /= rhs;

    return *this;
  }

  CIPoint3DT operator/(T rhs) const {
    CIPoint3DT t = *this;

    t /= rhs;

    return t;
  }

  friend bool operator==(const CIPoint3DT &lhs, const CIPoint3DT &rhs) {
    return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z);
  }

  friend bool operator!=(const CIPoint3DT &lhs, const CIPoint3DT &rhs) {
    return ! (lhs == rhs);
  }

  //-----

  T minComponent() {
    return std::min(std::min(x, y), z);
  }

  T maxComponent() {
    return std::max(std::max(x, y), z);
  }

  //-----

  static CIPoint3DT min(const CIPoint3DT &lhs, const CIPoint3DT &rhs) {
    return CIPoint3DT(std::min(lhs.x, rhs.x),
                      std::min(lhs.y, rhs.y),
                      std::min(lhs.z, rhs.z));
  }

  static CIPoint3DT max(const CIPoint3DT &lhs, const CIPoint3DT &rhs) {
    return CIPoint3DT(std::max(lhs.x, rhs.x),
                      std::max(lhs.y, rhs.y),
                      std::max(lhs.z, rhs.z));
  }

  //-----

  void print(std::ostream &os) const {
    os << "(" << x << "," << y << "," << z << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CIPoint3DT &point) {
    point.print(os);

    return os;
  }
};

typedef CIPoint3DT<int> CIPoint3D;

#endif
