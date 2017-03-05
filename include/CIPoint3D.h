#ifndef CIPOINT_3D_H
#define CIPOINT_3D_H

class CIPoint3D {
 public:
  CIPoint3D(int x1=0, int y1=0, int z1=0) :
   x(x1), y(y1), z(z1) {
  }

  CIPoint3D(const CIPoint3D &point) :
   x(point.x), y(point.y), z(point.z) {
  }

  const CIPoint3D &operator=(const CIPoint3D &point) {
    x = point.x;
    y = point.y;

    return *this;
  }

  int getX() const { return x; }
  int getY() const { return y; }
  int getZ() const { return z; }

  void setX(int x1) { x = x1; }
  void setY(int y1) { y = y1; }
  void setZ(int z1) { z = z1; }

  void setXYZ(int x1, int y1, int z1) {
    x = x1; y = y1; z = z1;
  }

  void getXYZ(int *x1, int *y1, int *z1) const {
    *x1 = x; *y1 = y; *z1 = z;
  }

  CIPoint3D &zero() {
    x = 0; y = 0; z = 0;

    return *this;
  }

  double modulus() const {
    return sqrt(x*x + y*y + z*z);
  }

  CIPoint3D &operator+=(const CIPoint3D &rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;

    return *this;
  }

  CIPoint3D operator+(const CIPoint3D &rhs) const {
    CIPoint3D t = *this;

    t += rhs;

    return t;
  }

  CIPoint3D &operator-=(const CIPoint3D &rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;

    return *this;
  }

  CIPoint3D operator-(const CIPoint3D &rhs) const {
    CIPoint3D t = *this;

    t -= rhs;

    return t;
  }

  CIPoint3D &operator*=(int rhs) {
    x *= rhs;
    y *= rhs;
    z *= rhs;

    return *this;
  }

  CIPoint3D operator*(int rhs) const {
    CIPoint3D t = *this;

    t *= rhs;

    return t;
  }

  CIPoint3D &operator/=(int rhs) {
    x /= rhs;
    y /= rhs;
    z /= rhs;

    return *this;
  }

  CIPoint3D operator/(int rhs) const {
    CIPoint3D t = *this;

    t /= rhs;

    return t;
  }

  friend bool operator==(const CIPoint3D &lhs, const CIPoint3D &rhs) {
    return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z);
  }

  friend bool operator!=(const CIPoint3D &lhs, const CIPoint3D &rhs) {
    return ! (lhs == rhs);
  }

  //-----

  int minComponent() {
    return std::min(std::min(x, y), z);
  }

  int maxComponent() {
    return std::max(std::max(x, y), z);
  }

  //-----

  static CIPoint3D min(const CIPoint3D &lhs, const CIPoint3D &rhs) {
    return CIPoint3D(std::min(lhs.x, rhs.x),
                      std::min(lhs.y, rhs.y),
                      std::min(lhs.z, rhs.z));
  }

  static CIPoint3D max(const CIPoint3D &lhs, const CIPoint3D &rhs) {
    return CIPoint3D(std::max(lhs.x, rhs.x),
                      std::max(lhs.y, rhs.y),
                      std::max(lhs.z, rhs.z));
  }

  //-----

  void print(std::ostream &os) const {
    os << "(" << x << "," << y << "," << z << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CIPoint3D &point) {
    point.print(os);

    return os;
  }

 public:
  int x { 0 }, y { 0 }, z { 0 };
};

#endif
