#ifndef CNPLANE_3D_H
#define CNPLANE_3D_H

#include <CPoint3D.h>
#include <CVector3D.h>
#include <CLine3D.h>

template<typename T>
class CNPlane3DT {
 private:
  typedef CNPlane3DT<T> Plane;
  typedef CVector3DT<T> Vector;
  typedef CLine3DT<T>   Line;
  typedef CPoint3DT<T>  Point;

  Vector direction_;
  T      scalar_;

 public:
  CNPlane3DT() :
   direction_(0.0, 0.0, 0.0),
   scalar_   (0.0) {
  }

  CNPlane3DT(const Vector &direction, T scalar=0.0) :
   direction_(direction),
   scalar_   (scalar) {
  }

  CNPlane3DT(const Point &direction, T scalar=0.0) :
   direction_(direction),
   scalar_   (scalar) {
  }

  CNPlane3DT(T x, T y, T z, T scalar=0.0) :
   direction_(x, y, z),
   scalar_   (scalar) {
  }

  CNPlane3DT(const Plane &normal) :
   direction_(normal.direction_),
   scalar_   (normal.scalar_) {
  }

  CNPlane3DT(const Point &v1, const Point &v2, const Point &v3) :
   direction_(), scalar_(0.0) {
    Vector diff1(v1, v2);
    Vector diff2(v2, v3);

    direction_ = CVector3D::crossProduct(diff1, diff2);

    scalar_ = CVector3D::dotProduct(direction_, Vector(v1));
  }

  const Vector &getDirection() const { return direction_; }

  T getDirectionX() const { return direction_.getX(); }
  T getDirectionY() const { return direction_.getY(); }
  T getDirectionZ() const { return direction_.getZ(); }

  T getScalar() const { return scalar_; }

  void setDirection(const Vector &direction) {
    direction_ = direction;
  }

  void setDirection(T x, T y, T z) {
    direction_.setXYZ(x, y, z);
  }

  void setScalar(T scalar) { scalar_ = scalar; }

  void zero() {
    direction_.zero();

    scalar_ = 0.0;
  }

  T modulus() const {
    return direction_.modulus();
  }

  T value(const Point &point) const {
    return (direction_.dotProduct(point) + scalar_);
  }

  void normalize() {
    direction_.normalize();
  }

  Plane &operator+=(const Plane &rhs) {
    direction_ += rhs.direction_;
    scalar_    += rhs.scalar_;

    return *this;
  }

  Plane operator+(const Plane &rhs) const {
    Plane t = *this;

    t += rhs;

    return t;
  }

  Plane &operator-=(const Plane &rhs) {
    direction_ -= rhs.direction_;
    scalar_    -= rhs.scalar_;

    return *this;
  }

  Plane operator-(const Plane &rhs) const {
    Plane t = *this;

    t -= rhs;

    return t;
  }

  Plane &operator*=(T rhs) {
    direction_ *= rhs;

    return *this;
  }

  Plane operator*(T rhs) const {
    Plane t = *this;

    t *= rhs;

    return t;
  }

  Plane &operator/=(T rhs) {
    direction_ /= rhs;

    return *this;
  }

  Plane operator/(T rhs) const {
    Plane t = *this;

    t /= rhs;

    return t;
  }

  bool intersect(const Line &line, Point &ipoint, double iparam) const;

  void print(std::ostream &os) const {
    os << direction_ << " " << scalar_;
  }

  friend Plane operator*(T lhs, const Plane &rhs) {
    return rhs*lhs;
  }

  friend std::ostream &operator<<(std::ostream &os, const Plane &normal) {
    normal.print(os);

    return os;
  }
};

typedef CNPlane3DT<double> CNPlane3D;

//------

#include <CMathGeom3D.h>

template<typename T>
bool
CNPlane3DT<T>::
intersect(const Line &line, Point &ipoint, double iparam) const
{
  return CMathGeom3D::LinePlaneIntersect(line.start(), line.end(), *this,
                                         ipoint, iparam);
}

#endif
