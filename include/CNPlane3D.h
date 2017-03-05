#ifndef CNPLANE_3D_H
#define CNPLANE_3D_H

#include <CPoint3D.h>
#include <CVector3D.h>
#include <CLine3D.h>

class CNPlane3D {
 public:
  CNPlane3D() :
   direction_(0.0, 0.0, 0.0),
   scalar_   (0.0) {
  }

  CNPlane3D(const CVector3D &direction, double scalar=0.0) :
   direction_(direction),
   scalar_   (scalar) {
  }

  CNPlane3D(const CPoint3D &direction, double scalar=0.0) :
   direction_(direction),
   scalar_   (scalar) {
  }

  CNPlane3D(double x, double y, double z, double scalar=0.0) :
   direction_(x, y, z),
   scalar_   (scalar) {
  }

  CNPlane3D(const CNPlane3D &normal) :
   direction_(normal.direction_),
   scalar_   (normal.scalar_) {
  }

  CNPlane3D(const CPoint3D &v1, const CPoint3D &v2, const CPoint3D &v3) :
   direction_(), scalar_(0.0) {
    CVector3D diff1(v1, v2);
    CVector3D diff2(v2, v3);

    direction_ = CVector3D::crossProduct(diff1, diff2);

    scalar_ = CVector3D::dotProduct(direction_, CVector3D(v1));
  }

  const CVector3D &getDirection() const { return direction_; }

  double getDirectionX() const { return direction_.getX(); }
  double getDirectionY() const { return direction_.getY(); }
  double getDirectionZ() const { return direction_.getZ(); }

  double getScalar() const { return scalar_; }

  void setDirection(const CVector3D &direction) {
    direction_ = direction;
  }

  void setDirection(double x, double y, double z) {
    direction_.setXYZ(x, y, z);
  }

  void setScalar(double scalar) { scalar_ = scalar; }

  void zero() {
    direction_.zero();

    scalar_ = 0.0;
  }

  double modulus() const {
    return direction_.modulus();
  }

  double value(const CPoint3D &point) const {
    return (direction_.dotProduct(point) + scalar_);
  }

  void normalize() {
    direction_.normalize();
  }

  CNPlane3D &operator+=(const CNPlane3D &rhs) {
    direction_ += rhs.direction_;
    scalar_    += rhs.scalar_;

    return *this;
  }

  CNPlane3D operator+(const CNPlane3D &rhs) const {
    CNPlane3D t = *this;

    t += rhs;

    return t;
  }

  CNPlane3D &operator-=(const CNPlane3D &rhs) {
    direction_ -= rhs.direction_;
    scalar_    -= rhs.scalar_;

    return *this;
  }

  CNPlane3D operator-(const CNPlane3D &rhs) const {
    CNPlane3D t = *this;

    t -= rhs;

    return t;
  }

  CNPlane3D &operator*=(double rhs) {
    direction_ *= rhs;

    return *this;
  }

  CNPlane3D operator*(double rhs) const {
    CNPlane3D t = *this;

    t *= rhs;

    return t;
  }

  CNPlane3D &operator/=(double rhs) {
    direction_ /= rhs;

    return *this;
  }

  CNPlane3D operator/(double rhs) const {
    CNPlane3D t = *this;

    t /= rhs;

    return t;
  }

  bool intersect(const CLine3D &line, CPoint3D &ipoint, double iparam) const;

  void print(std::ostream &os) const {
    os << direction_ << " " << scalar_;
  }

  friend CNPlane3D operator*(double lhs, const CNPlane3D &rhs) {
    return rhs*lhs;
  }

  friend std::ostream &operator<<(std::ostream &os, const CNPlane3D &normal) {
    normal.print(os);

    return os;
  }

 private:
  CVector3D direction_;
  double    scalar_;
};

//------

#include <CMathGeom3D.h>

inline bool CNPlane3D::intersect(const CLine3D &line, CPoint3D &ipoint, double iparam) const {
  return CMathGeom3D::LinePlaneIntersect(line.start(), line.end(), *this, ipoint, iparam);
}

#endif
