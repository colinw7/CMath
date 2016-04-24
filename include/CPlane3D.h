#ifndef CPLANE_3D_H
#define CPLANE_3D_H

#include <CMathGen.h>
#include <CPoint3D.h>
#include <CLine3D.h>
#include <CVector3D.h>

template<typename T>
class CPlane3DT {
 private:
  typedef CPoint3DT<T>  Point;
  typedef CVector3DT<T> Vector;
  typedef CPlane3DT<T>  Plane;
  typedef CLine3DT<T>   Line;

 private:
  Point  point_;
  Vector normal_;
  T      c_;

 public:
  // constructor/destructor
  CPlane3DT() { }

  CPlane3DT(const Point &point, const Vector &normal) :
   point_(point), normal_(normal) {
    c_ = normal_.dotProduct(CVector3D(point_));
  }

  CPlane3DT(const Vector &normal, T c) :
   normal_(normal), c_(c) {
    point_ = (normal_*c).point();
  }

  CPlane3DT(const Point &point1, const Point &point2, const Point &point3) {
    Vector v23(point3, point2); // point2 - point3
    Vector v31(point1, point3); // point3 - point1
    Vector v12(point2, point1); // point1 - point2

    T A = point1.y*v23.getZ() + point2.y*v31.getZ() + point3.y*v12.getZ();
    T B = point1.z*v23.getX() + point2.z*v31.getX() + point3.z*v12.getX();
    T C = point1.x*v23.getY() + point2.x*v31.getY() + point3.x*v12.getY();

    normal_ = Vector(A, B, C);
    point_  = point1;

    c_ = normal_.dotProduct(CVector3D(point_));
  }

 ~CPlane3DT() { }

  //------

  // copy operations
  CPlane3DT(const Plane &rhs) :
   point_(rhs.point_), normal_(rhs.normal_), c_(rhs.c_) {
  }

  Plane &operator=(const Plane &rhs) {
    point_  = rhs.point_;
    normal_ = rhs.normal_;
    c_      = rhs.c_;

    return *this;
  }

  //------

  // output
  void print(std::ostream &os) const {
    os << "(" << point_ << ") (" << normal_ << ") (" << c_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const Plane &plane) {
    plane.print(os);

    return os;
  }

  //------

  // accessors
  const Point  &getPoint   () const { return point_ ; }
  const Vector &getNormal  () const { return normal_; }
  T             getConstant() const { return c_     ; }

  T value(const Point &point) const {
    return (normal_.dotProduct(point - point_));
  }

  //------

  // intersect
  CMathGen::IntersectType
  intersectLine(const CLine3DT<T> &line, T *iparam) const {
    if (! intersect(line, iparam)) {
      if (fabs(value(line.start())) <= 1E-6)
        return CMathGen::INTERSECT_ALL;
      else
        return CMathGen::INTERSECT_NONE;
    }
    else {
      if (*iparam >= 0.0 && *iparam <= 1.0)
        return CMathGen::INTERSECT_INSIDE;
      else
        return CMathGen::INTERSECT_OUTSIDE;
    }
  }

  bool intersect(const CLine3DT<T> &line, T *iparam) const {
    const CVector3D &vector = line.vector();

    T dot_product1 = vector.dotProduct(normal_);

    if (fabs(dot_product1) < 1E-6)
      return false;

    T dot_product2 = CVector3D::dotProduct(line.start(), normal_);
    T dot_product3 = CVector3D::dotProduct(point_      , normal_);

    *iparam = (dot_product3 - dot_product2)/dot_product1;

    return true;
  }

  bool intersect(const Plane &plane, Line &iline) const {
    const CVector3D &normal1 =       getNormal();
    const CVector3D &normal2 = plane.getNormal();

    T n11 = normal1.dotProduct(normal1);
    T n12 = normal1.dotProduct(normal2);
    T n22 = normal2.dotProduct(normal2);

    T det = n11*n22 - n12*n12;

    if (fabs(det) < 1E-6)
      return false;

    T idet = 1.0/det;

    T c1 = (      getConstant()*n22 - plane.getConstant())*n12*idet;
    T c2 = (plane.getConstant()*n11 -       getConstant())*n12*idet;

    CVector3D p = normal1.crossProduct(normal2);

    CPoint3D p1 = (c1*normal1 + c2*normal2).point();
    CPoint3D p2 = p1 + p;

    iline = CLine3D(p1, p2);

    return true;
  }

  bool intersect(const Line &line, Point &ipoint, T iparam) const;
};

typedef CPlane3DT<double> CPlane3D;

//------

#include <CMathGeom3D.h>

template<typename T>
bool
CPlane3DT<T>::
intersect(const Line &line, Point &ipoint, T iparam) const
{
  return CMathGeom3D::LinePlaneIntersect(line, *this, ipoint, iparam);
}

#endif
