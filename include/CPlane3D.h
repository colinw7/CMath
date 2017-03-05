#ifndef CPLANE_3D_H
#define CPLANE_3D_H

#include <CMathGen.h>
#include <CPoint3D.h>
#include <CLine3D.h>
#include <CVector3D.h>

class CPlane3D {
 public:
  // constructor/destructor
  CPlane3D() { }

  CPlane3D(const CPoint3D &point, const CVector3D &normal) :
   point_(point), normal_(normal) {
    c_ = normal_.dotProduct(CVector3D(point_));
  }

  CPlane3D(const CVector3D &normal, double c) :
   normal_(normal), c_(c) {
    point_ = (normal_*c).point();
  }

  CPlane3D(const CPoint3D &point1, const CPoint3D &point2, const CPoint3D &point3) {
    CVector3D v23(point3, point2); // point2 - point3
    CVector3D v31(point1, point3); // point3 - point1
    CVector3D v12(point2, point1); // point1 - point2

    double A = point1.y*v23.getZ() + point2.y*v31.getZ() + point3.y*v12.getZ();
    double B = point1.z*v23.getX() + point2.z*v31.getX() + point3.z*v12.getX();
    double C = point1.x*v23.getY() + point2.x*v31.getY() + point3.x*v12.getY();

    normal_ = CVector3D(A, B, C);
    point_  = point1;

    c_ = normal_.dotProduct(CVector3D(point_));
  }

 ~CPlane3D() { }

  //------

  // copy operations
  CPlane3D(const CPlane3D &rhs) :
   point_(rhs.point_), normal_(rhs.normal_), c_(rhs.c_) {
  }

  CPlane3D &operator=(const CPlane3D &rhs) {
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

  friend std::ostream &operator<<(std::ostream &os, const CPlane3D &plane) {
    plane.print(os);

    return os;
  }

  //------

  // accessors
  const CPoint3D  &getPoint   () const { return point_ ; }
  const CVector3D &getNormal  () const { return normal_; }
  double             getConstant() const { return c_     ; }

  double value(const CPoint3D &point) const {
    return (normal_.dotProduct(point - point_));
  }

  //------

  // intersect
  CMathGen::IntersectType
  intersectLine(const CLine3D &line, double *iparam) const {
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

  bool intersect(const CLine3D &line, double *iparam) const {
    const CVector3D &vector = line.vector();

    double dot_product1 = vector.dotProduct(normal_);

    if (fabs(dot_product1) < 1E-6)
      return false;

    double dot_product2 = CVector3D::dotProduct(line.start(), normal_);
    double dot_product3 = CVector3D::dotProduct(point_      , normal_);

    *iparam = (dot_product3 - dot_product2)/dot_product1;

    return true;
  }

  bool intersect(const CPlane3D &plane, CLine3D &iline) const {
    const CVector3D &normal1 =       getNormal();
    const CVector3D &normal2 = plane.getNormal();

    double n11 = normal1.dotProduct(normal1);
    double n12 = normal1.dotProduct(normal2);
    double n22 = normal2.dotProduct(normal2);

    double det = n11*n22 - n12*n12;

    if (fabs(det) < 1E-6)
      return false;

    double idet = 1.0/det;

    double c1 = (      getConstant()*n22 - plane.getConstant())*n12*idet;
    double c2 = (plane.getConstant()*n11 -       getConstant())*n12*idet;

    CVector3D p = normal1.crossProduct(normal2);

    CPoint3D p1 = (c1*normal1 + c2*normal2).point();
    CPoint3D p2 = p1 + p;

    iline = CLine3D(p1, p2);

    return true;
  }

  bool intersect(const CLine3D &line, CPoint3D &ipoint, double iparam) const;

 private:
  CPoint3D  point_;
  CVector3D normal_;
  double    c_;
};

//------

#include <CMathGeom3D.h>

inline bool CPlane3D::intersect(const CLine3D &line, CPoint3D &ipoint, double iparam) const {
  return CMathGeom3D::LinePlaneIntersect(line, *this, ipoint, iparam);
}

#endif
