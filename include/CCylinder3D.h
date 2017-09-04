#ifndef CCYLINDER_3D_H
#define CCYLINDER_3D_H

#include <CShape3D.h>
#include <CVector3D.h>
#include <CLine3D.h>
#include <CPoint3D.h>
#include <CMathGen.h>
#include <COptVal.h>

/*! Cylinder of specified radius (r) and specified height (h)
 *
 *   x*x + y*y - r*r = 0
 *
 *  For (u,v)
 *
 *   phi = 2*u*PI
 *
 *   x = r*cos(phi)
 *   y = y*sin(phi)
 *   z = v*height_;
 */

class CCylinder3D : public CShape3D {
 public:
  CCylinder3D(double radius=double(1), double height=double(1)) :
   radius_(radius), height_(height) {
    setZLimit(0, height_);

    setPhiLimit(360.0);
  }

  double getRadius() const { return radius_; }

  void setRadius(double radius) {
    radius_ = radius;

    area_.setInvalid();
  }

  double getHeight() const { return height_; }

  void setHeight(double height) {
    height_ = height;

    area_.setInvalid();
  }

  void setZLimit(double zmin, double zmax) {
    zmin_ = std::min(std::max(std::min(zmin, zmax), 0.0), height_);
    zmax_ = std::min(std::max(std::max(zmin, zmax), 0.0), height_);

    area_.setInvalid();
  }

  void setPhiLimit(double phi_max) {
    phi_max_ = CMathGen::DegToRad(std::min(std::max(phi_max, 0.0), 360.0));

    area_.setInvalid();
  }

  CBBox3D getBBox() const {
    CPoint3D p1(-radius_, -radius_, zmin_);
    CPoint3D p2( radius_,  radius_, zmax_);

    return CBBox3D(CShape3D::transformFrom(p1), CShape3D::transformFrom(p2));
  }

  bool intersect(const CLine3D &line, double *tmin, double *tmax) const {
    // Solve
    //
    //  x^2 + y^2 - r^2 = 0 at ray (o + t*d)
    //
    // solving for t
    //
    //  ox^2 + 2*t*dx*ox + t^2*dx^2 +
    //  oy^2 + 2*t*dy*oy + t^2*dy^2 = r^2
    //
    // grouping by t
    //
    //  (dx^2 + dy^2)*t^2 + 2*(dx*ox + dy*oy)*t +
    //  (ox^2 + oy^2 - r^2) = 0
    //
    // solve quadratic
    //
    //  at^2 + bt + c = 0
    //
    // where
    //
    //  a = dx^2 + dy^2
    //  b = 2*(dx*ox + dy*oy)
    //  c = ox^2 + oy^2 - r^2
    //

    CPoint3D p1 = CShape3D::transformTo(line.start());
    CPoint3D p2 = CShape3D::transformTo(line.end  ());

    CLine3D l(p1, p2);

    CVector3D ro(l.start());

    const CVector3D &rd = l.vector();

    double rox = ro.getX(); double roy = ro.getY();
    double rdx = rd.getX(); double rdy = rd.getY();

    double a = rdx*rdx + rdy*rdy;
    double b = 2.0*(rdx*rox + rdy*roy);
    double c = rox*rox + roy*roy - radius_*radius_;

    if (! CMathGen::solveQuadratic(a, b, c, tmin, tmax))
      return false;

    bool b1 = checkPoint(l.point(*tmin));
    bool b2 = checkPoint(l.point(*tmax));

    if (! b1 && ! b2) return false;

    if      (! b1) *tmin = *tmax;
    else if (! b2) *tmax = *tmin;

    return true;
  }

  CVector3D pointNormal(const CPoint3D &point) const {
    CPoint3D p = CShape3D::transformTo(point);

    CVector3D dpdu, dpdv;

    pointDetails(p, nullptr, nullptr, &dpdu, &dpdv);

    dpdu = CShape3D::transformFrom(dpdu);
    dpdv = CShape3D::transformFrom(dpdv);

    CVector3D n = dpdu.crossProduct(dpdv).unit();

    return n;
  }

  CVector2D pointToSurfaceVector(const CPoint3D &point) const {
    CPoint3D p = CShape3D::transformTo(point);

    double u, v;

    pointDetails(p, &u, &v);

    return CVector2D(u, v);
  }

  double getArea() const {
    if (! area_.isValid()) {
      double area = (zmax_ - zmin_)*phi_max_*radius_;

      CCylinder3D *th = const_cast<CCylinder3D *>(this);

      th->area_.setValue(area);
    }

    return area_.getValue();
  }

 private:
  bool checkPoint(const CPoint3D &point) const {
    double z = point.z;

    if (z < zmin_ || z > zmax_) return false;

    double phi = getPhi(point);

    if (phi > phi_max_)
      return false;

    return true;
  }

  void pointDetails(const CPoint3D &point, double *u, double *v,
                    CVector3D *dpdu=0, CVector3D *dpdv=0) const {
    if (u) {
      double phi = getPhi(point);

      *u = phi/phi_max_;
    }

    if (v)
      *v = (point.z - zmin_)/(zmax_ - zmin_);

    if (dpdu)
      *dpdu = CVector3D(-phi_max_*point.y, phi_max_*point.x, 0.0);

    if (dpdv)
      *dpdv = CVector3D(0, 0, zmax_ - zmin_);
  }

  double getPhi(const CPoint3D &point) const {
    double phi = atan2(point.y, point.x);

    if (phi < 0.0) phi += 2.0*M_PI;

    return phi;
  }

 private:
  double   radius_ { 1.0 }; //! Radius r
  double   height_ { 1.0 }; //! Height h

  // limits
  double   zmin_    { 0.0 }, zmax_ { 0.0 }; //! Height Range
  double   phi_max_ { 0.0 };                //! Angle/Sweep Max

  COptReal area_;
};

#endif
