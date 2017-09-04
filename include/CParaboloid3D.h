#ifndef CPARABOLOID_3D_H
#define CPARABOLOID_3D_H

#include <CShape3D.h>
#include <CVector3D.h>
#include <CLine3D.h>
#include <CPoint3D.h>
#include <CMathGen.h>
#include <COptVal.h>

/*! Paraboloid of specified center (cx, cy, cz), radius (r) and height (h)
 *
 *   (h*x/r)^2 + (h*y/r)^2 - (z - h)^2 = 0
 *
 *  For (u,v)
 *
 *   phi = 2*u*PI
 *
 *   x = xc + r*(1 - v)*cos(phi)
 *   y = yc + r*(1 - v)*sin(phi)
 *   z = v*h
 */

class CParaboloid3D : public CShape3D {
 public:
  CParaboloid3D(double radius=double(1), double zmin=double(0), double zmax=double(1)) :
   radius_(radius) {
    setZRange(zmin, zmax);

    setPhiLimit(360.0);
  }

  double getRadius() const { return radius_; }

  void setRadius(double radius) {
    radius_ = radius;

    area_.setInvalid();
  }

  double getZMin() const { return zmin_; }
  double getZMax() const { return zmax_; }

  void setZRange(double zmin, double zmax) {
    zmin_ = std::min(zmin, zmax);
    zmax_ = std::max(zmin, zmax);

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
    CPoint3D p1 = CShape3D::transformTo(line.start());
    CPoint3D p2 = CShape3D::transformTo(line.end  ());

    CLine3D l(p1, p2);

    CVector3D ro(l.start());

    const CVector3D &rd = l.vector();

    double rox = ro.getX(); double roy = ro.getY(); double roz = ro.getZ();
    double rdx = rd.getX(); double rdy = rd.getY(); double rdz = rd.getZ();

    double k = zmax_/(radius_*radius_);

    double a = k*(rdx*rdx + rdy*rdy);
    double b = 2.0*k*(rdx*rox + rdy*roy) - rdz;
    double c = k*(rox*rox + roy*roy) - roz;

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
      double area = phi_max_/12.0*(pow(1 + 4*zmin_, 1.5) -
                              pow(1 + 4*zmax_, 1.5));

      CParaboloid3D *th = const_cast<CParaboloid3D *>(this);

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

    if (dpdv) {
      double ipz = 1.0/(2.0*point.z);

      *dpdv = (zmax_ - zmin_)*CVector3D(point.x*ipz, point.y*ipz, 1.0);
    }
  }

  double getPhi(const CPoint3D &point) const {
    double phi = atan2(point.y, point.x);

    if (phi < 0.0) phi += 2.0*M_PI;

    return phi;
  }

 private:
  double   radius_  { 0.0 };                //! Radius r
  double   zmin_    { 0.0 }, zmax_ { 0.0 }; //!
  double   phi_max_ { 0.0 };                //! Angle/Sweep Max
  COptReal area_;
};

#endif
