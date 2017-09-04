#ifndef CCONE_3D_H
#define CCONE_3D_H

#include <CShape3D.h>
#include <CVector3D.h>
#include <CLine3D.h>
#include <CPoint3D.h>
#include <CMathGen.h>
#include <COptVal.h>

/*! Cone of specified center (cx, cy, cz), radius (r) and height (h)
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

class CCone3D : public CShape3D {
 public:
  CCone3D(double radius=double(1), double height=double(1)) :
   radius_(radius), height_(height) {
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

  void setPhiLimit(double phi_max) {
    phi_max_ = CMathGen::DegToRad(std::min(std::max(phi_max, 0.0), 360.0));

    area_.setInvalid();
  }

  CBBox3D getBBox() const {
    CPoint3D p1(-radius_, -radius_,       0);
    CPoint3D p2( radius_,  radius_, height_);

    return CBBox3D(CShape3D::transformFrom(p1), CShape3D::transformFrom(p2));
  }

  bool intersect(const CLine3D &line, double *tmin, double *tmax) const {
    // solve
    //  k = (r/h)^2
    //
    //  x^2 + y^2 - k*(z - h)^2 = 0
    //
    // at
    //  p = ro + t*rd
    //
    // i.e.
    //  x = ro.x + t*rd.x
    //  y = ro.y + t*rd.y
    //  z = ro.z + t*rd.z

    CPoint3D p1 = CShape3D::transformTo(line.start());
    CPoint3D p2 = CShape3D::transformTo(line.end  ());

    CLine3D l(p1, p2);

    CVector3D ro(l.start());

    const CVector3D &rd = l.vector();

    double rox = ro.getX(); double roy = ro.getY(); double roz = ro.getZ();
    double rdx = rd.getX(); double rdy = rd.getY(); double rdz = rd.getZ();

    double rozh = roz - height_;

    double k = radius_/height_;

    k = k*k;

    double a = rdx*rdx + rdy*rdy - k*rdz*rdz;
    double b = 2.0*(rdx*rox + rdy*roy - k*rdz*rozh);
    double c = rox*rox + roy*roy - k*rozh*rozh;

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
      double h2 = height_*height_;
      double r2 = radius_*radius_;

      double area = phi_max_*h2*sqrt(h2 + r2)/(2.0*radius_);

      CCone3D *th = const_cast<CCone3D *>(this);

      th->area_.setValue(area);
    }

    return area_.getValue();
  }

 private:
  bool checkPoint(const CPoint3D &point) const {
    double z = point.z;

    if (z < 0 || z > height_) return false;

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

    double tv = point.z/height_;

    if (v)
      *v = tv;

    if (dpdu)
      *dpdu = CVector3D(-phi_max_*point.y, phi_max_*point.x, 0.0);

    if (dpdv) {
      double v1  = 1.0 - tv;
      double iv1 = 1.0/v1;

      *dpdv = CVector3D(-point.x*iv1, -point.y*iv1, height_);
    }
  }

  double getPhi(const CPoint3D &point) const {
    double phi = atan2(point.y, point.x);

    if (phi < 0.0) phi += 2.0*M_PI;

    return phi;
  }

 private:
  double   radius_; //! Radius r
  double   height_; //! Height h

  // limits
  double   phi_max_; //! Angle/Sweep Max

  COptReal area_;
};

#endif
