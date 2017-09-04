#ifndef CDISK_3D_H
#define CDISK_3D_H

#include <CShape3D.h>
#include <CVector3D.h>
#include <CLine3D.h>
#include <CPoint3D.h>
#include <CMathGen.h>
#include <COptVal.h>

/*! Disk of specified radius (r), specified z (h)
 *  and specified inner radius (ri) centered at origin
 *
 *   x*x + y*y - r*r = 0
 *
 *  For (u,v)
 *
 *   phi = u*PI
 *
 *   x = ((1 - v)*ri + v*r)*cos(phi)
 *   y = ((1 - v)*ri + v*r)*sin(phi)
 *   z = h
 */

class CDisk3D : public CShape3D {
 public:
  CDisk3D(double radius=double(1), double height=double(0), double inner_radius=double(0)) :
   radius_(radius), height_(height), inner_radius_(inner_radius) {
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

  double getInnerRadius() const { return inner_radius_; }

  void setInnerRadius(double radius) {
    inner_radius_ = radius;

    area_.setInvalid();
  }

  void setPhiLimit(double phi_max) {
    phi_max_ = CMathGen::DegToRad(std::min(std::max(phi_max, 0.0), 360.0));

    area_.setInvalid();
  }

  CBBox3D getBBox() const {
    CPoint3D p1(-radius_, -radius_, height_);
    CPoint3D p2( radius_,  radius_, height_);

    return CBBox3D(CShape3D::transformFrom(p1), CShape3D::transformFrom(p2));
  }

  bool intersect(const CLine3D &line, double *t) const {
    // Solve
    //
    //  z = h at line (o + t*d)
    //
    // i.e.
    //  (oz + t*dz) = h
    //
    // solving for t
    //
    //  t = (h - oz)/dz

    CPoint3D p1 = CShape3D::transformTo(line.start());
    CPoint3D p2 = CShape3D::transformTo(line.end  ());

    CLine3D l(p1, p2);

    CVector3D ro(l.start());

    const CVector3D &rd = l.vector();

    if (fabs(rd.getZ()) < 1E-7)
      return false;

    *t = (height_ - ro.getZ())/rd.getZ();

    if (! checkPoint(l.point(*t)))
      return false;

    return true;
  }

  CVector3D pointNormal(const CPoint3D &point) const {
    CPoint3D p = CShape3D::transformTo(point);

    CVector3D dpdu, dpdv;

    pointDetails(p, nullptr, nullptr, &dpdu, &dpdv);

    dpdu = CShape3D::transformFrom(dpdu);
    dpdv = CShape3D::transformFrom(dpdv);

    CVector3D n = dpdu.crossProduct(dpdv);

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
      double area = phi_max_*0.5*(radius_*radius_ - inner_radius_*inner_radius_);

      CDisk3D *th = const_cast<CDisk3D *>(this);

      th->area_.setValue(area);
    }

    return area_.getValue();
  }

 private:
  bool checkPoint(const CPoint3D &point) const {
    double dist2 = point.x*point.x + point.y*point.y;

    if (dist2 > radius_*radius_ || dist2 < inner_radius_*inner_radius_)
      return false;

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

    double dist2 = point.x*point.x + point.y*point.y;

    double tv = 1.0 - ((sqrt(dist2) - inner_radius_)/(radius_ - inner_radius_));

    if (v)
      *v = tv;

    if (dpdu) {
      *dpdu = CVector3D(-phi_max_*point.y, phi_max_*point.x, 0.0);

      *dpdu *= phi_max_/(2.0*M_PI);
    }

    if (dpdv) {
      double v1  = 1.0 - tv;
      double iv1 = 1.0/v1;

      *dpdv = CVector3D(-point.x*iv1, -point.y*iv1, 0.0);

      *dpdv *= (radius_ - inner_radius_)/radius_;
    }
  }

  double getPhi(const CPoint3D &point) const {
    double phi = atan2(point.y, point.x);

    if (phi < 0.0) phi += 2.0*M_PI;

    return phi;
  }

 private:
  double   radius_      ; //! Radius r
  double   height_      ; //! Height h
  double   inner_radius_; //! Inner Radius ri

  // limits
  double   phi_max_; //! Angle/Sweep Max

  COptReal area_;
};

#endif
