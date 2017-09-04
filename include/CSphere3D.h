#ifndef CSPHERE_3D_H
#define CSPHERE_3D_H

#include <CShape3D.h>
#include <CVector3D.h>
#include <CLine3D.h>
#include <CPoint3D.h>
#include <CVector2D.h>
#include <CMathGen.h>
#include <COptVal.h>

/*! Sphere of specified radius centered at origin
 *
 *  phi   = u*phi_max
 *  theta = theta_min + v(theta_max - theta_min)
 *
 *  x*x + y*y + z*z - r*r = 0
 *
 *  x = r*sin(theta)*cos(phi)
 *  y = r*sin(theta)*sin(phi)
 *  z = r*cos(theta)
*/

class CSphere3D : public CShape3D {
 public:
  CSphere3D(double radius=double(1)) :
   radius_(radius), zmin_(), zmax_(), phi_max_(), theta_min_(), theta_max_(), area_() {
    setZLimit(-radius_, radius_);

    setPhiLimit(360.0);
  }

  double getRadius() const { return radius_; }

  void setRadius(double radius) {
    radius_ = radius;

    area_.setInvalid();
  }

  void setZLimit(double zmin, double zmax) {
    zmin_ = std::min(std::max(std::min(zmin, zmax), -radius_), radius_);
    zmax_ = std::min(std::max(std::max(zmin, zmax), -radius_), radius_);

    theta_min_ = acos(std::min(std::max(zmin_/radius_, -1.0), 1.0));
    theta_max_ = acos(std::min(std::max(zmax_/radius_, -1.0), 1.0));

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
    // solve
    //  x^2 + y^2 + z^2 = r^2 at line (o + t*d)
    //
    // i.e.
    //  (ox + t*dx)^2 + (oy + t*dy)^2 + (oz + t*dz)^2 = r^2
    //
    // solving for t
    //
    //  ox^2 + 2*t*dx*ox + t^2*dx^2 +
    //  oy^2 + 2*t*dy*oy + t^2*dy^2 +
    //  oz^2 + 2*t*dz*oz + t^2*dz^2 = r^2
    //
    // grouping by t
    //
    //  (dx^2 + dy^2 + dz^2)*t^2 + 2*(dx*ox + dy*oy + dz*oz)*t +
    //  (ox^2 + oy^2 + oz^2 - r^2) = 0
    //
    // solve quadratic
    //
    //  at^2 + bt + c = 0
    //
    // where
    //
    //  a = dx^2 + dy^2 + dz^2
    //  b = 2*(dx*ox + dy*oy + dz*oz)
    //  c = ox^2 + oy^2 + oz^2 - r^2
    //

    CPoint3D p1 = CShape3D::transformTo(line.start());
    CPoint3D p2 = CShape3D::transformTo(line.end  ());

    CLine3D l(p1, p2);

    CVector3D ro(l.start());

    const CVector3D &rd = l.vector();

    double a = rd.dotProductSelf();
    double b = 2.0*rd.dotProduct(ro);
    double c = ro.dotProductSelf() - radius_*radius_;

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
      double area = phi_max_*radius_*(zmax_ - zmin_);

      CSphere3D *th = const_cast<CSphere3D *>(this);

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

    double theta     = acos(point.z/radius_);
    double theta_len = theta_max_ - theta_min_;

    if (v)
      *v = (theta - theta_min_)/theta_len;

    if (dpdu || dpdv) {
      double zradius2 = point.x*point.x + point.y*point.y;

      if (zradius2 != 0.0) {
        if (dpdu)
          *dpdu = CVector3D(-phi_max_*point.y, phi_max_*point.x, 0.0);

        if (dpdv) {
          double zradius  = sqrt(zradius2);
          double izradius = 1.0/zradius;

          double cos_phi = point.x*izradius;
          double sin_phi = point.y*izradius;

          *dpdv = theta_len*CVector3D( point.z*cos_phi,
                                    point.z*sin_phi,
                                   -radius_*sin(theta));
        }
      }
      else {
        double cos_phi = 0;
        double sin_phi = 1;

        CVector3D tdpdv = theta_len*CVector3D( point.z*cos_phi,
                                         point.z*sin_phi,
                                        -radius_*sin(theta));

        if (dpdv)
          *dpdv = tdpdv;

        if (dpdu)
          *dpdu = tdpdv.crossProduct(CVector3D(point));
      }
    }
  }

  double getPhi(const CPoint3D &point) const {
    double phi = atan2(point.y, point.x);

    if (phi < 0.0) phi += 2.0*M_PI;

    return phi;
  }

 private:
  double   radius_;      //! Radius R
  double   zmin_, zmax_; //! Height Limits
  double   phi_max_;     //! Angle/Sweep Limit

  double   theta_min_, theta_max_; // Temporaries

  COptReal area_;
};

#endif
