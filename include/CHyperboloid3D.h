#ifndef CHYPERBOLOID_3D_H
#define CHYPERBOLOID_3D_H

#include <CShape3D.h>
#include <CVector3D.h>
#include <CLine3D.h>
#include <CPoint3D.h>
#include <CMathGen.h>
#include <COptVal.h>
#include <NaN.h>

/*! Hyperboloid betweem specified points and sweep phi
 *
 * a*x^2 + a*y^2 - c*c^2 = 1
 *
 * phi = u*phi_max
 *
 * xr = (1 - v)*x1 + v*x2
 * yr = (1 - v)*y1 + v*y2
 *
 * x = xr*cos(phi) - yr*sin(phi)
 * y = xr*sin(phi) + yr*cos(phi)
 * z = (1 - v)*z1 + v*z2
 *
 * x = a*sinh(u)*cos(v)  ??
 * y = a*sinh(u)*sin(v)  ??
 * z = c*cosh(u)         ??
 */

class CHyperboloid3D : public CShape3D {
 public:
  CHyperboloid3D(const CPoint3D &point1=CPoint3D(0,0,0), const CPoint3D &point2=CPoint3D(1,1,1)) {
    setPoints(point1, point2);

    setPhiLimit(360.0);
  }

  void setPoints(const CPoint3D &point1, const CPoint3D &point2) {
    point1_ = point1;
    point2_ = point2;

    double rad1 = sqrt(point1.x*point1.x + point1.y*point1.y);
    double rad2 = sqrt(point2.x*point2.x + point2.y*point2.y);

    rmax_ = std::max(rad1, rad2);
    zmin_ = std::min(point1.z, point2.z);
    zmax_ = std::max(point1.z, point2.z);

    //---

    // Compute implicit function coefficients for hyperboloid
    if (point2_.z == 0.0)
      std::swap(point1_, point2_);

    CPoint3D pp = point1_;

    uint num_iters = 0;

    double xy1, xy2;

    do {
      pp += 2.0*(point2_ - point1_);

      xy1 = pp.x*pp.x + pp.y*pp.y;
      xy2 = point2_.x*point2_.x + point2_.y*point2_.y;

      a_ = (1.0/xy1 - (pp.z*pp.z)/(xy1*point2_.z*point2_.z))/
           (1.0 - (xy2*pp.z*pp.z)/(xy1*point2_.z*point2_.z));

      c_ = (a_*xy2 - 1.0)/(point2_.z*point2_.z);

      ++num_iters;
    } while (num_iters < 1000 && (IsInf(a_) || IsNaN(a_)));

    area_.setInvalid();
  }

  void setPhiLimit(double phi_max) {
    phi_max_ = CMathGen::DegToRad(std::min(std::max(phi_max, 0.0), 360.0));

    area_.setInvalid();
  }

  CBBox3D getBBox() const {
    CPoint3D p1(-rmax_, -rmax_, zmin_);
    CPoint3D p2( rmax_,  rmax_, zmax_);

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

    double a = a_*rdx*rdx + a_*rdy*rdy - c_*rdz*rdz;
    double b = 2.0*(a_*rdx*rox + a_*rdy*roy - c_*rdz*roz);
    double c = a_*rox*rox + a_*roy*roy - c_*roz*roz - 1.0;

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
      double p1x2 = point1_.x*point1_.x;
      double p1y2 = point1_.y*point1_.y;
      double p1z2 = point1_.z*point1_.z;

      double p2x2 = point2_.x*point2_.x;
      double p2y2 = point2_.y*point2_.y;
      double p2z2 = point2_.z*point2_.z;

      double dp12y = point1_.y - point2_.y;
      double dp12z = point1_.z - point2_.z;

      double p12x = point1_.x*point2_.x;
      double p12y = point1_.y*point2_.y;
      double p12z = point1_.z*point2_.z;

      double p1x4 = p1x2*p1x2;
      double p2x4 = p2x2*p2x2;

      double dp12y2 = dp12y*dp12y;
      double dp12z2 = dp12z*dp12z;

      double area = phi_max_/6.0*
        (2.0*p1x4 - 2.0*p1x2*p12x +
         2.0*p2x4 + 2.0*(p1y2 + p12y + p2y2)*(dp12y2 + dp12z2) +
         p2x2*( 5.0*p1y2 + 2.0*p12y - 4.0*p2y2 + 2.0*dp12z2) +
         p1x2*(-4.0*p1y2 + 2.0*p12y + 5.0*p2y2 + 2.0*dp12z2) -
         2.0*p12x*(p2x2 - p1y2 + 5.0*p12y - p2y2 - p1z2 + 2.0*p12z - p2z2));

      CHyperboloid3D *th = const_cast<CHyperboloid3D *>(this);

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
    double phi = getPhi(point);

    if (u)
      *u = phi/phi_max_;

    if (v)
      *v = (point.z - point1_.z)/(point2_.z - point1_.z);

    if (dpdu)
      *dpdu = CVector3D(-phi_max_*point.y, phi_max_*point.x, 0.0);

    if (dpdv) {
      double cosphi = cos(phi);
      double sinphi = sin(phi);

      *dpdv = CVector3D((point2_.x - point1_.x)*cosphi -
                     (point2_.y - point1_.y)*sinphi,
                     (point2_.x - point1_.x)*sinphi +
                     (point2_.y - point1_.y)*cosphi,
                     point2_.z - point1_.z);
    }
  }

  double getPhi(const CPoint3D &point) const {
    double v = (point.z - point1_.z)/(point2_.z - point1_.z);

    CPoint3D pr = (1.0 - v)*point1_ + v*point2_;

    double phi = atan2(pr.x*point.y - pr.y*point.x, pr.x*point.x + pr.y*point.y);

    if (phi < 0.0) phi += 2.0*M_PI;

    return phi;
  }

 private:
  CPoint3D point1_;                         //! CPoint3D 1
  CPoint3D point2_;                         //! CPoint3D 2
  double   phi_max_ { 0.0 };                //! Angle/Sweep Max
  double   zmin_    { 0.0 }, zmax_ { 0.0 }; //! Height Range
  double   rmax_    { 0.0 };                //!
  double   a_       { 0.0 }, c_ { 0.0 };
  COptReal area_;
};

#endif
