#ifndef CPARABOLOID_3D_H
#define CPARABOLOID_3D_H

#include <CShape3D.h>
#include <CVector3D.h>
#include <CLine3D.h>
#include <CPoint3D.h>
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

template<typename T>
class CParaboloid3DT : public CShape3DT<T> {
 private:
  typedef typename CShape3DT<T>::BBox BBox;
  typedef CVector3DT<T>               Vector;
  typedef CLine3DT<T>                 Line;
  typedef CPoint3DT<T>                Point;

  T        radius_;      //! Radius r
  T        zmin_, zmax_; //!

  T        phi_max_;     //! Angle/Sweep Max

  COptReal area_;

 public:
  CParaboloid3DT(T radius=T(1), T zmin=T(0), T zmax=T(1)) :
   radius_(radius) {
    setZRange(zmin, zmax);

    setPhiLimit(360.0);
  }

  T getRadius() const { return radius_; }

  void setRadius(T radius) {
    radius_ = radius;

    area_.setInvalid();
  }

  T getZMin() const { return zmin_; }
  T getZMax() const { return zmax_; }

  void setZRange(T zmin, T zmax) {
    zmin_ = std::min(zmin, zmax);
    zmax_ = std::max(zmin, zmax);

    area_.setInvalid();
  }

  void setPhiLimit(T phi_max) {
    phi_max_ = CMathGen::DegToRad(std::min(std::max(phi_max, 0.0), 360.0));

    area_.setInvalid();
  }

  BBox getBBox() const {
    Point p1(-radius_, -radius_, zmin_);
    Point p2( radius_,  radius_, zmax_);

    return BBox(CShape3D::transformFrom(p1), CShape3D::transformFrom(p2));
  }

  bool intersect(const Line &line, T *tmin, T *tmax) const {
    Point p1 = CShape3D::transformTo(line.start());
    Point p2 = CShape3D::transformTo(line.end  ());

    Line l(p1, p2);

    Vector ro(l.start());

    const Vector &rd = l.vector();

    T rox = ro.getX(); T roy = ro.getY(); T roz = ro.getZ();
    T rdx = rd.getX(); T rdy = rd.getY(); T rdz = rd.getZ();

    T k = zmax_/(radius_*radius_);

    T a = k*(rdx*rdx + rdy*rdy);
    T b = 2.0*k*(rdx*rox + rdy*roy) - rdz;
    T c = k*(rox*rox + roy*roy) - roz;

    if (! CMathGen::solveQuadratic(a, b, c, tmin, tmax))
      return false;

    bool b1 = checkPoint(l.point(*tmin));
    bool b2 = checkPoint(l.point(*tmax));

    if (! b1 && ! b2) return false;

    if      (! b1) *tmin = *tmax;
    else if (! b2) *tmax = *tmin;

    return true;
  }

  Vector pointNormal(const Point &point) const {
    Point p = CShape3D::transformTo(point);

    Vector dpdu, dpdv;

    pointDetails(p, NULL, NULL, &dpdu, &dpdv);

    dpdu = CShape3D::transformFrom(dpdu);
    dpdv = CShape3D::transformFrom(dpdv);

    Vector n = dpdu.crossProduct(dpdv).unit();

    return n;
  }

  CVector2D pointToSurfaceVector(const Point &point) const {
    Point p = CShape3D::transformTo(point);

    T u, v;

    pointDetails(p, &u, &v);

    return CVector2D(u, v);
  }

  T getArea() const {
    if (! area_.isValid()) {
      T area = phi_max_/12.0*(pow(1 + 4*zmin_, 1.5) -
                              pow(1 + 4*zmax_, 1.5));

      CParaboloid3DT *th = const_cast<CParaboloid3DT *>(this);

      th->area_.setValue(area);
    }

    return area_.getValue();
  }

 private:
  bool checkPoint(const Point &point) const {
    T z = point.z;

    if (z < zmin_ || z > zmax_) return false;

    T phi = getPhi(point);

    if (phi > phi_max_)
      return false;

    return true;
  }

  void pointDetails(const Point &point, T *u, T *v,
                    Vector *dpdu=0, Vector *dpdv=0) const {
    if (u) {
      T phi = getPhi(point);

      *u = phi/phi_max_;
    }

    if (v)
      *v = (point.z - zmin_)/(zmax_ - zmin_);

    if (dpdu)
      *dpdu = Vector(-phi_max_*point.y, phi_max_*point.x, 0.0);

    if (dpdv) {
      T ipz = 1.0/(2.0*point.z);

      *dpdv = (zmax_ - zmin_)*Vector(point.x*ipz, point.y*ipz, 1.0);
    }
  }

  T getPhi(const Point &point) const {
    T phi = atan2(point.y, point.x);

    if (phi < 0.0) phi += 2.0*M_PI;

    return phi;
  }
};

typedef CParaboloid3DT<double> CParaboloid3D;

#endif
