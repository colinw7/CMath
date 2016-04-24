#ifndef CDISK_3D_H
#define CDISK_3D_H

#include <CShape3D.h>
#include <CVector3D.h>
#include <CLine3D.h>
#include <CPoint3D.h>
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

template<typename T>
class CDisk3DT : public CShape3DT<T> {
 private:
  typedef typename CShape3DT<T>::BBox BBox;
  typedef CVector3DT<T>               Vector;
  typedef CLine3DT<T>                 Line;
  typedef CPoint3DT<T>                Point;

  T        radius_      ; //! Radius r
  T        height_      ; //! Height h
  T        inner_radius_; //! Inner Radius ri

  // limits
  T        phi_max_; //! Angle/Sweep Max

  COptReal area_;

 public:
  CDisk3DT(T radius=T(1), T height=T(0), T inner_radius=T(0)) :
   radius_(radius), height_(height), inner_radius_(inner_radius) {
    setPhiLimit(360.0);
  }

  T getRadius() const { return radius_; }

  void setRadius(T radius) {
    radius_ = radius;

    area_.setInvalid();
  }

  T getHeight() const { return height_; }

  void setHeight(T height) {
    height_ = height;

    area_.setInvalid();
  }

  T getInnerRadius() const { return inner_radius_; }

  void setInnerRadius(T radius) {
    inner_radius_ = radius;

    area_.setInvalid();
  }

  void setPhiLimit(T phi_max) {
    phi_max_ = CMathGen::DegToRad(std::min(std::max(phi_max, 0.0), 360.0));

    area_.setInvalid();
  }

  BBox getBBox() const {
    Point p1(-radius_, -radius_, height_);
    Point p2( radius_,  radius_, height_);

    return BBox(CShape3D::transformFrom(p1), CShape3D::transformFrom(p2));
  }

  bool intersect(const Line &line, T *t) const {
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

    Point p1 = CShape3D::transformTo(line.start());
    Point p2 = CShape3D::transformTo(line.end  ());

    Line l(p1, p2);

    Vector ro(l.start());

    const Vector &rd = l.vector();

    if (fabs(rd.getZ()) < 1E-7)
      return false;

    *t = (height_ - ro.getZ())/rd.getZ();

    if (! checkPoint(l.point(*t)))
      return false;

    return true;
  }

  Vector pointNormal(const Point &point) const {
    Point p = CShape3D::transformTo(point);

    Vector dpdu, dpdv;

    pointDetails(p, NULL, NULL, &dpdu, &dpdv);

    dpdu = CShape3D::transformFrom(dpdu);
    dpdv = CShape3D::transformFrom(dpdv);

    Vector n = dpdu.crossProduct(dpdv);

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
      T area = phi_max_*0.5*(radius_*radius_ - inner_radius_*inner_radius_);

      CDisk3DT *th = const_cast<CDisk3DT *>(this);

      th->area_.setValue(area);
    }

    return area_.getValue();
  }

 private:
  bool checkPoint(const Point &point) const {
    T dist2 = point.x*point.x + point.y*point.y;

    if (dist2 > radius_*radius_ || dist2 < inner_radius_*inner_radius_)
      return false;

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

    T dist2 = point.x*point.x + point.y*point.y;

    T tv = 1.0 - ((sqrt(dist2) - inner_radius_)/(radius_ - inner_radius_));

    if (v)
      *v = tv;

    if (dpdu) {
      *dpdu = Vector(-phi_max_*point.y, phi_max_*point.x, 0.0);

      *dpdu *= phi_max_/(2.0*M_PI);
    }

    if (dpdv) {
      T v1  = 1.0 - tv;
      T iv1 = 1.0/v1;

      *dpdv = Vector(-point.x*iv1, -point.y*iv1, 0.0);

      *dpdv *= (radius_ - inner_radius_)/radius_;
    }
  }

  T getPhi(const Point &point) const {
    T phi = atan2(point.y, point.x);

    if (phi < 0.0) phi += 2.0*M_PI;

    return phi;
  }
};

typedef CDisk3DT<double> CDisk3D;

#endif
