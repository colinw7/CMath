#ifndef CCONE_3D_H
#define CCONE_3D_H

#include <CShape3D.h>
#include <CVector3D.h>
#include <CLine3D.h>
#include <CPoint3D.h>
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

template<typename T>
class CCone3DT : public CShape3DT<T> {
 private:
  typedef typename CShape3DT<T>::BBox BBox;
  typedef CVector3DT<T>               Vector;
  typedef CLine3DT<T>                 Line;
  typedef CPoint3DT<T>                Point;

  T        radius_; //! Radius r
  T        height_; //! Height h

  // limits
  T        phi_max_; //! Angle/Sweep Max

  COptReal area_;

 public:
  CCone3DT(T radius=T(1), T height=T(1)) :
   radius_(radius), height_(height) {
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

  void setPhiLimit(T phi_max) {
    phi_max_ = CMathGen::DegToRad(std::min(std::max(phi_max, 0.0), 360.0));

    area_.setInvalid();
  }

  BBox getBBox() const {
    Point p1(-radius_, -radius_,       0);
    Point p2( radius_,  radius_, height_);

    return BBox(CShape3D::transformFrom(p1), CShape3D::transformFrom(p2));
  }

  bool intersect(const Line &line, T *tmin, T *tmax) const {
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

    Point p1 = CShape3D::transformTo(line.start());
    Point p2 = CShape3D::transformTo(line.end  ());

    Line l(p1, p2);

    Vector ro(l.start());

    const Vector &rd = l.vector();

    T rox = ro.getX(); T roy = ro.getY(); T roz = ro.getZ();
    T rdx = rd.getX(); T rdy = rd.getY(); T rdz = rd.getZ();

    T rozh = roz - height_;

    T k = radius_/height_;

    k = k*k;

    T a = rdx*rdx + rdy*rdy - k*rdz*rdz;
    T b = 2.0*(rdx*rox + rdy*roy - k*rdz*rozh);
    T c = rox*rox + roy*roy - k*rozh*rozh;

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
    if (! area_.getValid()) {
      T h2 = height_*height_;
      T r2 = radius_*radius_;

      T area = phi_max_*h2*sqrt(h2 + r2)/(2.0*radius_);

      CCone3DT *th = const_cast<CCone3DT *>(this);

      th->area_.setValue(area);
    }

    return area_.getValue();
  }

 private:
  bool checkPoint(const Point &point) const {
    T z = point.z;

    if (z < 0 || z > height_) return false;

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

    T tv = point.z/height_;

    if (v)
      *v = tv;

    if (dpdu)
      *dpdu = Vector(-phi_max_*point.y, phi_max_*point.x, 0.0);

    if (dpdv) {
      T v1  = 1.0 - tv;
      T iv1 = 1.0/v1;

      *dpdv = Vector(-point.x*iv1, -point.y*iv1, height_);
    }
  }

  T getPhi(const Point &point) const {
    T phi = atan2(point.y, point.x);

    if (phi < 0.0) phi += 2.0*M_PI;

    return phi;
  }
};

typedef CCone3DT<double> CCone3D;

#endif
