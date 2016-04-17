#ifndef CCYLINDER_3D_H
#define CCYLINDER_3D_H

#include <CShape3D.h>
#include <CVector3D.h>
#include <CLine3D.h>
#include <CPoint3D.h>
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

template<typename T>
class CCylinder3DT : public CShape3DT<T> {
 private:
  typedef typename CShape3DT<T>::BBox BBox;
  typedef CVector3DT<T>               Vector;
  typedef CLine3DT<T>                 Line;
  typedef CPoint3DT<T>                Point;

  T        radius_; //! Radius r
  T        height_; //! Height h

  // limits
  T        zmin_, zmax_; //! Height Range
  T        phi_max_;     //! Angle/Sweep Max

  COptReal area_;

 public:
  CCylinder3DT(T radius=T(1), T height=T(1)) :
   radius_(radius), height_(height) {
    setZLimit(0, height_);

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

  void setZLimit(T zmin, T zmax) {
    zmin_ = std::min(std::max(std::min(zmin, zmax), 0.0), height_);
    zmax_ = std::min(std::max(std::max(zmin, zmax), 0.0), height_);

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

    Point p1 = CShape3D::transformTo(line.start());
    Point p2 = CShape3D::transformTo(line.end  ());

    Line l(p1, p2);

    Vector ro(l.start());

    const Vector &rd = l.vector();

    T rox = ro.getX(); T roy = ro.getY();
    T rdx = rd.getX(); T rdy = rd.getY();

    T a = rdx*rdx + rdy*rdy;
    T b = 2.0*(rdx*rox + rdy*roy);
    T c = rox*rox + roy*roy - radius_*radius_;

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
      T area = (zmax_ - zmin_)*phi_max_*radius_;

      CCylinder3DT *th = const_cast<CCylinder3DT *>(this);

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

    if (dpdv)
      *dpdv = Vector(0, 0, zmax_ - zmin_);
  }

  T getPhi(const Point &point) const {
    T phi = atan2(point.y, point.x);

    if (phi < 0.0) phi += 2.0*M_PI;

    return phi;
  }
};

typedef CCylinder3DT<double> CCylinder3D;

#endif
