#ifndef CSPHERE_3D_H
#define CSPHERE_3D_H

#include <CShape3D.h>
#include <CVector3D.h>
#include <CLine3D.h>
#include <CPoint3D.h>
#include <CVector2D.h>
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

template<typename T>
class CSphere3DT : public CShape3DT<T> {
 private:
  typedef typename CShape3DT<T>::BBox BBox;
  typedef CVector3DT<T>               Vector;
  typedef CLine3DT<T>                 Line;
  typedef CPoint3DT<T>                Point;

  T        radius_;      //! Radius R
  T        zmin_, zmax_; //! Height Limits
  T        phi_max_;     //! Angle/Sweep Limit

  T        theta_min_, theta_max_; // Temporaries

  COptReal area_;

 public:
  CSphere3DT(T radius=T(1)) :
   radius_(radius), zmin_(), zmax_(), phi_max_(), theta_min_(), theta_max_(), area_() {
    setZLimit(-radius_, radius_);

    setPhiLimit(360.0);
  }

  T getRadius() const { return radius_; }

  void setRadius(T radius) {
    radius_ = radius;

    area_.setInvalid();
  }

  void setZLimit(T zmin, T zmax) {
    zmin_ = std::min(std::max(std::min(zmin, zmax), -radius_), radius_);
    zmax_ = std::min(std::max(std::max(zmin, zmax), -radius_), radius_);

    theta_min_ = acos(std::min(std::max(zmin_/radius_, -1.0), 1.0));
    theta_max_ = acos(std::min(std::max(zmax_/radius_, -1.0), 1.0));

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

    Point p1 = CShape3D::transformTo(line.start());
    Point p2 = CShape3D::transformTo(line.end  ());

    Line l(p1, p2);

    Vector ro(l.start());

    const Vector &rd = l.vector();

    T a = rd.dotProductSelf();
    T b = 2.0*rd.dotProduct(ro);
    T c = ro.dotProductSelf() - radius_*radius_;

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
      T area = phi_max_*radius_*(zmax_ - zmin_);

      CSphere3DT *th = const_cast<CSphere3DT *>(this);

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

    T theta     = acos(point.z/radius_);
    T theta_len = theta_max_ - theta_min_;

    if (v)
      *v = (theta - theta_min_)/theta_len;

    if (dpdu || dpdv) {
      T zradius2 = point.x*point.x + point.y*point.y;

      if (zradius2 != 0.0) {
        if (dpdu)
          *dpdu = Vector(-phi_max_*point.y, phi_max_*point.x, 0.0);

        if (dpdv) {
          T zradius  = sqrt(zradius2);
          T izradius = 1.0/zradius;

          T cos_phi = point.x*izradius;
          T sin_phi = point.y*izradius;

          *dpdv = theta_len*Vector( point.z*cos_phi,
                                    point.z*sin_phi,
                                   -radius_*sin(theta));
        }
      }
      else {
        T cos_phi = 0;
        T sin_phi = 1;

        Vector tdpdv = theta_len*Vector( point.z*cos_phi,
                                         point.z*sin_phi,
                                        -radius_*sin(theta));

        if (dpdv)
          *dpdv = tdpdv;

        if (dpdu)
          *dpdu = tdpdv.crossProduct(Vector(point));
      }
    }
  }

  T getPhi(const Point &point) const {
    T phi = atan2(point.y, point.x);

    if (phi < 0.0) phi += 2.0*M_PI;

    return phi;
  }
};

typedef CSphere3DT<double> CSphere3D;

#endif
