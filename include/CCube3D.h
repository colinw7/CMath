#ifndef CCUBE_3D_H
#define CCUBE_3D_H

#include <CShape3D.h>
#include <CPlane3D.h>

/* Cube of specified size centered at origin
*/

template<typename T>
class CCube3DT : public CShape3DT<T> {
 private:
  typedef typename CShape3DT<T>::BBox BBox;
  typedef CPlane3DT<T>                Plane;
  typedef CVector3DT<T>               Vector;
  typedef CLine3DT<T>                 Line;
  typedef CPoint3DT<T>                Point;

  T     size_, radius_;
  Point pmin_, pmax_;
  Plane plane1_, plane2_, plane3_, plane4_, plane5_, plane6_;

 public:
  CCube3DT(T size=T(1)) :
   size_(size) {
    radius_ = 0.5*size_;

    pmin_ = Point(-radius_, -radius_, -radius_);
    pmax_ = Point( radius_,  radius_,  radius_);

    plane1_ = createPlane(Vector(-1,  0,  0));
    plane2_ = createPlane(Vector( 1,  0,  0));
    plane3_ = createPlane(Vector( 0, -1,  0));
    plane4_ = createPlane(Vector( 0,  1,  0));
    plane5_ = createPlane(Vector( 0,  0, -1));
    plane6_ = createPlane(Vector( 0,  0,  1));
  }

  T getSize() const { return size_; }

  BBox getBBox() const {
    return BBox(CShape3D::transformFrom(pmin_), CShape3D::transformFrom(pmax_));
  }

  bool intersect(const Line &line, T *tmin, T *tmax) const {
    Point p1 = CShape3D::transformTo(line.start());
    Point p2 = CShape3D::transformTo(line.end  ());

    Line l(p1, p2);

    typename CShape3DT<T>::TRange trange;

    T t;

    CBBox2D bbox_xy(pmin_.x, pmin_.y, pmax_.x, pmax_.y);
    CBBox2D bbox_xz(pmin_.x, pmin_.z, pmax_.x, pmax_.z);
    CBBox2D bbox_yz(pmin_.y, pmin_.z, pmax_.y, pmax_.z);

    if (plane1_.intersect(l, &t) && trange.isOutside(t)) {
      Point p = l.point(t);

      CPoint2D p1(p.y, p.z);

      if (CMathGeom2D::PointInsideRect(p1, bbox_yz))
        trange.update(t);
    }

    if (plane2_.intersect(l, &t) && trange.isOutside(t)) {
      Point p = l.point(t);

      CPoint2D p1(p.y, p.z);

      if (CMathGeom2D::PointInsideRect(p1, bbox_yz))
        trange.update(t);
    }

    if (plane3_.intersect(l, &t) && trange.isOutside(t)) {
      Point p = l.point(t);

      CPoint2D p1(p.x, p.z);

      if (CMathGeom2D::PointInsideRect(p1, bbox_xz))
        trange.update(t);
    }

    if (plane4_.intersect(l, &t) && trange.isOutside(t)) {
      Point p = l.point(t);

      CPoint2D p1(p.x, p.z);

      if (CMathGeom2D::PointInsideRect(p1, bbox_xz))
        trange.update(t);
    }

    if (plane5_.intersect(l, &t) && trange.isOutside(t)) {
      Point p = l.point(t);

      CPoint2D p1(p.x, p.y);

      if (CMathGeom2D::PointInsideRect(p1, bbox_xy))
        trange.update(t);
    }

    if (plane6_.intersect(l, &t) && trange.isOutside(t)) {
      Point p = l.point(t);

      CPoint2D p1(p.x, p.y);

      if (CMathGeom2D::PointInsideRect(p1, bbox_xy))
        trange.update(t);
    }

    if (! trange.set)
      return false;

    *tmin = trange.tmin;
    *tmax = trange.tmax;

    return true;
  }

  Vector pointNormal(const Point &point) const {
    Point p = CShape3D::transformTo(point);

    Vector n;

    if      (fabs(plane1_.value(p)) <= 1E-6) n = plane1_.getNormal();
    else if (fabs(plane2_.value(p)) <= 1E-6) n = plane2_.getNormal();
    else if (fabs(plane3_.value(p)) <= 1E-6) n = plane3_.getNormal();
    else if (fabs(plane4_.value(p)) <= 1E-6) n = plane4_.getNormal();
    else if (fabs(plane5_.value(p)) <= 1E-6) n = plane5_.getNormal();
    else if (fabs(plane6_.value(p)) <= 1E-6) n = plane6_.getNormal();
    else                                     n = Vector(1,0,0);

    return CShape3D::transformFrom(n);
  }

  CVector2D pointToSurfaceVector(const Point &point) const {
    Point p = CShape3D::transformTo(point);

    if (fabs(plane1_.value(p)) <= 1E-6)
      return CVector2D((p.y - pmin_.y)/(pmax_.y - pmin_.y), (p.z - pmin_.z)/(pmax_.z - pmin_.z));
    if (fabs(plane2_.value(p)) <= 1E-6)
      return CVector2D((p.y - pmin_.y)/(pmax_.y - pmin_.y), (p.z - pmin_.z)/(pmax_.z - pmin_.z));
    if (fabs(plane3_.value(p)) <= 1E-6)
      return CVector2D((p.x - pmin_.x)/(pmax_.x - pmin_.x), (p.z - pmin_.z)/(pmax_.z - pmin_.z));
    if (fabs(plane4_.value(p)) <= 1E-6)
      return CVector2D((p.x - pmin_.x)/(pmax_.x - pmin_.x), (p.z - pmin_.z)/(pmax_.z - pmin_.z));
    if (fabs(plane5_.value(p)) <= 1E-6)
      return CVector2D((p.x - pmin_.x)/(pmax_.x - pmin_.x), (p.y - pmin_.y)/(pmax_.y - pmin_.y));
    if (fabs(plane6_.value(p)) <= 1E-6)
      return CVector2D((p.x - pmin_.x)/(pmax_.x - pmin_.x), (p.y - pmin_.y)/(pmax_.y - pmin_.y));

    return CVector2D(0,0);
  }

 private:
  Plane createPlane(const Vector &vector) {
    T r = radius_;

    return Plane((vector*r).point(), vector);
  }
};

typedef CCube3DT<double> CCube3D;

#endif
