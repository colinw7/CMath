#ifndef CBOX_3D_H
#define CBOX_3D_H

#include <CShape3D.h>
#include <CPlane3D.h>

// Box of specified x, y, z sizes centered at origin

template<typename T>
class CBox3DT : public CShape3DT<T> {
 private:
  typedef typename CShape3DT<T>::BBox BBox;
  typedef CVector3DT<T>               Vector;
  typedef CLine3DT<T>                 Line;
  typedef CPoint3DT<T>                Point;

  T     rx_, ry_, rz_;
  Point pmin_, pmax_;

 public:
  CBox3DT(T rx, T ry, T rz) :
   rx_(rx), ry_(ry), rz_(rz) {
    pmin_ = Point(-rx_/2.0, -ry_/2.0, -rz_/2.0);
    pmax_ = Point( rx_/2.0,  ry_/2.0,  rz_/2.0);
  }

  T getXSize() const { return rx_; }
  T getYSize() const { return ry_; }
  T getZSize() const { return rz_; }

  BBox getBBox() const {
    return BBox(CShape3D::transformFrom(pmin_), CShape3D::transformFrom(pmax_));
  }

  bool intersect(const Line &line, T *tmin, T *tmax) const {
    Point p1 = CShape3D::transformTo(line.start());
    Point p2 = CShape3D::transformTo(line.end  ());

    Point v = Vector(p1, p2).point();

    typename CShape3DT<T>::TRange trange;

    T t;

    if (v.x != 0) {
      // Plane (x = 0)
      if (v.x > 0)
        t = (    - p1.x)/v.x;
      // Plane (x = rx)
      else
        t = (rx_ - p1.x)/v.x;

      if (trange.isOutside(t)) {
        T y = p1.y + t*v.y;
        T z = p1.z + t*v.z;

        if (y >= 0 && y <= ry_ && z >= 0 && z <= rz_)
          trange.update(t);
      }
    }

    if (v.y != 0) {
      // Plane (y = 0)
      if      (v.y > 0)
        t = (    - p1.y)/v.y;
      // Plane (y = ry)
      else
        t = (ry_ - p1.y)/v.y;

      if (trange.isOutside(t)) {
        T x = p1.x + t*v.x;
        T z = p1.z + t*v.z;

        if (x >= 0 && x <= rx_ && z >= 0 && z <= rz_)
          trange.update(t);
      }
    }

    if (v.z != 0) {
      // Plane (z = 0)
      if      (v.z > 0)
        t = (    - p1.z)/v.z;
      // Plane (z = rz)
      else if (v.z < 0)
        t = (rz_ - p1.z)/v.z;

      if (trange.isOutside(t)) {
        T x = p1.x + t*v.x;
        T y = p1.y + t*v.y;

        if (x >= 0 && x <= rx_ && y >= 0 && y <= ry_)
          trange.update(t);
      }
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

    if      (REAL_EQ(p.x, 0.0)) n = Vector(-1, 0, 0);
    else if (REAL_EQ(p.y, 0.0)) n = Vector( 0,-1, 0);
    else if (REAL_EQ(p.z, 0.0)) n = Vector( 0, 0,-1);
    else if (REAL_EQ(p.x, rx_)) n = Vector( 1, 0, 0);
    else if (REAL_EQ(p.y, ry_)) n = Vector( 0, 1, 0);
    else if (REAL_EQ(p.z, rz_)) n = Vector( 0, 0, 1);
    else                        n = Vector( 1, 0, 0);

    return CShape3D::transformFrom(n);
  }

  CVector2D pointToSurfaceVector(const Point &point) const {
    Point p = CShape3D::transformTo(point);

    if (REAL_EQ(p.x, 0.0)) return CVector2D(p.y/ry_, p.z/rz_);
    if (REAL_EQ(p.y, 0.0)) return CVector2D(p.y/rx_, p.z/rz_);
    if (REAL_EQ(p.z, 0.0)) return CVector2D(p.y/rx_, p.z/ry_);
    if (REAL_EQ(p.x, rx_)) return CVector2D(p.y/ry_, p.z/rz_);
    if (REAL_EQ(p.y, ry_)) return CVector2D(p.y/rx_, p.z/rz_);
    if (REAL_EQ(p.z, rz_)) return CVector2D(p.y/rx_, p.z/ry_);

    return CVector2D(1, 0);
  }
};

typedef CBox3DT<double> CBox3D;

#endif
