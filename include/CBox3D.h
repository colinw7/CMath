#ifndef CBOX_3D_H
#define CBOX_3D_H

#include <CShape3D.h>
#include <CPlane3D.h>

// Box of specified x, y, z sizes centered at origin

class CBox3D : public CShape3D {
 public:
  CBox3D(double rx, double ry, double rz) :
   rx_(rx), ry_(ry), rz_(rz) {
    //pmin_ = CPoint3D(-rx_/2.0, -ry_/2.0, -rz_/2.0);
    //pmax_ = CPoint3D( rx_/2.0,  ry_/2.0,  rz_/2.0);

    pmin_ = CPoint3D(0.0, 0.0, 0.0);
    pmax_ = CPoint3D(rx_, ry_, rz_);
  }

  double getXSize() const { return rx_; }
  double getYSize() const { return ry_; }
  double getZSize() const { return rz_; }

  CBBox3D getBBox() const {
    return CBBox3D(CShape3D::transformFrom(pmin_), CShape3D::transformFrom(pmax_));
  }

  bool intersect(const CLine3D &line, double *tmin, double *tmax) const {
    CPoint3D p1 = CShape3D::transformTo(line.start());
    CPoint3D p2 = CShape3D::transformTo(line.end  ());

    CPoint3D v = CVector3D(p1, p2).point();

    typename CShape3D::TRange trange;

    double t;

    if (v.x != 0) {
      // Plane (x = 0)
      if (v.x > 0)
        t = (    - p1.x)/v.x;
      // Plane (x = rx)
      else
        t = (rx_ - p1.x)/v.x;

      if (trange.isOutside(t)) {
        double y = p1.y + t*v.y;
        double z = p1.z + t*v.z;

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
        double x = p1.x + t*v.x;
        double z = p1.z + t*v.z;

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
        double x = p1.x + t*v.x;
        double y = p1.y + t*v.y;

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

  CVector3D pointNormal(const CPoint3D &point) const {
    CPoint3D p = CShape3D::transformTo(point);

    CVector3D n;

    if      (REAL_EQ(p.x, 0.0)) n = CVector3D(-1, 0, 0);
    else if (REAL_EQ(p.y, 0.0)) n = CVector3D( 0,-1, 0);
    else if (REAL_EQ(p.z, 0.0)) n = CVector3D( 0, 0,-1);
    else if (REAL_EQ(p.x, rx_)) n = CVector3D( 1, 0, 0);
    else if (REAL_EQ(p.y, ry_)) n = CVector3D( 0, 1, 0);
    else if (REAL_EQ(p.z, rz_)) n = CVector3D( 0, 0, 1);
    else                        n = CVector3D( 1, 0, 0);

    return CShape3D::transformFrom(n);
  }

  CVector2D pointToSurfaceVector(const CPoint3D &point) const {
    CPoint3D p = CShape3D::transformTo(point);

    if (REAL_EQ(p.x, 0.0)) return CVector2D(p.y/ry_, p.z/rz_);
    if (REAL_EQ(p.y, 0.0)) return CVector2D(p.y/rx_, p.z/rz_);
    if (REAL_EQ(p.z, 0.0)) return CVector2D(p.y/rx_, p.z/ry_);
    if (REAL_EQ(p.x, rx_)) return CVector2D(p.y/ry_, p.z/rz_);
    if (REAL_EQ(p.y, ry_)) return CVector2D(p.y/rx_, p.z/rz_);
    if (REAL_EQ(p.z, rz_)) return CVector2D(p.y/rx_, p.z/ry_);

    return CVector2D(1, 0);
  }

 private:
  double   rx_, ry_, rz_;
  CPoint3D pmin_, pmax_;
};

#endif
