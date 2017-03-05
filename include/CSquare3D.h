#ifndef CSQUARE_3D_H
#define CSQUARE_3D_H

#include <CShape3D.h>
#include <CPlane3D.h>

/* Square of specified size centered at origin
*/

class CSquare3D : public CShape3D {
 public:
  typedef std::vector<CPoint3D> PointList;

  struct Size {
    double s1, s2;

    Size(double s1_, double s2_) :
     s1(s1_), s2(s2_) {
    }
  };

 public:
  CSquare3D(const CPoint3D &p1=CPoint3D(0,0,1),
            const CPoint3D &p2=CPoint3D(0,0,0),
            const CPoint3D &p3=CPoint3D(1,0,0)) :
   plane_(p1, p2, p3) {
    points_.resize(4);

    points_[0] = p1;
    points_[1] = p2;
    points_[2] = p3;
    points_[3] = points_[0] + CVector3D(points_[1], points_[2]);

    xp_[0] = points_[0].x; yp_[0] = points_[0].y; zp_[0] = points_[0].z;
    xp_[1] = points_[1].x; yp_[1] = points_[1].y; zp_[1] = points_[1].z;
    xp_[2] = points_[2].x; yp_[2] = points_[2].y; zp_[2] = points_[2].z;
    xp_[3] = points_[3].x; yp_[3] = points_[3].y; zp_[3] = points_[3].z;
  }

  Size getSize() {
    const CPoint3D &p1 = points_[0];
    const CPoint3D &p2 = points_[1];
    const CPoint3D &p3 = points_[2];

    double dx = std::max(fabs(p1.x - p2.x), fabs(p1.x - p3.x));
    double dy = std::max(fabs(p1.y - p2.y), fabs(p1.y - p3.y));
    double dz = std::max(fabs(p1.z - p2.z), fabs(p1.x - p3.z));

    if      (dx > dy && dx > dz) return Size(dx, dz > dy ? dz : dy);
    else if (dy > dz           ) return Size(dy, dx > dz ? dx : dz);
    else                         return Size(dz, dy > dx ? dy : dx);
  }

  const PointList &getPoints() const { return points_; }

  const CPoint3D &getPoint1() const { return points_[0]; }
  const CPoint3D &getPoint2() const { return points_[1]; }
  const CPoint3D &getPoint3() const { return points_[2]; }
  const CPoint3D &getPoint4() const { return points_[3]; }

  const CPlane3D &getPlane() const { return plane_; }

  CBBox3D getBBox() const {
    CPoint3D p1(std::min(std::min(std::min(xp_[0], xp_[1]), xp_[2]), xp_[3]),
                std::min(std::min(std::min(yp_[0], yp_[1]), yp_[2]), yp_[3]),
                std::min(std::min(std::min(zp_[0], zp_[1]), zp_[2]), zp_[3]));
    CPoint3D p2(std::max(std::max(std::max(xp_[0], xp_[1]), xp_[2]), xp_[3]),
                std::max(std::max(std::max(yp_[0], yp_[1]), yp_[2]), yp_[3]),
                std::max(std::max(std::max(zp_[0], zp_[1]), zp_[2]), zp_[3]));

    return CBBox3D(CShape3D::transformFrom(p1), CShape3D::transformFrom(p2));
  }

  bool intersect(const CLine3D &line, double *t) const {
    CPoint3D p1 = CShape3D::transformTo(line.start());
    CPoint3D p2 = CShape3D::transformTo(line.end  ());

    CLine3D l(p1, p2);

    if (! plane_.intersect(l, t))
      return false;

    CPoint3D p = l.point(*t);

    if (! CMathGeom2D::PointInsideConvex(p.x, p.y, xp_, yp_, 4))
      return false;

    return true;
  }

  CVector3D pointNormal(const CPoint3D &) const {
    CVector3D n = plane_.getNormal();

    return CShape3D::transformFrom(n);
  }

  CVector2D pointToSurfaceVector(const CPoint3D &point) const {
    double ux = points_[1].x - points_[0].x;
    double uy = points_[1].y - points_[0].y;
    double vx = points_[3].x - points_[0].x;
    double vy = points_[3].y - points_[0].y;

    double u = 0.0, v = 0.0;

    if      (ux > 0)
      u = (point.x - points_[0].x)/ux;
    else if (vx > 0)
      u = (point.x - points_[0].x)/vx;

    if      (vy > 0)
      v = (point.y - points_[0].y)/vy;
    else if (uy > 0)
      v = (point.y - points_[0].y)/uy;

    u = std::min(std::max(u, 0.0), 1.0);
    v = std::min(std::max(v, 0.0), 1.0);

    return CVector2D(u, v);
  }

 private:
  PointList points_;
  CPlane3D  plane_;
  double    xp_[4], yp_[4], zp_[4];
};

#endif
