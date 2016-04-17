#ifndef CSQUARE_3D_H
#define CSQUARE_3D_H

#include <CShape3D.h>
#include <CPlane3D.h>

/* Square of specified size centered at origin
*/

template<typename T>
class CSquare3DT : public CShape3DT<T> {
 private:
  typedef typename CShape3DT<T>::BBox BBox;
  typedef CPlane3DT<T>                Plane;
  typedef CVector3DT<T>               Vector;
  typedef CLine3DT<T>                 Line;
  typedef CPoint3DT<T>                Point;
  typedef std::vector<Point>          PointList;

  PointList points_;
  Plane     plane_;
  T         xp_[4], yp_[4], zp_[4];

 public:
  struct Size {
    double s1, s2;

    Size(double s1_, double s2_) :
     s1(s1_), s2(s2_) {
    }
  };

  CSquare3DT(const Point &p1=Point(0,0,1),
             const Point &p2=Point(0,0,0),
             const Point &p3=Point(1,0,0)) :
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
    const Point &p1 = points_[0];
    const Point &p2 = points_[1];
    const Point &p3 = points_[2];

    T dx = std::max(fabs(p1.x - p2.x), fabs(p1.x, p3.x));
    T dy = std::max(fabs(p1.y - p2.y), fabs(p1.y, p3.y));
    T dz = std::max(fabs(p1.z - p2.z), fabs(p1.x, p3.z));

    if      (dx > dy && dx > dz) return Size(dx, dz > dy ? dz : dy);
    else if (dy > dz           ) return Size(dy, dx > dz ? dx : dz);
    else                         return Size(dz, dy > dx ? dy : dx);
  }

  const PointList &getPoints() const { return points_; }

  const Point &getPoint1() const { return points_[0]; }
  const Point &getPoint2() const { return points_[1]; }
  const Point &getPoint3() const { return points_[2]; }
  const Point &getPoint4() const { return points_[3]; }

  const Plane &getPlane() const { return plane_; }

  BBox getBBox() const {
    Point p1(std::min(std::min(std::min(xp_[0], xp_[1]), xp_[2]), xp_[3]),
             std::min(std::min(std::min(yp_[0], yp_[1]), yp_[2]), yp_[3]),
             std::min(std::min(std::min(zp_[0], zp_[1]), zp_[2]), zp_[3]));
    Point p2(std::max(std::max(std::max(xp_[0], xp_[1]), xp_[2]), xp_[3]),
             std::max(std::max(std::max(yp_[0], yp_[1]), yp_[2]), yp_[3]),
             std::max(std::max(std::max(zp_[0], zp_[1]), zp_[2]), zp_[3]));

    return BBox(CShape3D::transformFrom(p1), CShape3D::transformFrom(p2));
  }

  bool intersect(const Line &line, T *t) const {
    Point p1 = CShape3D::transformTo(line.start());
    Point p2 = CShape3D::transformTo(line.end  ());

    Line l(p1, p2);

    if (! plane_.intersect(l, t))
      return false;

    Point p = l.point(*t);

    if (! CMathGeom2D::PointInsideConvex(p.x, p.y, xp_, yp_, 4))
      return false;

    return true;
  }

  Vector pointNormal(const Point &) const {
    Vector n = plane_.getNormal();

    return CShape3D::transformFrom(n);
  }

  CVector2D pointToSurfaceVector(const Point &point) const {
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
};

typedef CSquare3DT<double> CSquare3D;

#endif
