#ifndef CLINE_2D_H
#define CLINE_2D_H

#include <CShape2D.h>
#include <CPoint2D.h>
#include <CVector2D.h>
#include <CMathGen.h>
#include <CBBox2D.h>

class CLine2D : public CShape2D {
 public:
  static void setInsideTolerance(double t) {
    *CLine2D::insideToleranceP() = t;
  }

  static double getInsideTolerance(const double &tol=1E-6) {
    double tolerance = *insideToleranceP();

    if (tolerance < 0)
      tolerance = tol;

    return tolerance;
  }

  CLine2D(double x1=0, double y1=0, double x2=0, double y2=0) :
   p1_(x1, y1), p2_(x2, y2), v_(x2 - x1, y2 - y1) {
  }

  CLine2D(const CPoint2D &p1, const CPoint2D &p2) :
   p1_(p1), p2_(p2), v_(p2_ - p1_) {
  }

  const CPoint2D  &start () const { return p1_; }
  const CPoint2D  &end   () const { return p2_; }
  const CVector2D &vector() const { return v_ ; }

  void setStart(const CPoint2D &p1) {
    p1_ = p1;
    v_  = p2_ - p1_;
  }

  void setEnd(const CPoint2D &p2) {
    p2_ = p2;
    v_  = p2_ - p1_;
  }

  CPoint2D getMid() {
    return CPoint2D((p1_.x + p2_.x)/2, (p1_.y + p2_.y)/2);
  }

  CBBox2D getBBox() const {
    return CBBox2D(p1_, p2_);
  }

  void setBBox(const CBBox2D &bbox) {
    CBBox2D obbox = getBBox();

    if (p1_.x < p2_.x) {
      p1_.x += bbox.getXMin() - obbox.getXMin();
      p2_.x += bbox.getXMax() - obbox.getXMax();
    }
    else {
      p1_.x += bbox.getXMax() - obbox.getXMax();
      p2_.x += bbox.getXMin() - obbox.getXMin();
    }

    if (p1_.y < p2_.y) {
      p1_.y += bbox.getYMin() - obbox.getYMin();
      p2_.y += bbox.getYMax() - obbox.getYMax();
    }
    else {
      p1_.y += bbox.getYMax() - obbox.getYMax();
      p2_.y += bbox.getXMin() - obbox.getYMin();
    }
  }

  bool inside(const CPoint2D &p) const {
    return insideTol(p, 1E-6);
  }

  bool insideTol(const CPoint2D &p, const double &tol) const {
    if ((p.x < std::min(p1_.x, p2_.x) || p.x > std::max(p1_.x, p2_.x)) &&
        (p.y < std::min(p1_.y, p2_.y) || p.y > std::max(p1_.y, p2_.y)))
      return false;

    return pointOn(p, getInsideTolerance(tol));
  }

  void moveBy(const CPoint2D &p) {
    p1_ += p; p2_ += p;
  }

  void resizeBy(const CPoint2D &ll, const CPoint2D &ur) {
    if (p1_.x < p2_.x) {
      p1_.x += ll.x;
      p2_.x += ur.x;
    }
    else {
      p1_.x += ur.x;
      p2_.x += ll.x;
    }

    if (p1_.y < p2_.y) {
      p1_.y += ll.y;
      p2_.y += ur.y;
    }
    else {
      p1_.y += ur.y;
      p2_.y += ll.y;
    }
  }

  void rotateBy(double da, const CPoint2D &o) {
    p1_ = CShape2D::rotatePoint(p1_, da, o);
    p2_ = CShape2D::rotatePoint(p2_, da, o);
  }

  CMathGen::IntersectType
  intersectParms(const CLine2D &line, double *t1, double *t2) const {
    double det = (v_.getX()*line.v_.getY() - v_.getY()*line.v_.getX());

    if (::fabs(det) == 0.0)
      return CMathGen::INTERSECT_NONE;

    double idet = 1.0/det;

    double dx = p1_.x - line.p1_.x;
    double dy = p1_.y - line.p1_.y;

    *t1 = (line.v_.getX()*dy - line.v_.getY()*dx)*idet;
    *t2 = (     v_.getX()*dy -      v_.getY()*dx)*idet;

    if (*t1 >= 0.0 && *t1 <= 1.0 && *t2 >= 0.0 && *t2 <= 1.0)
      return CMathGen::INTERSECT_INSIDE;
    else
      return CMathGen::INTERSECT_OUTSIDE;
  }

  CMathGen::IntersectType
  intersect(const CLine2D &line, CPoint2D &point) const {
    double t1, t2;

    CMathGen::IntersectType type = intersectParms(line, &t1, &t2);

    if (type != CMathGen::INTERSECT_NONE) {
      point.x = p1_.x + v_.getX()*t1;
      point.y = p1_.y + v_.getY()*t1;
    }

    return type;
  }

  CMathGen::IntersectType
  intersect(const CLine2D &line, CPoint2D &point, double *mu1, double *mu2) const {
    CMathGen::IntersectType type = intersectParms(line, mu1, mu2);

    if (type != CMathGen::INTERSECT_NONE) {
      point.x = p1_.x + v_.getX()*(*mu1);
      point.y = p1_.y + v_.getY()*(*mu1);
    }

    return type;
  }

  bool intersects(const CLine2D &line) const;

  bool leftOrOn(const CPoint2D &point) const;

  bool left(const CPoint2D &point) const;

#if 0
  bool intersectsProperly(const CLine2D &line) {
    return CMathGeom2D::IntersectsProperly(p1_, p2_, line.p1_, line.p2_);
  }

  // TODO: tolerance
  bool pointLeft(const CPoint2D &point) {
    return triangleArea2(p1_, p2_, point) > 0.0;
  }

  bool pointLeftOn(const CPoint2D &point) {
    return triangleArea2(p1_, p2_, point) >= 0.0;
  }
#endif

  bool pointOn(const CPoint2D &point, double tol=0.0) const {
    double a = fabs(triangleArea2(p1_, p2_, point));

    if (a <= tol)
      return true;

    return false;
  }

#if 0
  bool pointIn(const CPoint2D &point) {
    if (! pointOn(point))
      return FALSE;

    return pointBetween(point);
  }

  bool pointBetween(const CPoint2D &point) {
    return CMathGeom2D::PointBetween(p1_, p2_, point);
  }
#endif

  bool clip(const CBBox2D &bbox, CLine2D &line) const;

  double lengthSqr() const {
    return v_.lengthSqr();
  }

  double length() const {
    return ::sqrt(lengthSqr());
  }

  CPoint2D interp(double t) const {
    return p1_ + t*v_;
  }

  void split(double t, CLine2D &line1, CLine2D &line2) const {
    CPoint2D pi = interp(t);

    line1 = CLine2D(p1_, pi);
    line2 = CLine2D(pi, p2_);
  }

  // assume already checked point on line
  double getParam(const CPoint2D &p) const {
    double dx = p2_.x - p1_.x;
    double dy = p2_.y - p1_.y;

    if (fabs(dx) > fabs(dy))
      return (p.x - p1_.x)/dx;
    else
      return (p.y - p1_.y)/dy;
  }

  void print(std::ostream &os) const {
    os << p1_ << " " << p2_;
  }

  friend std::ostream &operator<<(std::ostream &os, const CLine2D &line) {
    line.print(os);

    return os;
  }

 private:
  static double *insideToleranceP() {
    static double tolerance = -1;

    return &tolerance;
  }

  static double triangleArea2(const CPoint2D &point1, const CPoint2D &point2,
                              const CPoint2D &point3) {
    return (point2.x - point1.x)*(point3.y - point1.y) -
           (point3.x - point1.x)*(point2.y - point1.y);

  }

 private:
  CPoint2D  p1_, p2_;
  CVector2D v_;
};

//---------

#include <CMathGeom2D.h>

inline bool CLine2D::intersects(const CLine2D &line) const {
  return CMathGeom2D::Intersects(p1_, p2_, line.p1_, line.p2_);
}

inline bool CLine2D::leftOrOn(const CPoint2D &point) const
{
  return CMathGeom2D::PointLineLeftOn(p1_, p2_, point);
}

inline bool CLine2D::left(const CPoint2D &point) const {
  return CMathGeom2D::PointLineLeft(p1_, p2_, point);
}

inline bool CLine2D::clip(const CBBox2D &bbox, CLine2D &line) const {
  double x1 = p1_.x;
  double y1 = p1_.y;
  double x2 = p2_.x;
  double y2 = p2_.y;

  if (CMathGeom2D::clipLine(bbox.getXMin(), bbox.getYMin(), bbox.getXMax(), bbox.getYMax(),
                            &x1, &y1, &x2, &y2)) {
    line = CLine2D(x1, y1, x2, y2);
    return true;
  }

  return false;
}

#endif
