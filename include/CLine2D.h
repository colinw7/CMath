#ifndef CLINE_2D_H
#define CLINE_2D_H

#include <CShape2D.h>
#include <CPoint2D.h>
#include <CVector2D.h>
#include <CMathGen.h>
#include <CBBox2D.h>

template<class T>
class CLine2DT : public CShape2DT<T> {
 private:
  typedef CLine2DT<T>   Line;
  typedef CPoint2DT<T>  Point;
  typedef CVector2DT<T> Vector;
  typedef CBBox2DT<T>   BBox;

 public:
  static void setInsideTolerance(T t) {
    *CLine2DT<T>::insideToleranceP() = t;
  }

  static T getInsideTolerance() {
    return *CLine2DT<T>::insideToleranceP();
  }

  CLine2DT(T x1=0, T y1=0, T x2=0, T y2=0) :
   p1_(x1, y1), p2_(x2, y2), v_(x2 - x1, y2 - y1) {
  }

  CLine2DT(const Point &p1, const Point &p2) :
   p1_(p1), p2_(p2), v_(p2_ - p1_) {
  }

  const Point  &start () const { return p1_; }
  const Point  &end   () const { return p2_; }
  const Vector &vector() const { return v_ ; }

  void setStart(const Point &p1) {
    p1_ = p1;
    v_  = p2_ - p1_;
  }

  void setEnd(const Point &p2) {
    p2_ = p2;
    v_  = p2_ - p1_;
  }

  Point getMid() {
    return Point((p1_.x + p2_.x)/2, (p1_.y + p2_.y)/2);
  }

  BBox getBBox() const {
    return BBox(p1_, p2_);
  }

  void setBBox(const BBox &bbox) {
    BBox obbox = getBBox();

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

  bool inside(const Point &p) const {
    if ((p.x < std::min(p1_.x, p2_.x) || p.x > std::max(p1_.x, p2_.x)) &&
        (p.y < std::min(p1_.y, p2_.y) || p.y > std::max(p1_.y, p2_.y)))
      return false;

    return pointOn(p, getInsideTolerance());
  }

  void moveBy(const Point &p) {
    p1_ += p; p2_ += p;
  }

  void resizeBy(const Point &ll, const Point &ur) {
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

  void rotateBy(double da, const Point &o) {
    p1_ = CShape2D::rotatePoint(p1_, da, o);
    p2_ = CShape2D::rotatePoint(p2_, da, o);
  }

  CMathGen::IntersectType
  intersectParms(const Line &line, T *t1, T *t2) const {
    T det = (v_.getX()*line.v_.getY() - v_.getY()*line.v_.getX());

    if (::fabs(det) == 0.0)
      return CMathGen::INTERSECT_NONE;

    T idet = 1.0/det;

    T dx = p1_.x - line.p1_.x;
    T dy = p1_.y - line.p1_.y;

    *t1 = (line.v_.getX()*dy - line.v_.getY()*dx)*idet;
    *t2 = (     v_.getX()*dy -      v_.getY()*dx)*idet;

    if (*t1 >= 0.0 && *t1 <= 1.0 && *t2 >= 0.0 && *t2 <= 1.0)
      return CMathGen::INTERSECT_INSIDE;
    else
      return CMathGen::INTERSECT_OUTSIDE;
  }

  CMathGen::IntersectType
  intersect(const Line &line, Point &point) const {
    T t1, t2;

    CMathGen::IntersectType type = intersectParms(line, &t1, &t2);

    if (type != CMathGen::INTERSECT_NONE) {
      point.x = p1_.x + v_.getX()*t1;
      point.y = p1_.y + v_.getY()*t1;
    }

    return type;
  }

  CMathGen::IntersectType
  intersect(const Line &line, Point &point, T *mu1, T *mu2) const {
    CMathGen::IntersectType type = intersectParms(line, mu1, mu2);

    if (type != CMathGen::INTERSECT_NONE) {
      point.x = p1_.x + v_.getX()*(*mu1);
      point.y = p1_.y + v_.getY()*(*mu1);
    }

    return type;
  }

  bool intersects(const Line &line) const;

  bool leftOrOn(const Point &point) const;

  bool left(const Point &point) const;

#if 0
  bool intersectsProperly(const Line &line) {
    return CMathGeom2D::IntersectsProperly(p1_, p2_, line.p1_, line.p2_);
  }

  // TODO: tolerance
  bool pointLeft(const Point &point) {
    return triangleArea2(p1_, p2_, point) > 0.0;
  }

  bool pointLeftOn(const Point &point) {
    return triangleArea2(p1_, p2_, point) >= 0.0;
  }
#endif

  bool pointOn(const Point &point, T tol=0.0) const {
    return triangleArea2(p1_, p2_, point) <= tol;
  }

#if 0
  bool pointIn(const Point &point) {
    if (! pointOn(point))
      return FALSE;

    return pointBetween(point);
  }

  bool pointBetween(const Point &point) {
    return CMathGeom2D::PointBetween(p1_, p2_, point);
  }
#endif

  bool clip(const BBox &bbox, Line &line) const;

  T lengthSqr() const {
    return v_.lengthSqr();
  }

  T length() const {
    return ::sqrt(lengthSqr());
  }

  Point interp(double t) const {
    return p1_ + t*v_;
  }

  void split(double t, Line &line1, Line &line2) const {
    Point pi = interp(t);

    line1 = Line(p1_, pi);
    line2 = Line(pi, p2_);
  }

  // assume already checked point on line
  double getParam(const Point &p) const {
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

  friend std::ostream &operator<<(std::ostream &os, const Line &line) {
    line.print(os);

    return os;
  }

 private:
  static T *insideToleranceP() {
    static T tolerance;

    return &tolerance;
  }

  static double triangleArea2(const Point &point1, const Point &point2, const Point &point3) {
    return (point2.x - point1.x)*(point3.y - point1.y) -
           (point3.x - point1.x)*(point2.y - point1.y);

  }

 private:
  Point  p1_, p2_;
  Vector v_;
};

typedef CLine2DT<double> CLine2D;

//---------

#include <CMathGeom2D.h>

template<class T>
bool
CLine2DT<T>::
intersects(const Line &line) const
{
  return CMathGeom2D::Intersects(p1_, p2_, line.p1_, line.p2_);
}

template<class T>
bool
CLine2DT<T>::
leftOrOn(const Point &point) const
{
  return CMathGeom2D::PointLineLeftOn(p1_, p2_, point);
}

template<class T>
bool
CLine2DT<T>::
left(const Point &point) const
{
  return CMathGeom2D::PointLineLeft(p1_, p2_, point);
}

template<class T>
bool
CLine2DT<T>::
clip(const BBox &bbox, Line &line) const
{
  double x1 = p1_.x;
  double y1 = p1_.y;
  double x2 = p2_.x;
  double y2 = p2_.y;

  if (CMathGeom2D::clipLine(bbox.getXMin(), bbox.getYMin(), bbox.getXMax(), bbox.getYMax(),
                            &x1, &y1, &x2, &y2)) {
    line = Line(x1, y1, x2, y2);
    return true;
  }

  return false;
}

#endif
