#ifndef C2BEZIER_2D_H
#define C2BEZIER_2D_H

#include <CPoint2D.h>
#include <CMathGen.h>
#include <CBBox2D.h>

template<typename T>
class C2Bezier2DT {
 private:
  typedef CPoint2DT<T>   Point;
  typedef C2Bezier2DT<T> Bezier;

 public:
  C2Bezier2DT() :
   p1_(), p2_(), p3_() {
  }

  C2Bezier2DT(T x1, T y1, T x2, T y2, T x3, T y3) :
   p1_(x1, y1), p2_(x2, y2), p3_(x3, y3) {
  }

  C2Bezier2DT(const Point &p1, const Point &p2, const Point &p3) :
   p1_(p1), p2_(p2), p3_(p3) {
  }

  const Point &getFirstPoint  () const { return p1_; }
  const Point &getControlPoint() const { return p2_; }
  const Point &getLastPoint   () const { return p3_; }

  void setFirstPoint  (const Point &p1) { p1_ = p1; }
  void setControlPoint(const Point &p2) { p2_ = p2; };
  void setLastPoint   (const Point &p3) { p3_ = p3; };

  void getFirstPoint  (T *x, T *y) const { *x = p1_.x; *y = p1_.y; }
  void getControlPoint(T *x, T *y) const { *x = p2_.x; *y = p2_.y; }
  void getLastPoint   (T *x, T *y) const { *x = p3_.x; *y = p3_.y; }

  void setFirstPoint  (T x, T y) { setFirstPoint  (Point(x, y)); }
  void setControlPoint(T x, T y) { setControlPoint(Point(x, y)); }
  void setLastPoint   (T x, T y) { setLastPoint   (Point(x, y)); }

  void setPoints(T x1, T y1, T x2, T y2, T x3, T y3) {
    setPoints(Point(x1, y1), Point(x2, y2), Point(x3, y3));
  }

  void setPoints(const Point &p1, const Point &p2, const Point &p3) {
    p1_ = p1; p2_ = p2; p3_ = p3;
  }

  void calc(T t, T *x, T *y) const {
    Point p;

    calc(t, p);

    *x = p.x;
    *y = p.y;
  }

  void calc(T t, Point &p) const {
    T u = (1.0 - t);

    T tt = t*t;
    T uu = u*u;

    p = p1_*uu + 2.0*p2_*t*u + p3_*tt;
  }

  bool interp(T x, T y, T *t) const {
    return interp(Point(x, y), t);
  }

  bool interp(const Point &p, T *t) const {
    T t1 = (::fabs(p.x   - p1_.x) + ::fabs(p.y   - p1_.y))/
           (::fabs(p3_.x - p1_.x) + ::fabs(p3_.y - p1_.y));

    Point pp;

    calc(t1, pp);

    T dx1 = ::fabs(p.x - pp.x);
    T dy1 = ::fabs(p.y - pp.y);

    while (dx1 > 1E-5 || dy1 > 1E-5) {
      if ((pp.x < p.x && pp.x < p3_.x) || (pp.x > p.x && pp.x > p3_.x)) {
        if (pp.x != p3_.x)
          t1 = (1.0 - t1)*(p.x - pp.x)/(p3_.x - pp.x) + t1;
        else
          t1 = 1.0;
      }
      else {
        if (pp.x != p1_.x)
          t1 = t1*(p.x - p1_.x)/(pp.x - p1_.x);
        else
          t1 = 0.0;
      }

      calc(t1, pp);

      T dx2 = ::fabs(p.x - pp.x);
      T dy2 = ::fabs(p.y - pp.y);

      if (dx2 < dx1 && dy2 < dy1)
        goto next;

      if ((pp.y < p.y && pp.y < p3_.y) || (pp.y > p.y && pp.y > p3_.y)) {
        if (pp.y != p3_.y)
          t1 = (1.0 - t1)*(p.y - pp.y)/(p3_.y - pp.y) + t1;
        else
          t1 = 1.0;
      }
      else {
        if (pp.y != p1_.y)
          t1 = t1*(p.y - p1_.y)/(pp.y - p1_.y);
        else
          t1 = 0.0;
      }

      calc(t1, pp);

      dx2 = ::fabs(p.x - pp.x);
      dy2 = ::fabs(p.y - pp.y);

      if (dx2 >= dx1 || dy2 >= dy1)
        return false;

   next:
      dx1 = dx2;
      dy1 = dy2;
    }

    *t = t1;

    return true;
  }

  T gradientStart() const {
    return atan2(p2_.y - p1_.y, p2_.x - p1_.x);
  }

  T gradientEnd() const {
    return atan2(p3_.y - p2_.y, p3_.x - p2_.x);
  }

  T gradient(T t) const {
    T u = 1.0 - t;

    Point p = (p2_ - p1_)*u + (p3_ - p2_)*t;

    T g = atan2(p.y, p.x);

    return g;
  }

  void getHullPolygon(std::vector<Point> &points) const {
    points.push_back(p1_);
    points.push_back(p2_);
    points.push_back(p3_);
  }

  void getHullBBox(CBBox2D &bbox) const {
    bbox.reset();

    bbox.add(p1_);
    bbox.add(p2_);
    bbox.add(p3_);
  }

  void split(Bezier &bezier1, Bezier &bezier2) const {
    // split at control point
    Point p12 = (p1_ + p2_)/2.0;
    Point p23 = (p2_ + p3_)/2.0;

    Point pm = (p12 + p23)/2.0;

    bezier1 = Bezier(p1_, p12, pm );
    bezier2 = Bezier(pm , p23, p3_);
  }

  void print(std::ostream &os) const {
    os << "[[" << p1_.x << ", " << p1_.y << "] [" <<
                  p2_.x << ", " << p2_.y << "] [" <<
                  p3_.x << ", " << p3_.y << "]]";
  }

  friend std::ostream &operator<<(std::ostream &os, const C2Bezier2DT &bezier) {
    bezier.print(os);

    return os;
  }

 private:
  Point p1_, p2_, p3_;
};

typedef C2Bezier2DT<double> C2Bezier2D;

#endif
