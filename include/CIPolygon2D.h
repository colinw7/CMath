#ifndef CIPOLYGON_2D_H
#define CIPOLYGON_2D_H

#include <CIPoint2D.h>
#include <CILine2D.h>
#include <CMathGeom2D.h>
#include <CIBBox2D.h>
#include <CIterator.h>

template<typename T>
class CIPolygon2DT {
 private:
  typedef CIPolygon2DT<T> Polygon;
  typedef CIPoint2DT<T>   Point;
  typedef CILine2DT<T>    Line;

 public:
  typedef std::vector<Point> PointList;

  //------

 public:
  // iterate polygon points
  class PointIteratorState {
   public:
    PointIteratorState(const Polygon *poly) :
     poly_(poly), len_(poly->getNumPoints()), pos_(0), end_(false) {
      end_ = (pos_ >= len_);
    }

    PointIteratorState() :
     poly_(0), len_(0), pos_(0), end_(true) {
    }

    void next() {
      assert(pos_ < len_);

      ++pos_;

      end_ = (pos_ >= len_);
    }

    const Point &contents() const {
      assert(pos_ < len_);

      return poly_->getPoint(pos_);;
    }

    friend bool operator==(const PointIteratorState &lhs, const PointIteratorState &rhs) {
      if (lhs.end_ == rhs.end_) return true;
      if (lhs.end_ != rhs.end_) return false;

      return (lhs.pos_ == rhs.pos_);
    }

   private:
    const Polygon *poly_;
    int            len_, pos_;
    bool           end_;
  };

  typedef CInputIterator<PointIteratorState, Point> PointIterator;

  //------

 public:
  // iterate polygon lines
  class LineIteratorState {
   public:
    LineIteratorState(const Polygon *poly) :
     poly_(poly), len_(poly->getNumPoints()), pos_(0), end_(false) {
      end_ = (pos_ >= len_);
    }

    LineIteratorState() :
     poly_(0), len_(0), pos_(0), end_(true) {
    }

    void next() {
      assert(pos_ < len_);

      ++pos_;

      end_ = (pos_ >= len_);
    }

    const Line &contents() const {
      assert(pos_ < len_);

      line_.setStart(poly_->getPoint(pos_));
      line_.setEnd  (poly_->getPoint(pos_ < len_ - 1 ? pos_ + 1 : 0));

      return line_;
    }

    friend bool operator==(const LineIteratorState &lhs, const LineIteratorState &rhs) {
      if (lhs.end_ == rhs.end_) return true;
      if (lhs.end_ != rhs.end_) return false;

      return (lhs.pos_ == rhs.pos_);
    }

   private:
    const Polygon *poly_;
    int            len_, pos_;
    bool           end_;
    mutable Line   line_;
  };

  typedef CInputIterator<LineIteratorState, Line> LineIterator;

  //------

 public:
  CIPolygon2DT() :
   array_set_(false), x_(), y_() {
  }

  CIPolygon2DT(const CIPolygon2DT &rhs) :
   points_(rhs.points_), array_set_(false), x_(), y_() {
  }

  CIPolygon2DT &operator=(const CIPolygon2DT &rhs) {
    points_ = rhs.points_;

    array_set_ = false;

    return *this;
  }

  CIPolygon2DT(const PointList &points) :
   points_(points), array_set_(false), x_(), y_() {
  }

  CIPolygon2DT(const Point *points, uint num_points) :
   points_(&points[0], &points[num_points]), array_set_(false), x_(), y_() {
  }

  CIPolygon2DT(const T *x, const T *y, uint num_xy) :
   points_(), array_set_(false), x_(), y_() {
    for (uint i = 0; i < num_xy; ++i)
      points_.push_back(Point(x[i], y[i]));
  }

 ~CIPolygon2DT() { }

  const PointList &getPoints() const { return points_; }

  uint getNumPoints() const { return points_.size(); }

  PointIterator getPointsBegin() const { return PointIterator(PointIteratorState(this)); }
  PointIterator getPointsEnd  () const { return PointIterator();}

  LineIterator getLinesBegin() const { return LineIterator(LineIteratorState(this)); }
  LineIterator getLinesEnd  () const { return LineIterator(); }

  const Point &getPoint(uint i) const {
    return points_[i];
  }

  Point &getPoint(uint i) {
    return points_[i];
  }

  void getPoint(uint i, T *x, T *y) const {
    *x = points_[i].x;
    *y = points_[i].y;
  }

  void clearPoints() {
    points_.clear();

    array_set_ = false;
  }

  void setPoint(uint i, T x, T y) {
    points_[i] = CIPoint2D(x, y);

    array_set_ = false;
  }

  void setPoint(uint i, const Point &point) {
    points_[i] = point;

    array_set_ = false;
  }

  void addPoint(const Point &point) {
    points_.push_back(point);

    array_set_ = false;
  }

  CIBBox2D getBBox() const {
    CIBBox2D bbox;

    typename PointList::const_iterator ps = points_.begin();
    typename PointList::const_iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps)
      bbox.add(*ps);

    return bbox;
  }

  void setBBox(const CIBBox2D &bbox) {
    CIBBox2D obbox = getBBox();

    T sx = bbox.getWidth () / obbox.getWidth ();
    T sy = bbox.getHeight() / obbox.getHeight();

    typename PointList::iterator ps = points_.begin();
    typename PointList::iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps) {
      Point &p = *ps;

      p.x = sx*(p.x - obbox.getXMin()) + bbox.getXMin();
      p.y = sy*(p.y - obbox.getYMin()) + bbox.getYMin();
    }

    array_set_ = false;
  }

  Point centroid() const {
    initArray();

    T cx, cy;

    uint n = points_.size();

    CMathGeom2D::PolygonCentroid(&x_[0], &y_[0], n, &cx, &cy);

    return Point(cx, cy);
  }

  void centroid(T *cx, T *cy) const {
    Point c = centroid();

    *cx = c.x;
    *cy = c.y;
  }

  bool intersect(const CIPolygon2DT &polygon, CIPolygon2DT &ipolygon) const {
    initArray();

    polygon.initArray();

    uint  ni;
    T    *xi, *yi;

    uint n  = points_.size();
    uint n1 = polygon.points_.size();

    if (! CMathGeom2D::IntersectPolygons(&x_[0], &y_[0], n,
                                         &polygon.x_[0], &polygon.y_[0], n1,
                                         &xi, &yi, &ni))
      return false;

    ipolygon = CIPolygon2DT(xi, yi, ni);

    return true;
  }

  bool intersect(const Line &line, std::vector<Point> &ipoints) const {
    return CMathGeom2D::PolygonLineIntersect(points_, line, ipoints);
  }

  bool insideConvex(const Point &point) const {
    return CMathGeom2D::PointInsideConvex(point, points_);
  }

  bool insideEvenOdd(const Point &point) const {
    return CMathGeom2D::PointInsideEvenOdd(point, points_);
  }

  bool inside(const Point &point) const {
    return insideEvenOdd(point);
  }

  void moveBy(const Point &p) {
    typename PointList::iterator ps = points_.begin();
    typename PointList::iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps)
      *ps += p;

    array_set_ = false;
  }

  void resizeBy(const Point &ll, const Point &ur) {
    CIBBox2D bbox = getBBox();

    T w = bbox.getWidth ();
    T h = bbox.getHeight();

    T dw = ur.x - ll.x;
    T dh = ur.y - ll.y;

    T sx = (w > 0 ? (w + dw) / w : 1);
    T sy = (h > 0 ? (h + dh) / h : 1);

    typename PointList::iterator ps = points_.begin();
    typename PointList::iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps) {
      Point &p = *ps;

      p.x = sx*(p.x - bbox.getXMin()) + bbox.getXMin() + ll.x;
      p.y = sy*(p.y - bbox.getYMin()) + bbox.getYMin() + ll.y;
    }

    array_set_ = false;
  }

  void rotateBy(double da, const Point &o) {
    typename PointList::iterator ps = points_.begin();
    typename PointList::iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps)
      *ps = rotatePoint(*ps, da, o);

    array_set_ = false;
  }

  int orientation() const {
    initArray();

    uint n = points_.size();

    return CMathGeom2D::PolygonOrientation(&x_[0], &y_[0], n);
  }

  T area() const {
    initArray();

    uint n = points_.size();

    return CMathGeom2D::PolygonArea(&x_[0], &y_[0], n);
  }

  bool distanceTo(const Point &point, T *dist) const {
    return CMathGeom2D::PointLineDistance(point, *this, dist);
  }

  bool isConvex() const {
    initArray();

    uint n = points_.size();

    return CMathGeom2D::PolygonIsConvex(x_, y_, n);
  }

  const int *getX() const {
    initArray();

    return &x_[0];
  }

  const int *getY() const {
    initArray();

    return &y_[0];
  }

 private:
  void initArray() const {
    if (array_set_) return;

    uint n = points_.size();

    if (n > x_.size()) {
      x_.resize(n);
      y_.resize(n);
    }

    typename PointList::const_iterator ps = points_.begin();
    typename PointList::const_iterator pe = points_.end  ();

    for (uint i = 0; ps != pe; ++ps, ++i)
      (*ps).getXY(&x_[i], &y_[i]);

    array_set_ = true;
  }

 private:
  PointList points_;

  mutable bool           array_set_;
  mutable std::vector<T> x_, y_;
};

typedef CIPolygon2DT<int> CIPolygon2D;

#endif
