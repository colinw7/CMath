#ifndef CIPOLYGON_2D_H
#define CIPOLYGON_2D_H

#include <CIPoint2D.h>
#include <CILine2D.h>
#include <CMathGeom2D.h>
#include <CIBBox2D.h>
#include <CIterator.h>

class CIPolygon2D {
 private:
  typedef CIPolygon2D Polygon;
  typedef CIPoint2D   Point;
  typedef CILine2D    Line;

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
    const Polygon *poly_ { nullptr };
    int            len_ { 0 };
    int            pos_ { 0 };
    bool           end_ { true };
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
     poly_(nullptr), len_(0), pos_(0), end_(true) {
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
    const Polygon *poly_ { nullptr };
    int            len_ { 0 };
    int            pos_ { 0 };
    bool           end_ { true };
    mutable Line   line_;
  };

  typedef CInputIterator<LineIteratorState, Line> LineIterator;

  //------

 public:
  CIPolygon2D() :
   array_set_(false), x_(), y_() {
  }

  CIPolygon2D(const CIPolygon2D &rhs) :
   points_(rhs.points_), array_set_(false), x_(), y_() {
  }

  CIPolygon2D &operator=(const CIPolygon2D &rhs) {
    points_ = rhs.points_;

    array_set_ = false;

    return *this;
  }

  CIPolygon2D(const PointList &points) :
   points_(points), array_set_(false), x_(), y_() {
  }

  CIPolygon2D(const Point *points, uint num_points) :
   points_(&points[0], &points[num_points]), array_set_(false), x_(), y_() {
  }

  CIPolygon2D(const int *x, const int *y, uint num_xy) :
   points_(), array_set_(false), x_(), y_() {
    for (uint i = 0; i < num_xy; ++i)
      points_.push_back(Point(x[i], y[i]));
  }

 ~CIPolygon2D() { }

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

  void getPoint(uint i, int *x, int *y) const {
    *x = points_[i].x;
    *y = points_[i].y;
  }

  void clearPoints() {
    points_.clear();

    array_set_ = false;
  }

  void setPoint(uint i, int x, int y) {
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

    int sx = bbox.getWidth () / obbox.getWidth ();
    int sy = bbox.getHeight() / obbox.getHeight();

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

    int cx, cy;

    uint n = points_.size();

    CMathGeom2D::PolygonCentroid(&x_[0], &y_[0], n, &cx, &cy);

    return Point(cx, cy);
  }

  void centroid(int *cx, int *cy) const {
    Point c = centroid();

    *cx = c.x;
    *cy = c.y;
  }

  bool intersect(const CIPolygon2D &polygon, CIPolygon2D &ipolygon) const {
    initArray();

    polygon.initArray();

    uint  ni { 0 };
    int  *xi { nullptr }, *yi { nullptr };

    uint n  = points_.size();
    uint n1 = polygon.points_.size();

    if (! CMathGeom2D::IntersectPolygons(&x_[0], &y_[0], n,
                                         &polygon.x_[0], &polygon.y_[0], n1,
                                         &xi, &yi, &ni))
      return false;

    ipolygon = CIPolygon2D(xi, yi, ni);

    delete [] xi;
    delete [] yi;

    return true;
  }

#if 0
  bool intersect(const Line &line, std::vector<Point> &ipoints) const {
    return CMathGeom2D::PolygonLineIntersect(points_, line, ipoints);
  }
#endif

#if 0
  bool insideConvex(const Point &point) const {
    return CMathGeom2D::PointInsideConvex(point, points_);
  }

  bool insideEvenOdd(const Point &point) const {
    return CMathGeom2D::PointInsideEvenOdd(point, points_);
  }

  bool inside(const Point &point) const {
    return insideEvenOdd(point);
  }
#endif

  void moveBy(const Point &p) {
    typename PointList::iterator ps = points_.begin();
    typename PointList::iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps)
      *ps += p;

    array_set_ = false;
  }

  void resizeBy(const Point &ll, const Point &ur) {
    CIBBox2D bbox = getBBox();

    int w = bbox.getWidth ();
    int h = bbox.getHeight();

    int dw = ur.x - ll.x;
    int dh = ur.y - ll.y;

    int sx = (w > 0 ? (w + dw) / w : 1);
    int sy = (h > 0 ? (h + dh) / h : 1);

    typename PointList::iterator ps = points_.begin();
    typename PointList::iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps) {
      Point &p = *ps;

      p.x = sx*(p.x - bbox.getXMin()) + bbox.getXMin() + ll.x;
      p.y = sy*(p.y - bbox.getYMin()) + bbox.getYMin() + ll.y;
    }

    array_set_ = false;
  }

#if 0
  void rotateBy(double da, const Point &o) {
    typename PointList::iterator ps = points_.begin();
    typename PointList::iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps)
      *ps = rotatePoint(*ps, da, o);

    array_set_ = false;
  }
#endif

  int orientation() const {
    initArray();

    uint n = points_.size();

    return CMathGeom2D::PolygonOrientation(&x_[0], &y_[0], n);
  }

  int area() const {
    initArray();

    uint n = points_.size();

    return CMathGeom2D::PolygonArea(&x_[0], &y_[0], n);
  }

#if 0
  bool distanceTo(const Point &point, int *dist) const {
    return CMathGeom2D::PointLineDistance(point, *this, dist);
  }
#endif

  bool isConvex() const {
    initArray();

    uint n = points_.size();

    return CMathGeom2D::PolygonIsConvex(&x_[0], &y_[0], n);
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

  mutable bool             array_set_ { false };
  mutable std::vector<int> x_, y_;
};

#endif
