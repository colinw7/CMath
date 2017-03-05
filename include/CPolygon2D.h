#ifndef CPOLYGON_2D_H
#define CPOLYGON_2D_H

#include <CShape2D.h>
#include <CPoint2D.h>
#include <CLine2D.h>
#include <CTriangle2D.h>
#include <CMathGeom2D.h>
#include <CBBox2D.h>
#include <CStateIterator.h>

class CPolygon2D : public CShape2D {
 private:
  struct EarPoint {
    const CPoint2D &point;
    bool            is_ear { false };

    EarPoint(const CPoint2D &p) :
     point(p), is_ear(false) {
    }

    bool operator==(const EarPoint &ep) {
      return (&point == &ep.point);
    }
  };

 public:
  typedef std::vector<CPoint2D> PointList;
  typedef std::vector<CLine2D>  LineList;
  typedef std::list<EarPoint>   EarPointList;

  //------

 public:
  // iterate polygon points
  class PointIteratorState {
   public:
    PointIteratorState(const CPolygon2D *poly) :
     poly_(poly), len_(poly->getNumPoints()), pos_(0), end_(false) {
      end_ = (pos_ >= len_);
    }

    PointIteratorState() :
     poly_(nullptr), len_(0), pos_(0), end_(true) {
    }

    void next() {
      assert(pos_ < len_);

      ++pos_;

      end_ = (pos_ >= len_);
    }

    const CPoint2D &contents() const {
      assert(pos_ < len_);

      return poly_->getPoint(pos_);;
    }

    friend bool operator==(const PointIteratorState &lhs, const PointIteratorState &rhs) {
      if (lhs.end_ == rhs.end_) return true;
      if (lhs.end_ != rhs.end_) return false;

      return (lhs.pos_ == rhs.pos_);
    }

   private:
    const CPolygon2D *poly_ { nullptr };
    int               len_ { 0 };
    int               pos_ { 0 };
    bool              end_ { true };
  };

  typedef CInputIterator<PointIteratorState, CPoint2D> PointIterator;

  //------

 public:
  // iterate polygon lines
  class LineIteratorState {
   public:
    LineIteratorState(const CPolygon2D *poly) :
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

    const CLine2D &contents() const {
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
    const CPolygon2D *poly_;
    int               len_, pos_;
    bool              end_;
    mutable CLine2D   line_;
  };

  typedef CInputIterator<LineIteratorState, CLine2D> LineIterator;

  //------

 public:
  CPolygon2D() :
   array_set_(false), x_(), y_() {
  }

  CPolygon2D(const CPolygon2D &rhs) :
   CShape2D(rhs), points_(rhs.points_), array_set_(false), x_(), y_() {
  }

  CPolygon2D &operator=(const CPolygon2D &rhs) {
    points_ = rhs.points_;

    array_set_ = false;

    return *this;
  }

  CPolygon2D(const PointList &points) :
   points_(points), array_set_(false), x_(), y_() {
  }

  CPolygon2D(const CPoint2D *points, uint num_points) :
   points_(&points[0], &points[num_points]), array_set_(false), x_(), y_() {
  }

  CPolygon2D(const double *x, const double *y, uint num_xy) :
   points_(), array_set_(false), x_(), y_() {
    for (uint i = 0; i < num_xy; ++i)
      points_.push_back(CPoint2D(x[i], y[i]));
  }

 ~CPolygon2D() { }

  const PointList &getPoints() const { return points_; }

  uint getNumPoints() const { return points_.size(); }

  PointIterator getPointsBegin() const { return PointIterator(PointIteratorState(this)); }
  PointIterator getPointsEnd  () const { return PointIterator();}

  LineIterator getLinesBegin() const { return LineIterator(LineIteratorState(this)); }
  LineIterator getLinesEnd  () const { return LineIterator(); }

  const CPoint2D &getPoint(uint i) const {
    return points_[i];
  }

  CPoint2D &getPoint(uint i) {
    return points_[i];
  }

  void getPoint(uint i, double *x, double *y) const {
    *x = points_[i].x;
    *y = points_[i].y;
  }

  void clearPoints() {
    points_.clear();

    array_set_ = false;
  }

  void setPoint(uint i, double x, double y) {
    points_[i] = CPoint2D(x, y);

    array_set_ = false;
  }

  void setPoint(uint i, const CPoint2D &point) {
    points_[i] = point;

    array_set_ = false;
  }

  void addPoint(const CPoint2D &point) {
    points_.push_back(point);

    array_set_ = false;
  }

  CBBox2D getBBox() const {
    CBBox2D bbox;

    typename PointList::const_iterator ps = points_.begin();
    typename PointList::const_iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps)
      bbox.add(*ps);

    return bbox;
  }

  void setBBox(const CBBox2D &bbox) {
    CBBox2D obbox = getBBox();

    double sx = bbox.getWidth () / obbox.getWidth ();
    double sy = bbox.getHeight() / obbox.getHeight();

    typename PointList::iterator ps = points_.begin();
    typename PointList::iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps) {
      CPoint2D &p = *ps;

      p.x = sx*(p.x - obbox.getXMin()) + bbox.getXMin();
      p.y = sy*(p.y - obbox.getYMin()) + bbox.getYMin();
    }

    array_set_ = false;
  }

  CPoint2D centroid() const {
    initArray();

    double cx, cy;

    uint n = points_.size();

    CMathGeom2D::PolygonCentroid(&x_[0], &y_[0], n, &cx, &cy);

    return CPoint2D(cx, cy);
  }

  void centroid(double *cx, double *cy) const {
    CPoint2D c = centroid();

    *cx = c.x;
    *cy = c.y;
  }

  bool intersect(const CPolygon2D &polygon, CPolygon2D &ipolygon) const {
    initArray();

    polygon.initArray();

    uint    ni { 0 };
    double *xi { nullptr }, *yi { nullptr };

    uint n  = points_.size();
    uint n1 = polygon.points_.size();

    if (! CMathGeom2D::IntersectPolygons(&x_[0], &y_[0], n,
                                         &polygon.x_[0], &polygon.y_[0], n1,
                                         &xi, &yi, &ni))
      return false;

    ipolygon = CPolygon2D(xi, yi, ni);

    delete [] xi;
    delete [] yi;

    return true;
  }

  bool intersect(const CLine2D &line, std::vector<CPoint2D> &ipoints) const {
    return CMathGeom2D::PolygonLineIntersect(points_, line, ipoints);
  }

  bool insideConvex(const CPoint2D &point) const {
    return CMathGeom2D::PointInsideConvex(point, points_);
  }

  bool insideEvenOdd(const CPoint2D &point) const {
    return CMathGeom2D::PointInsideEvenOdd(point, points_);
  }

  bool inside(const CPoint2D &point) const {
    return insideEvenOdd(point);
  }

  void moveBy(const CPoint2D &p) {
    typename PointList::iterator ps = points_.begin();
    typename PointList::iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps)
      *ps += p;

    array_set_ = false;
  }

  void resizeBy(const CPoint2D &ll, const CPoint2D &ur) {
    CBBox2D bbox = getBBox();

    double w = bbox.getWidth ();
    double h = bbox.getHeight();

    double dw = ur.x - ll.x;
    double dh = ur.y - ll.y;

    double sx = (w > 0 ? (w + dw) / w : 1);
    double sy = (h > 0 ? (h + dh) / h : 1);

    typename PointList::iterator ps = points_.begin();
    typename PointList::iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps) {
      CPoint2D &p = *ps;

      p.x = sx*(p.x - bbox.getXMin()) + bbox.getXMin() + ll.x;
      p.y = sy*(p.y - bbox.getYMin()) + bbox.getYMin() + ll.y;
    }

    array_set_ = false;
  }

  void rotateBy(double da, const CPoint2D &o) {
    typename PointList::iterator ps = points_.begin();
    typename PointList::iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps)
      *ps = CShape2D::rotatePoint(*ps, da, o);

    array_set_ = false;
  }

  int orientation() const {
    initArray();

    uint n = points_.size();

    return CMathGeom2D::PolygonOrientation(&x_[0], &y_[0], n);
  }

  double area() const {
    initArray();

    uint n = points_.size();

    return CMathGeom2D::PolygonArea(&x_[0], &y_[0], n);
  }

  bool distanceTo(const CPoint2D &point, double *dist) const {
    int n = points_.size();

    double minDist { 0.0 };
    bool   minSet  { false };

    for (int i1 = n - 1, i2 = 0; i2 < n; i1 = i2++) {
      CLine2D line(points_[i1], points_[i2]);

      double dist1;

      if (CMathGeom2D::PointLineDistance(point, line, &dist1)) {
        if (minSet) {
          minDist = dist1;
          minSet  = true;
        }
        else
          minDist = std::min(minDist, dist1);
      }
    }

    *dist = minDist;

    return minSet;
  }

  bool isConvex() const {
    initArray();

    uint n = points_.size();

    return CMathGeom2D::PolygonIsConvex(&x_[0], &y_[0], n);
  }

  void triangulate(std::vector<CTriangle2D> &triangle_list) const {
    EarPointList ear_points;

    typename PointList::const_iterator ps = points_.begin();
    typename PointList::const_iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps)
      ear_points.push_back(EarPoint(*ps));

    typename EarPointList::iterator eps = ear_points.begin();
    typename EarPointList::iterator epe = ear_points.end  ();

    typename EarPointList::iterator ep1, ep2, ep3, ep4, ep5;

    ep3 = eps;
    ep1 = ep3++;
    ep2 = ep3++;

    do {
      (*ep2).is_ear = isDiagonal(ear_points, ep1, ep3);

      ep1 = ep2;
      ep2 = ep3++;

      if (ep3 == epe) ep3 = eps;
    } while (ep1 != eps);

    while (ear_points.size() > 3) {
      ep5 = eps;
      ep1 = ep5++;
      ep2 = ep5++;
      ep3 = ep5++;
      ep4 = ep5++;

      if (ep5 == epe) ep5 = eps;

      do {
        if ((*ep3).is_ear) {
          addTriangle(triangle_list, ep3, ep2, ep4);

          (*ep2).is_ear = isDiagonal(ear_points, ep1, ep4);
          (*ep4).is_ear = isDiagonal(ear_points, ep2, ep5);

          ear_points.erase(ep3);

          eps = ear_points.begin();
          epe = ear_points.end  ();

          break;
        }

        ep1 = ep2;
        ep2 = ep3;
        ep3 = ep4;
        ep4 = ep5++;

        if (ep5 == epe) ep5 = eps;
      } while (ep1 != eps);
    }

    ep3 = eps;
    ep1 = ep3++;
    ep2 = ep3++;

    addTriangle(triangle_list, ep1, ep2, ep3);
  }

  const double *getX() const {
    initArray();

    return &x_[0];
  }

  const double *getY() const {
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
  static void addTriangle(std::vector<CTriangle2D> &triangle_list,
                          const typename EarPointList::const_iterator &ep1,
                          const typename EarPointList::const_iterator &ep2,
                          const typename EarPointList::const_iterator &ep3) {
    //std::cout << "Diagonal " <<
    //             (*ep2).point << "->" << (*ep3).point << std::endl;

    triangle_list.push_back(CTriangle2D((*ep1).point, (*ep2).point, (*ep3).point));
  }

  static bool isDiagonal(const EarPointList &ear_points,
                         const typename EarPointList::const_iterator &epa,
                         const typename EarPointList::const_iterator &epb) {
    return inCone(ear_points, epa, epb) &&
           inCone(ear_points, epb, epa) &&
           isDiagonalInOut(ear_points, epa, epb);
  }

  static bool isDiagonalInOut(const EarPointList &ear_points,
                              const typename EarPointList::const_iterator &epa,
                              const typename EarPointList::const_iterator &epb) {
    typename EarPointList::const_iterator eps = ear_points.begin();
    typename EarPointList::const_iterator epe = ear_points.end  ();

    typename EarPointList::const_iterator ep2 = eps;
    typename EarPointList::const_iterator ep1 = ep2++;

    do {
      if (ep1 != epa && ep2 != epa && ep1 != epb && ep2 != epb) {
        if (Intersects(epa, epb, ep1, ep2))
          return false;
      }

      ep1 = ep2++;

      if (ep2 == epe) ep2 = eps;
    } while (ep1 != eps);

    return true;
  }

  static bool inCone(const EarPointList &ear_points,
                     const typename EarPointList::const_iterator &ep1,
                     const typename EarPointList::const_iterator &epb) {
    typename EarPointList::const_iterator eps = ear_points.begin();
    typename EarPointList::const_iterator epe = ear_points.end();

    typename EarPointList::const_iterator ep0 = ep1; --ep0;
    typename EarPointList::const_iterator ep2 = ep1; ++ep2;

    if (ep0 == epe) ep0 = (++ear_points.rbegin()).base(); // last element
    if (ep2 == epe) ep2 = eps;

    if (PointLineLeftOn(ep1, ep2, ep0))
      return    PointLineLeft(ep1, epb, ep0) && PointLineLeft(epb, ep1, ep2);
    else
      return ! (PointLineLeft(ep1, epb, ep2) && PointLineLeft(epb, ep1, ep0));
  }

  static bool Intersects(const typename EarPointList::const_iterator &ep1,
                         const typename EarPointList::const_iterator &ep2,
                         const typename EarPointList::const_iterator &ep3,
                         const typename EarPointList::const_iterator &ep4) {
    return CMathGeom2D::Intersects((*ep1).point, (*ep2).point, (*ep3).point, (*ep4).point);
  }

  static bool PointLineLeftOn(const typename EarPointList::const_iterator &ep1,
                              const typename EarPointList::const_iterator &ep2,
                              const typename EarPointList::const_iterator &ep3) {
    return CMathGeom2D::PointLineLeftOn((*ep1).point, (*ep2).point, (*ep3).point);
  }

  static bool PointLineLeft(const typename EarPointList::const_iterator &ep1,
                            const typename EarPointList::const_iterator &ep2,
                            const typename EarPointList::const_iterator &ep3) {
    return CMathGeom2D::PointLineLeft((*ep1).point, (*ep2).point, (*ep3).point);
  }

 private:
  PointList points_;

  mutable bool                array_set_ { false };
  mutable std::vector<double> x_, y_;
};

#endif
