#ifndef CLINE_LIST_2D_H
#define CLINE_LIST_2D_H

#include <CShape2D.h>
#include <CPoint2D.h>
#include <CLine2D.h>
#include <CTriangle2D.h>
#include <CMathGeom2D.h>
#include <CBBox2D.h>

class CLineList2D {
 public:
  typedef std::vector<CPoint2D> PointList;

 public:
  CLineList2D() { }

  CLineList2D(const PointList &points) :
   points_(points) {
  }

  CLineList2D(const CPoint2D *points, uint num_points) :
   points_(&points[0], &points[num_points]) {
  }

  CLineList2D(const double *x, const double *y, uint num_xy) :
   points_() {
    for (uint i = 0; i < num_xy; ++i)
      points_.push_back(CPoint2D(x[i], y[i]));
  }

 ~CLineList2D() { }

  const PointList &getPoints() const { return points_; }

  uint getNumPoints() const {
    return points_.size();
  }

  const CPoint2D &getPoint(uint i) const {
    return points_[i];
  }

  CPoint2D getPoint(uint i) {
    return points_[i];
  }

  void getPoint(uint i, double *x, double *y) const {
    *x = points_[i].x;
    *y = points_[i].y;
  }

  void setPoint(uint i, double x, double y) {
    points_[i] = CPoint2D(x, y);
  }

  void setPoint(uint i, const CPoint2D &point) {
    points_[i] = point;
  }

  void addPoint(const CPoint2D &point) {
    points_.push_back(point);
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
  }

  void moveBy(const CPoint2D &p) {
    typename PointList::iterator ps = points_.begin();
    typename PointList::iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps)
      *ps += p;
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
  }

  void rotateBy(double da, const CPoint2D &o) {
    typename PointList::iterator ps = points_.begin();
    typename PointList::iterator pe = points_.end  ();

    for ( ; ps != pe; ++ps)
      *ps = rotatePoint(*ps, da, o);
  }

  CPoint2D rotatePoint(const CPoint2D &point, double da, const CPoint2D &o) {
    double s = sin(da), c = cos(da);

    double x1 = point.x - o.x, y1 = point.y - o.y;
    double x2 = x1*c - y1*s  , y2 = x1*s + y1*c  ;

    return CPoint2D(x2 + o.x, y2 + o.y);
  }

  double includedAngle() const {
    uint np = points_.size();

    assert(np >= 3);

    return CMathGeom2D::IncludedAngle(points_[0], points_[1], points_[2]);
  }

  bool arcThrough(double xr, double yr, double *xc, double *yc,
                  double *xt1, double *yt1, double *xt2, double *yt2) {
    uint np = points_.size();

    assert(np >= 3);

    return CMathGeom2D::ArcThrough(points_[0].x, points_[0].y,
                                   points_[1].x, points_[1].y,
                                   points_[2].x, points_[2].y, xr, yr,
                                   xc, yc, xt1, yt1, xt2, yt2);
  }

 private:
  PointList points_;
};

#endif
