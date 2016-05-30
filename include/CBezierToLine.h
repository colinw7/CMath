#ifndef CBEZIER_TO_LINE_H
#define CBEZIER_TO_LINE_H

#include <C3Bezier2D.h>
#include <C2Bezier2D.h>
#include <vector>

class CBezierToLine {
 public:
  typedef std::vector<CPoint2D> PointList;

  enum { MAX_BEZIER_DEPTH = 9 };

 public:
  CBezierToLine() { }

  virtual ~CBezierToLine() { }

  void setTolerance(double tol) { tol_ = tol; }

  void calc(const C3Bezier2D &bezier) {
    init(bezier);

    //----

    if (bezierPoints_.size() <= MAX_BEZIER_DEPTH)
      bezierPoints_.resize(MAX_BEZIER_DEPTH + 1);

    uint num_points = 0;

    bezier.getFirstPoint(&bezierPoints_[0].x, &bezierPoints_[0].y);

    ++num_points;

    calc1(bezier, &num_points, 1);

    //----

    bezierPoints_.resize(num_points);
  }

  void calc(const C2Bezier2D &bezier) {
    init(bezier);

    //----

    if (bezierPoints_.size() <= MAX_BEZIER_DEPTH)
      bezierPoints_.resize(MAX_BEZIER_DEPTH + 1);

    uint num_points = 0;

    bezier.getFirstPoint(&bezierPoints_[0].x, &bezierPoints_[0].y);

    ++num_points;

    calc1(bezier, &num_points, 1);

    //----

    bezierPoints_.resize(num_points);
  }

  virtual void init(const C3Bezier2D &bezier) {
    double x1, y1, x4, y4;

    bezier.getFirstPoint(&x1, &y1);
    bezier.getLastPoint (&x4, &y4);

    if (tol_ < 0.0)
      tol_ = std::max(fabs(x4 - x1), fabs(y4 - y1))/20;
  }

  virtual void init(const C2Bezier2D &bezier) {
    double x1, y1, x3, y3;

    bezier.getFirstPoint(&x1, &y1);
    bezier.getLastPoint (&x3, &y3);

    if (tol_ < 0.0)
      tol_ = std::max(fabs(x3 - x1), fabs(y3 - y1))/20;
  }

  virtual bool checkLength(double s) {
    return (s <= tol_);
  }

  uint getNumPoints() const { return bezierPoints_.size(); }

  const CPoint2D &getPoint(uint i) const { return bezierPoints_[i]; }

  void toLines(const C2Bezier2D &bezier, std::vector<CPoint2D> &points) {
    calc(bezier);

    uint num_points = getNumPoints();

    points.resize(num_points);

    for (uint i = 0; i < num_points; ++i)
      points[i] = getPoint(i);
  }

  void toLines(const C3Bezier2D &bezier, std::vector<CPoint2D> &points) {
    calc(bezier);

    uint num_points = getNumPoints();

    points.resize(num_points);

    for (uint i = 0; i < num_points; ++i)
      points[i] = getPoint(i);
  }

 private:
  void calc1(const C3Bezier2D &bezier, uint *num_points, uint depth) {
    if (*num_points >= bezierPoints_.size())
      bezierPoints_.resize(2*bezierPoints_.size() + *num_points + 1);

    //-----

    if (depth >= MAX_BEZIER_DEPTH) {
      bezier.getLastPoint(&bezierPoints_[*num_points].x, &bezierPoints_[*num_points].y);

      ++(*num_points);

      return;
    }

    //-----

    double x1, y1, x4, y4;

    bezier.getFirstPoint(&x1, &y1);
    bezier.getLastPoint (&x4, &y4);

    double a =  y4   - y1;
    double b =  x1   - x4;
    double c = -x1*a - y1*b;

    double s = a*a + b*b;

    if (s == 0.0)
      return;

    double x2, y2, x3, y3;

    bezier.getControlPoint1(&x2, &y2);
    bezier.getControlPoint2(&x3, &y3);

    double f2 = a*x2 + b*y2 + c;
    double f3 = a*x3 + b*y3 + c;

    double s2 = f2*f2/s;
    double s3 = f3*f3/s;

    if (! checkLength(s2) || ! checkLength(s3)) {
      C3Bezier2D bezier1, bezier2;

      bezier.split(bezier1, bezier2);

      calc1(bezier1, num_points, depth + 1);
      calc1(bezier2, num_points, depth + 1);
    }
    else {
      bezier.getLastPoint(&bezierPoints_[*num_points].x, &bezierPoints_[*num_points].y);

      ++(*num_points);
    }
  }

  void calc1(const C2Bezier2D &bezier, uint *num_points, uint depth) {
    if (*num_points >= bezierPoints_.size())
      bezierPoints_.resize(2*bezierPoints_.size() + *num_points + 1);

    //-----

    if (depth >= MAX_BEZIER_DEPTH) {
      bezier.getLastPoint(&bezierPoints_[*num_points].x, &bezierPoints_[*num_points].y);

      ++(*num_points);

      return;
    }

    //-----

    double x1, y1, x3, y3;

    bezier.getFirstPoint(&x1, &y1);
    bezier.getLastPoint (&x3, &y3);

    double a =  y3   - y1;
    double b =  x1   - x3;
    double c = -x1*a - y1*b;

    double s = a*a + b*b;

    if (s == 0.0)
      return;

    double x2, y2;

    bezier.getControlPoint(&x2, &y2);

    double f2 = a*x2 + b*y2 + c;

    double s2 = f2*f2/s;

    if (! checkLength(s2)) {
      double x12 = (x1 + x2 )/2.0;
      double y12 = (y1 + y2 )/2.0;

      double x23 = (x2 + x3 )/2.0;
      double y23 = (y2 + y3 )/2.0;

      double xm = (x12 + x23)/2.0;
      double ym = (y12 + y23)/2.0;

      C2Bezier2D bezier1(x1, y1, x12, y12, xm, ym);

      calc1(bezier1, num_points, depth + 1);

      C2Bezier2D bezier2(xm, ym, x23, y23, x3, y3);

      calc1(bezier2, num_points, depth + 1);
    }
    else {
      bezier.getLastPoint(&bezierPoints_[*num_points].x, &bezierPoints_[*num_points].y);

      ++(*num_points);
    }
  }

 protected:
  double    tol_ { -1 };
  PointList bezierPoints_;
};

#endif
