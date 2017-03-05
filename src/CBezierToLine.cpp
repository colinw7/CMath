#include <CBezierToLine.h>

void
CBezierToLine::
calc(const C3Bezier2D &bezier)
{
  init(bezier);

  //----

  if (bezierPoints_.size() <= maxPoints_)
    bezierPoints_.resize(maxPoints_ + 1);

  uint numPoints = 0;

  bezier.getFirstPoint(&bezierPoints_[0].x, &bezierPoints_[0].y);

  ++numPoints;

  calc1(bezier, &numPoints, 1);

  //----

  bezierPoints_.resize(numPoints);
}

void
CBezierToLine::
calc1(const C3Bezier2D &bezier, uint *numPoints, uint depth)
{
  if (*numPoints >= bezierPoints_.size())
    bezierPoints_.resize(2*bezierPoints_.size() + *numPoints + 1);

  //-----

  if (depth >= maxPoints_) {
    bezier.getLastPoint(&bezierPoints_[*numPoints].x, &bezierPoints_[*numPoints].y);

    ++(*numPoints);

    return;
  }

  //-----

  // get end point distance
  double x1, y1, x4, y4;

  bezier.getFirstPoint(&x1, &y1);
  bezier.getLastPoint (&x4, &y4);

  double a = y4 - y1;
  double b = x1 - x4;

  double s = a*a + b*b;

  if (s == 0.0)
    return;

  // get distance of control points to line
  double x2, y2, x3, y3;

  bezier.getControlPoint1(&x2, &y2);
  bezier.getControlPoint2(&x3, &y3);

  double c = -x1*a - y1*b;

  double f2 = a*x2 + b*y2 + c;
  double f3 = a*x3 + b*y3 + c;

  double s2 = f2*f2/s;
  double s3 = f3*f3/s;

  if (! checkLength(s2) || ! checkLength(s3)) {
    C3Bezier2D bezier1, bezier2;

    bezier.split(bezier1, bezier2);

    calc1(bezier1, numPoints, depth + 1);
    calc1(bezier2, numPoints, depth + 1);
  }
  else {
    bezier.getLastPoint(&bezierPoints_[*numPoints].x, &bezierPoints_[*numPoints].y);

    ++(*numPoints);
  }
}

void
CBezierToLine::
calc(const C2Bezier2D &bezier)
{
  init(bezier);

  //----

  if (bezierPoints_.size() <= maxPoints_)
    bezierPoints_.resize(maxPoints_ + 1);

  uint numPoints = 0;

  bezier.getFirstPoint(&bezierPoints_[0].x, &bezierPoints_[0].y);

  ++numPoints;

  calc1(bezier, &numPoints, 1);

  //----

  bezierPoints_.resize(numPoints);
}

void
CBezierToLine::
calc1(const C2Bezier2D &bezier, uint *numPoints, uint depth)
{
  if (*numPoints >= bezierPoints_.size())
    bezierPoints_.resize(2*bezierPoints_.size() + *numPoints + 1);

  //-----

  if (depth >= maxPoints_) {
    bezier.getLastPoint(&bezierPoints_[*numPoints].x, &bezierPoints_[*numPoints].y);

    ++(*numPoints);

    return;
  }

  //-----

  // get end point distance
  double x1, y1, x3, y3;

  bezier.getFirstPoint(&x1, &y1);
  bezier.getLastPoint (&x3, &y3);

  double a = y3 - y1;
  double b = x1 - x3;

  double s = a*a + b*b;

  if (s == 0.0)
    return;

  // get distance of control point
  double x2, y2;

  bezier.getControlPoint(&x2, &y2);

  double c = -x1*a - y1*b;

  double f2 = a*x2 + b*y2 + c;

  double s2 = f2*f2/s;

  if (! checkLength(s2)) {
    C2Bezier2D bezier1, bezier2;

    bezier.split(bezier1, bezier2);

    calc1(bezier1, numPoints, depth + 1);
    calc1(bezier2, numPoints, depth + 1);
  }
  else {
    bezier.getLastPoint(&bezierPoints_[*numPoints].x, &bezierPoints_[*numPoints].y);

    ++(*numPoints);
  }
}

void
CBezierToLine::
init(const C3Bezier2D &bezier)
{
  double x1, y1, x4, y4;

  bezier.getFirstPoint(&x1, &y1);
  bezier.getLastPoint (&x4, &y4);

  if (tol_ < 0.0)
    tol_ = std::max(fabs(x4 - x1), fabs(y4 - y1))/20;
}

void
CBezierToLine::
init(const C2Bezier2D &bezier)
{
  double x1, y1, x3, y3;

  bezier.getFirstPoint(&x1, &y1);
  bezier.getLastPoint (&x3, &y3);

  if (tol_ < 0.0)
    tol_ = std::max(fabs(x3 - x1), fabs(y3 - y1))/20;
}

bool
CBezierToLine::
checkLength(double s)
{
  return (s <= tol_);
}
