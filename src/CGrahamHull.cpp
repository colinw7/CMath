#include <CGrahamHull.h>

#include <cstdlib>

static CGrahamHull *s_hull = 0;

CGrahamHull::
CGrahamHull()
{
}

void
CGrahamHull::
addPoint(const CPoint2D &point)
{
  points_.push_back(point);
}

bool
CGrahamHull::
calc()
{
  int num_points = points_.size();

  if (num_points < 3)
    return false;

  // set point indices
  ipoints_.resize(num_points);

  for (int i = 0; i < num_points; ++i)
    ipoints_[i] = i;

  // find lowest point
  findLowest();

  // sort against lowest clockwise
  sortLowestClockwise();

  // remove colinear
  squash();

  int num_ipoints = ipoints_.size();

  if (num_ipoints < 3)
    return false;

  return doScan();
}

bool
CGrahamHull::
doScan()
{
  // do graham scan
  IPoints ipoint_stack;

  ipoint_stack.push_back(ipoints_[0]);
  ipoint_stack.push_back(ipoints_[1]);

  int i = 2;

  int num_ipoints = ipoints_.size();

  while (i < num_ipoints) {
    uint ns = ipoint_stack.size();

    if (ns < 2) {
      ipoints_.clear();
      return false;
    }

    int i1 = ipoint_stack[ns - 2];
    int i2 = ipoint_stack[ns - 1];
    int i3 = ipoints_[i];

    if (pointLineLeft(points_[i1], points_[i2], points_[i3])) {
      ipoint_stack.push_back(i3);

      ++i;
    }
    else
      ipoint_stack.pop_back();
  }

  // set return points
  ipoints_ = ipoint_stack;

  return true;
}

void
CGrahamHull::
getHull(Points &points) const
{
  int num_ipoints = ipoints_.size();

  points.resize(num_ipoints);

  for (int i = 0; i < num_ipoints; ++i)
    points[i] = points_[ipoints_[i]];
}

void
CGrahamHull::
findLowest()
{
  int num_points = ipoints_.size();

  int      min_i = 0;
  CPoint2D min_p = points_[ipoints_[min_i]];

  for (int i = 1; i < num_points; ++i) {
    const CPoint2D &p = points_[ipoints_[i]];

    if (p.y < min_p.y || (p.y == min_p.y && p.x > min_p.x)) {
      min_i = i;
      min_p = p;
    }
  }

  if (min_i > 0)
    std::swap(ipoints_[0], ipoints_[min_i]);
}

void
CGrahamHull::
sortLowestClockwise()
{
  // sort points by angle they make with lowest point and x-axis
  del_points_.clear();

  s_hull = this;

  qsort(&ipoints_[1], ipoints_.size() - 1, sizeof(int), sortLowestClockwiseCmp);
}

int
CGrahamHull::
sortLowestClockwiseCmp(const void *s1, const void *s2)
{
  int i0 = s_hull->ipoints_[0];
  int i1 = *(int *) s1;
  int i2 = *(int *) s2;

  const CPoint2D &p0 = s_hull->points_[i0];
  const CPoint2D &p1 = s_hull->points_[i1];
  const CPoint2D &p2 = s_hull->points_[i2];

  int as = areaSign(p0, p1, p2);

  if (as > 0) return -1;
  if (as < 0) return  1;

  // zero area (p1 and/or p2 are colinear)
  double x = fabs(p1.x - p0.x) - fabs(p2.x - p0.x);
  double y = fabs(p1.y - p0.y) - fabs(p2.y - p0.y);

  // remove colinear points

  if (x < 0 || y < 0) {
    s_hull->del_points_.insert(i1);
    return -1;
  }

  if (x > 0 || y > 0) {
    s_hull->del_points_.insert(i2);
    return 1;
  }

  // p1 and p2 are coincident
  if (i1 > i2)
    s_hull->del_points_.insert(i2);
  else
    s_hull->del_points_.insert(i1);

  return 0;
}

void
CGrahamHull::
squash()
{
  if (del_points_.empty())
    return;

  // remove colinear points
  int num_ipoints = ipoints_.size();

  int j = 0;

  for (int i = 0; i < num_ipoints; ++i) {
    int pi = ipoints_[i];

    if (del_points_.find(pi) == del_points_.end())
      ipoints_[j++] = pi;
  }

  assert(j == (num_ipoints - int(del_points_.size())));

  ipoints_.resize(j);
}

bool
CGrahamHull::
pointLineLeft(const CPoint2D &a, const CPoint2D &b, const CPoint2D &c)
{
  return areaSign(a, b, c) > 0;
}

int
CGrahamHull::
areaSign(const CPoint2D &a, const CPoint2D &b, const CPoint2D &c)
{
  double area2 = (b.x - a.x)*(c.y - a.y) - (c.x - a.x)*(b.y - a.y);

  if      (area2 > 0.0) return  1;
  else if (area2 < 0.0) return -1;
  else                  return  0;
}
