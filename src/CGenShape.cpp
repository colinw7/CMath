#include <CGenShape.h>

bool
Boundary::
intersect(const IRect &rect, BoundaryP &iboundary) const
{
  // TODO: optimize
  BoundaryP boundary(new RectBoundary(rect));

  return intersectImpl(boundary.get(), iboundary);
}

bool giftWrap(const std::vector<RPoint> &points, std::vector<uint> &poly_inds);

//------

void
RectBoundary::
moveImpl(int dx, int dy)
{
  rect_.moveBy(dx, dy);
}

void
RectBoundary::
expandByImpl(int dx, int dy)
{
  rect_.expandBy(dx, dy);
}

bool
RectBoundary::
overlapsImpl(const Boundary *boundary) const
{
  return IRect::overlaps(rect_, static_cast<const RectBoundary *>(boundary)->rect_);
}

bool
RectBoundary::
intersectImpl(const Boundary *boundary, BoundaryP &iboundary) const
{
  IRect irect;

  if (! IRect::intersect(rect_, static_cast<const RectBoundary *>(boundary)->rect_, irect))
    return false;

  iboundary = BoundaryP(new RectBoundary(irect));

  return true;
}

bool
RectBoundary::
contains(const IPoint &p) const
{
  return rect_.contains(p);
}

void
RectBoundary::
getCentroid(IPoint &c) const
{
  rect_.getCenter(c);
}

double
RectBoundary::
getArea() const
{
  return getRectArea(rect_);
}

void
RectBoundary::
getRects(std::vector<IRect> &rects) const
{
  rects.push_back(rect_);
}

void
RectBoundary::
calcEnclosure()
{
  enclosure_.push_back(rect_.ll());
  enclosure_.push_back(rect_.lr());
  enclosure_.push_back(rect_.ur());
  enclosure_.push_back(rect_.ul());
}

void
RectBoundary::
print(std::ostream &os) const
{
  os << "{" << rect_.left () << " " << rect_.bottom() << " " <<
               rect_.right() << " " << rect_.top   () << "}";
}

//------

void
RectListBoundary::
moveImpl(int dx, int dy)
{
  uint n = rects_.size();

  for (uint i = 0; i < n; ++i) {
    IRect &rect = rects_[i];

    rect.moveBy(dx, dy);
  }
}

void
RectListBoundary::
expandByImpl(int dx, int dy)
{
  uint n = rects_.size();

  for (uint i = 0; i < n; ++i) {
    IRect &rect = rects_[i];

    rect.expandBy(dx, dy);
  }
}

bool
RectListBoundary::
overlapsImpl(const Boundary *boundary) const
{
  if      (boundary->type() == RECT) {
    const IRect &rect = static_cast<const RectBoundary *>(boundary)->rect_;

    uint n = rects_.size();

    for (uint i = 0; i < n; ++i) {
      const IRect &rect1 = rects_[i];

      if (IRect::overlaps(rect, rect1))
        return true;
    }
  }
  else { // RECT_LIST
    const IRects &rects = static_cast<const RectListBoundary *>(boundary)->rects_;

    uint n1 = rects_.size();
    uint n2 = rects .size();

    for (uint i = 0; i < n1; ++i) {
      const IRect &rect1 = rects_[i];

      for (uint j = 0; j < n2; ++j) {
        const IRect &rect2 = rects[j];

        if (IRect::overlaps(rect1, rect2))
          return true;
      }
    }
  }

  return false;
}

bool
RectListBoundary::
intersectImpl(const Boundary *boundary, BoundaryP &iboundary) const
{
  if      (boundary->type() == RECT) {
    IRects irects;

    const IRect &rect = static_cast<const RectBoundary *>(boundary)->rect_;

    uint n = rects_.size();

    for (uint i = 0; i < n; ++i) {
      const IRect &rect1 = rects_[i];

      IRect irect;

      if (IRect::intersect(rect, rect1, irect))
        irects.push_back(irect);
    }

    if (irects.empty())
      return false;

    iboundary = BoundaryP(new RectListBoundary(irects));

    return true;
  }
  else { // RECT_LIST
    IRects irects;

    const IRects &rects = static_cast<const RectListBoundary *>(boundary)->rects_;

    uint n1 = rects_.size();
    uint n2 = rects .size();

    for (uint i = 0; i < n1; ++i) {
      const IRect &rect1 = rects_[i];

      for (uint j = 0; j < n2; ++j) {
        const IRect &rect2 = rects[j];

        IRect irect;

        if (IRect::intersect(rect1, rect2, irect))
          irects.push_back(irect);
      }
    }

    if (rects.empty())
      return false;

    iboundary = BoundaryP(new RectListBoundary(irects));

    return true;
  }

  return true;
}

bool
RectListBoundary::
contains(const IPoint &p) const
{
  uint n = rects_.size();

  for (uint i = 0; i < n; ++i) {
    const IRect &rect = rects_[i];

    if (rect.contains(p))
      return true;
  }

  return false;
}

void
RectListBoundary::
getCentroid(IPoint &c) const
{
  getRect().getCenter(c);
}

double
RectListBoundary::
getArea() const
{
  double area = 0.0;

  uint n = rects_.size();

  for (uint i = 0; i < n; ++i) {
    const IRect &rect = rects_[i];

    area += getRectArea(rect);
  }

  return area;
}

void
RectListBoundary::
getRects(std::vector<IRect> &rects) const
{
  uint n = rects_.size();

  for (uint i = 0; i < n; ++i) {
    const IRect &rect = rects_[i];

    rects.push_back(rect);
  }
}

void
RectListBoundary::
calcEnclosure()
{
  // calc all points
  std::vector<RPoint> points;

  int nr = rects_.size();

  points.resize(nr*4);

  for (int i = 0, k = 0; i < nr; ++i, k += 4) {
    const IRect &rect = rects_[i];

    std::vector<IPoint> ipoints;

    rect.toPoints(ipoints);

    for (int j = 0; j < 4; ++j)
      points[k + j] = RPoint(ipoints[j]);
  }

  // build convex boundary
  std::vector<uint> inds;

  (void) giftWrap(points, inds);

  // TODO: remove colinear points (polygon orientation == 0)
  uint num_inds = inds.size();

  // store in enclosure
  enclosure_.resize(num_inds);

  for (uint i = 0; i < num_inds; ++i)
    enclosure_[i] = points[inds[i]].ipoint();
}

void
RectListBoundary::
print(std::ostream &os) const
{
  uint n = rects_.size();

  os << "{";

  for (uint i = 0; i < n; ++i) {
    const IRect &rect = rects_[i];

    os << "{" << rect.left () << " " << rect.bottom() << " " <<
                 rect.right() << " " << rect.top   () << "}";
  }

  os << "}";
}

//------

void
PolyBoundary::
moveImpl(int dx, int dy)
{
  poly_.moveBy(dx, dy);
}

void
PolyBoundary::
expandByImpl(int dx, int dy)
{
  poly_.expandBy(dx, dy);
}

bool
PolyBoundary::
overlapsImpl(const Boundary *boundary) const
{
  if      (boundary->type() == RECT) {
    Polygon poly(static_cast<const RectBoundary *>(boundary)->rect_);

    Polygon ipoly;

    if (! poly_.intersect(poly, ipoly))
      return false;

    return true;
  }
  else if (boundary->type() == RECT_LIST) {
    const IRects &rects = static_cast<const RectListBoundary *>(boundary)->rects_;

    Polygons ipolys;

    uint n = rects.size();

    for (uint i = 0; i < n; ++i) {
      const IRect &rect1 = rects[i];

      Polygon ipoly;

      if (poly_.intersect(rect1, ipoly))
        return true;
    }

    return false;
  }
  else { // POLY
    Polygon ipoly;

    if (! poly_.intersect(static_cast<const PolyBoundary *>(boundary)->poly_, ipoly))
      return false;

    return true;
  }
}

bool
PolyBoundary::
intersectImpl(const Boundary *boundary, BoundaryP &iboundary) const
{
  if      (boundary->type() == RECT) {
    Polygon poly(static_cast<const RectBoundary *>(boundary)->rect_);

    Polygon ipoly;

    if (! poly_.intersect(poly, ipoly))
      return false;

    iboundary = BoundaryP(new PolyBoundary(ipoly));
  }
  else if (boundary->type() == RECT_LIST) {
    const IRects &rects = static_cast<const RectListBoundary *>(boundary)->rects_;

    Polygons ipolys;

    uint n = rects.size();

    for (uint i = 0; i < n; ++i) {
      const IRect &rect1 = rects[i];

      Polygon ipoly;

      if (poly_.intersect(rect1, ipoly))
        ipolys.push_back(ipoly);
    }

    if (ipolys.empty())
      return false;

    iboundary = BoundaryP(new PolyListBoundary(ipolys));
  }
  else {
    Polygon ipoly;

    if (! poly_.intersect(static_cast<const PolyBoundary *>(boundary)->poly_, ipoly))
      return false;

    iboundary = BoundaryP(new PolyBoundary(ipoly));
  }

  return true;
}

bool
PolyBoundary::
contains(const IPoint &p) const
{
  return poly_.inside(RPoint(p.x(), p.y()));
}

void
PolyBoundary::
getCentroid(IPoint &c) const
{
  RPoint rc;

  poly_.getCentroid(rc);

  c = rc.ipoint();
}

double
PolyBoundary::
getArea() const
{
  return poly_.getArea();
}

void
PolyBoundary::
getRects(std::vector<IRect> &rects) const
{
  rects.push_back(getRect());
}

void
PolyBoundary::
calcEnclosure()
{
  // get polygon points
  const std::vector<RPoint> &points = poly_.getPoints();

  // build convex boundary
  std::vector<uint> inds;

  (void) giftWrap(points, inds);

  // TODO: remove colinear points (polygon orientation == 0)
  uint num_inds = inds.size();

  // store in enclosure
  enclosure_.resize(num_inds);

  for (uint i = 0; i < num_inds; ++i)
    enclosure_[i] = points[inds[i]].ipoint();
}

void
PolyBoundary::
print(std::ostream &os) const
{
  uint n = poly_.size();

  os << "{";

  for (uint i = 0; i < n; ++i) {
    const RPoint &p = poly_.point(i);

    os << "{" << p.x() << " " << p.y() << "}";
  }

  os << "}";
}

//------

void
PolyListBoundary::
moveImpl(int dx, int dy)
{
  uint n = polys_.size();

  for (uint i = 0; i < n; ++i) {
    Polygon &poly = polys_[i];

    poly.moveBy(dx, dy);
  }
}

void
PolyListBoundary::
expandByImpl(int dx, int dy)
{
  uint n = polys_.size();

  for (uint i = 0; i < n; ++i) {
    Polygon &poly = polys_[i];

    poly.expandBy(dx, dy);
  }
}

bool
PolyListBoundary::
overlapsImpl(const Boundary *boundary) const
{
  if      (boundary->type() == RECT) {
    const IRect &rect = static_cast<const RectBoundary *>(boundary)->rect_;

    uint n = polys_.size();

    for (uint i = 0; i < n; ++i) {
      const Polygon &poly1 = polys_[i];

      Polygon ipoly;

      if (poly1.intersect(rect, ipoly))
        return true;
    }
  }
  else if (boundary->type() == POLY) {
    const Polygon &poly = static_cast<const PolyBoundary *>(boundary)->poly_;

    uint n = polys_.size();

    for (uint i = 0; i < n; ++i) {
      const Polygon &poly1 = polys_[i];

      Polygon ipoly;

      if (poly1.intersect(poly, ipoly))
        return true;
    }
  }
  else if (boundary->type() == RECT_LIST) {
    const IRects &rects = static_cast<const RectListBoundary *>(boundary)->rects_;

    uint n1 = polys_.size();
    uint n2 = rects .size();

    for (uint i = 0; i < n1; ++i) {
      const Polygon &poly1 = polys_[i];

      for (uint j = 0; j < n2; ++j) {
        const IRect &rect2 = rects[j];

        Polygon ipoly;

        if (poly1.intersect(rect2, ipoly))
          return true;
      }
    }
  }
  else {
    const Polygons &polys = static_cast<const PolyListBoundary *>(boundary)->polys_;

    uint n1 = polys_.size();
    uint n2 = polys .size();

    for (uint i = 0; i < n1; ++i) {
      const Polygon &poly1 = polys_[i];

      for (uint j = 0; j < n2; ++j) {
        const Polygon &poly2 = polys[j];

        Polygon ipoly;

        if (poly1.intersect(poly2, ipoly))
          return true;
      }
    }
  }

  return false;
}

bool
PolyListBoundary::
intersectImpl(const Boundary *boundary, BoundaryP &iboundary) const
{
  Polygons ipolys;

  if      (boundary->type() == RECT) {
    const IRect &rect = static_cast<const RectBoundary *>(boundary)->rect_;

    uint n = polys_.size();

    for (uint i = 0; i < n; ++i) {
      const Polygon &poly1 = polys_[i];

      Polygon ipoly;

      if (poly1.intersect(rect, ipoly))
        ipolys.push_back(ipoly);
    }
  }
  else if (boundary->type() == POLY) {
    const Polygon &poly = static_cast<const PolyBoundary *>(boundary)->poly_;

    uint n = polys_.size();

    for (uint i = 0; i < n; ++i) {
      const Polygon &poly1 = polys_[i];

      Polygon ipoly;

      if (poly1.intersect(poly, ipoly))
        ipolys.push_back(ipoly);
    }
  }
  else if (boundary->type() == RECT_LIST) {
    const IRects &rects = static_cast<const RectListBoundary *>(boundary)->rects_;

    uint n1 = polys_.size();
    uint n2 = rects .size();

    for (uint i = 0; i < n1; ++i) {
      const Polygon &poly1 = polys_[i];

      for (uint j = 0; j < n2; ++j) {
        const IRect &rect2 = rects[j];

        Polygon ipoly;

        if (poly1.intersect(rect2, ipoly))
          ipolys.push_back(ipoly);
      }
    }
  }
  else {
    const Polygons &polys = static_cast<const PolyListBoundary *>(boundary)->polys_;

    uint n1 = polys_.size();
    uint n2 = polys .size();

    for (uint i = 0; i < n1; ++i) {
      const Polygon &poly1 = polys_[i];

      for (uint j = 0; j < n2; ++j) {
        const Polygon &poly2 = polys[j];

        Polygon ipoly;

        if (poly1.intersect(poly2, ipoly))
          ipolys.push_back(ipoly);
      }
    }
  }

  if (ipolys.empty())
    return false;

  iboundary = BoundaryP(new PolyListBoundary(ipolys));

  return true;
}

bool
PolyListBoundary::
contains(const IPoint &p) const
{
  RPoint rp(p.x(), p.y());

  uint n = polys_.size();

  for (uint i = 0; i < n; ++i) {
    const Polygon &poly = polys_[i];

    if (poly.inside(rp))
      return true;
  }

  return false;
}

void
PolyListBoundary::
getCentroid(IPoint &c) const
{
  getRect().getCenter(c);
}

double
PolyListBoundary::
getArea() const
{
  double area = 0.0;

  uint n = polys_.size();

  for (uint i = 0; i < n; ++i) {
    const Polygon &poly = polys_[i];

    area += poly.getArea();
  }

  return area;
}

void
PolyListBoundary::
getRects(std::vector<IRect> &rects) const
{
  uint n = polys_.size();

  for (uint i = 0; i < n; ++i) {
    const Polygon &poly = polys_[i];

    rects.push_back(poly.getRect().irect());
  }
}

void
PolyListBoundary::
calcEnclosure()
{
  // get all polygon points
  std::vector<RPoint> points;

  uint n = polys_.size();

  for (uint i = 0; i < n; ++i) {
    const Polygon &poly = polys_[i];

    const std::vector<RPoint> &points1 = poly.getPoints();

    std::copy(points1.begin(), points1.end(), back_inserter(points));
  }

  // build convex boundary
  std::vector<uint> inds;

  (void) giftWrap(points, inds);

  // TODO: remove colinear points (polygon orientation == 0)
  uint num_inds = inds.size();

  // store in enclosure
  enclosure_.resize(num_inds);

  for (uint i = 0; i < num_inds; ++i)
    enclosure_[i] = points[inds[i]].ipoint();
}

void
PolyListBoundary::
print(std::ostream &os) const
{
  uint n = polys_.size();

  os << "{";

  for (uint i = 0; i < n; ++i) {
    const Polygon &poly = polys_[i];

    uint n1 = poly.size();

    os << "{";

    for (uint j = 0; j < n1; ++j) {
      const RPoint &p = poly.point(j);

      os << "{" << p.x() << " " << p.y() << "}";
    }

    os << "}";
  }

  os << "}";
}

//----------------

static uint            gift_wrap_num = 0;
static uint            gift_wrap_max = 0;
static uint           *gift_wrap_ind = 0;
static std::set<uint>  gift_wrap_set;

void
giftWrapAddPoint(uint ind)
{
  if (gift_wrap_num >= gift_wrap_max) {
    uint  max1;
    uint *ind1;

    max1 = 2*gift_wrap_max + 32;
    ind1 = new uint [max1];

    memcpy(ind1, gift_wrap_ind, gift_wrap_num*sizeof(int));

    delete [] gift_wrap_ind;

    gift_wrap_max = max1;
    gift_wrap_ind = ind1;
  }

  gift_wrap_ind[gift_wrap_num++] = ind;

  gift_wrap_set.insert(ind);
}

//! create a convex polygon from an array of points
bool
giftWrap(const std::vector<RPoint> &points, std::vector<uint> &poly_inds)
{
  uint num_points = points.size();

  if (num_points <= 0)
    return false;

  // find lowest, rightmost point
  uint   ind1 = 0;
  double xmin = points[ind1].x();
  double ymin = points[ind1].y();

  for (uint i = 1; i < num_points; ++i) {
    if      (points[i].y() < ymin) {
      ind1 = i;
      xmin = points[ind1].x();
      ymin = points[ind1].y();
    }
    else if (points[i].y() == ymin && points[i].x() > xmin) {
      ind1 = i;
      xmin = points[ind1].x();
      ymin = points[ind1].y();
    }
  }

  gift_wrap_num = 0;

  gift_wrap_set.clear();

  // add first point
  giftWrapAddPoint(ind1);

  // add remaining points
  uint ind2, ind3;

  while (gift_wrap_num != num_points) {
    // process each next free second point
    for (ind2 = 0; ind2 < num_points; ++ind2) {
      // skip, already done
      if (gift_wrap_set.find(ind2) != gift_wrap_set.end()) continue;

      // process each next free third point
      for (ind3 = 0; ind3 < num_points; ++ind3) {
        if (ind3 == ind2 || ind3 == ind1) continue;

        // if three points make non-clockwise shape don't use this point
        int orient = pointOrientation(points[ind1].x(), points[ind1].y(),
                                      points[ind2].x(), points[ind2].y(),
                                      points[ind3].x(), points[ind3].y());

        if (orient < 0)
          break;
      }

      // if all points are clockwise or parallel then second point is good
      if (ind3 >= num_points)
        break;
    }

    // didn't find a point (fail ?)
    if (ind2 >= num_points)
      break;

    // add point to return poly indices
    giftWrapAddPoint(ind2);

    // move to next point
    ind1 = ind2;
  }

  // return poly indices
  poly_inds = std::vector<uint>(gift_wrap_ind, gift_wrap_ind + gift_wrap_num);

  return true;
}
