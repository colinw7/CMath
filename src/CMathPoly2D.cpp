#include <CMathPoly2D.h>
#include <algorithm>

struct IntersectPoint {
  double   mu;
  CPoint2D p;
  int      id;

  IntersectPoint(double mu1, const CPoint2D &p1, int id1) :
   mu(mu1), p(p1), id(id1) {
  }

  bool equal(const IntersectPoint &p) const {
    return (id == p.id);
  }

  friend bool operator<(const IntersectPoint &lhs, const IntersectPoint &rhs) {
    return (lhs.mu < rhs.mu);
  }
};

struct NodePoint {
  enum Type {
    NONE      = 0,
    INSIDE    = (1<<0),
    OUTSIDE   = (1<<1),
    INTERSECT = (1<<2)
  };

  Type       type;
  int        i;
  int        id;
  CPoint2D   p;
  NodePoint *ip;
  bool       used;

  NodePoint(Type type1=NONE, int i1=0, int id1=0, const CPoint2D &p1=CPoint2D()) :
   type(type1), i(i1), id(id1), p(p1), ip(NULL), used(false) {
  }
};

#define FIX_IND(i, n) { if (i >= (n)) i = 0; else if (i < 0) i = (n) - 1; }

#define CALC_ID(i, j, n) ((j)*(n) + (i))

class NodePointArray {
 private:
  typedef std::vector<NodePoint *> NodePoints;

 public:
  NodePointArray() { }

 ~NodePointArray() {
    uint n = size();

    for (uint i = 0; i < n; ++i)
      delete points_[i];
  }

  void push_back(NodePoint *np) {
    points_.push_back(np);
  }

  uint size() const { return points_.size(); }

  const NodePoint *operator[](int i) const { return points_[i]; }

  NodePoint *operator[](int i) { return points_[i]; }

 private:
  NodePoints points_;
};

typedef std::map<int,NodePoint *> NodePointMap;

static bool CPolygonOrOp
             (NodePointArray &points1, NodePointArray &points2, CPolygon2D &rpoly);
static bool CPolygonAndOp
             (NodePointArray &points1, NodePointArray &points2, CPolygon2D &rpoly);
static bool CPolygonXorOp
             (NodePointArray &points1, NodePointArray &points2, CPolygon2D &rpoly);
static bool CPolygonFindPoint
             (NodePointArray &points1, uint type, int &i);

bool
CPolygonBinaryOp(const CPolygon2D &poly1, const CPolygon2D &poly2, CBinaryOp op,
                 std::vector<CPolygon2D> &rpolys)
{
  int n1 = poly1.getNumPoints();
  int n2 = poly2.getNumPoints();

  if (n1 < 3 || n2 < 3) return false;

  //------

  // build list of node points for first polygon
  NodePointArray points1;

  NodePointMap pointMap1;

  {
  // find first point of poly1 not inside poly2
  int ii = 0;

  for ( ; ii < n1; ++ii) {
    const CPoint2D &p = poly1.getPoint(ii);

    if (! poly2.insideEvenOdd(p))
      break;
  }

  // if all inside then we are done
  if (ii >= n1) {
    if      (op == CBINARY_OP_OR)
      rpolys.push_back(poly1);
    else if (op == CBINARY_OP_AND)
      rpolys.push_back(poly1);

    return true;
  }

  //----

  bool inside1 = false;
  bool inside2 = inside1;

  int np = 0;

  for (int i = 0, i1 = ii, i2 = i1 + 1; i < n1; ++i, i1 = i2++, inside1 = inside2) {
    FIX_IND(i2, n1);

    const CPoint2D &p1 = poly1.getPoint(i1);
    const CPoint2D &p2 = poly1.getPoint(i2);

    inside2 = poly2.insideEvenOdd(p2);

    // find intersection points of this edge with other polygon's edges
    std::vector<IntersectPoint> idata;

    for (int j1 = n2 - 1, j2 = 0; j2 < n2; j1 = j2++) {
      // calc intersection (if any)
      const CPoint2D &p21 = poly2.getPoint(j1);
      const CPoint2D &p22 = poly2.getPoint(j2);

      CPoint2D pi;
      double   mu1, mu2;

      if (! CMathGeom2D::IntersectLine(p1, p2, p21, p22, &pi, &mu1, &mu2))
        continue;

      if (mu1 <= 0.0 || mu1 >= 1.0 || mu2 <= 0.0 || mu2 >= 1.0)
        continue;

      // save intersection data
      int id = CALC_ID(i1, j1, n1);

      idata.push_back(IntersectPoint(mu1, pi, id));
    }

    // sort multiple intersects by mu1
    std::sort(idata.begin(), idata.end());

    // add intersect points
    uint n = idata.size();

    for (uint k = 0; k < n; ++k) {
      NodePoint *np1 = new NodePoint(NodePoint::INTERSECT, np++, idata[k].id, idata[k].p);

      points1.push_back(np1);

      pointMap1[idata[k].id] = np1;
    }

    // add end of line
    NodePoint *np1 = new NodePoint((inside2 ? NodePoint::INSIDE : NodePoint::OUTSIDE), np++, 0, p2);

    points1.push_back(np1);
  }
  }

  //-----

  // build list of node points for second polygon
  NodePointArray points2;

  NodePointMap pointMap2;

  {
  // find first point of poly2 not inside poly1
  int ii = 0;

  for ( ; ii < n1; ++ii) {
    const CPoint2D &p = poly2.getPoint(ii);

    if (! poly1.insideEvenOdd(p))
      break;
  }

  // if all inside then we are done
  if (ii >= n1) {
    if      (op == CBINARY_OP_OR)
      rpolys.push_back(poly2);
    else if (op == CBINARY_OP_AND)
      rpolys.push_back(poly2);

    return true;
  }

  //----

  bool inside1 = false;
  bool inside2 = inside1;

  int np = 0;

  for (int i = 0, i1 = ii, i2 = i1 + 1; i < n2; ++i, i1 = i2++, inside1 = inside2) {
    FIX_IND(i2, n2);

    const CPoint2D &p1 = poly2.getPoint(i1);
    const CPoint2D &p2 = poly2.getPoint(i2);

    inside2 = poly1.insideEvenOdd(p2);

    // find intersection points of this edge with other polygon's edges
    std::vector<IntersectPoint> idata;

    for (int j1 = n1 - 1, j2 = 0; j2 < n1; j1 = j2++) {
      int id = CALC_ID(j1, i1, n1);

      // skip if no intersection (faster ?)
      NodePointMap::const_iterator p = pointMap1.find(id);

      if (p == pointMap1.end()) continue;

      // calc intersection
      const CPoint2D &p21 = poly1.getPoint(j1);
      const CPoint2D &p22 = poly1.getPoint(j2);

      CPoint2D pi;
      double   mu1, mu2;

      bool rc = CMathGeom2D::IntersectLine(p1, p2, p21, p22, &pi, &mu1, &mu2);

      assert(rc && mu1 > 0.0 && mu1 < 1.0 && mu2 > 0.0 && mu2 < 1.0);

      idata.push_back(IntersectPoint(mu1, pi, id));
    }

    // sort multiple intersects by mu1
    std::sort(idata.begin(), idata.end());

    // add intersect points
    uint n = idata.size();

    for (uint k = 0; k < n; ++k) {
      NodePoint *np1 = new NodePoint(NodePoint::INTERSECT, np++, idata[k].id, idata[k].p);

      points2.push_back(np1);

      pointMap2[idata[k].id] = np1;
    }

    // add end of line
    NodePoint *np1 = new NodePoint((inside2 ? NodePoint::INSIDE : NodePoint::OUTSIDE), np++, 0, p2);

    points2.push_back(np1);
  }
  }

  //-----

  // Connect Intersections between shapes
  assert(pointMap1.size() == pointMap2.size());

  NodePointMap::iterator pm1, pm2;

  for (pm1 = pointMap1.begin(), pm2 = pointMap1.end(); pm1 != pm2; ++pm1) {
    NodePoint *np1 = (*pm1).second;

    NodePointMap::iterator pm = pointMap2.find(np1->id);

    assert(pm != pointMap2.end());

    NodePoint *np2 = (*pm).second;

    np1->ip = np2;
    np2->ip = np1;
  }

  //-----

  // build list of points for boolean op
  if      (op == CBINARY_OP_OR) {
    CPolygon2D rpoly;

    while (CPolygonOrOp(points1, points2, rpoly))
      rpolys.push_back(rpoly);
  }
  else if (op == CBINARY_OP_AND) {
    CPolygon2D rpoly;

    while (CPolygonAndOp(points1, points2, rpoly))
      rpolys.push_back(rpoly);
  }
  else if (op == CBINARY_OP_XOR) {
    CPolygon2D rpoly;

    while (CPolygonXorOp(points1, points2, rpoly))
      rpolys.push_back(rpoly);
  }
  else if (op == CBINARY_OP_NOT) {
    CPolygon2D rpoly;

    while (CPolygonXorOp(points2, points1, rpoly))
      rpolys.push_back(rpoly);
  }

  return true;
}

static bool
CPolygonOrOp(NodePointArray &points1, NodePointArray &points2, CPolygon2D &rpoly)
{
  rpoly.clearPoints();

  // all outside and intersect
  NodePointArray *points[2];

  points[0] = &points1;
  points[1] = &points2;

  int np1 = points1.size();
  int np2 = points2.size();

  int ind = 0;

  //----

  // find first outside
  int i = 0;

  if (! CPolygonFindPoint(*points[ind], NodePoint::OUTSIDE, i)) {
    ind = 1 - ind;

    if (! CPolygonFindPoint(*points[ind], NodePoint::OUTSIDE, i))
      return false;
  }

  int d = 1;

  while (true) {
    // add outside points
    while ((*points[ind])[i]->type == NodePoint::OUTSIDE && ! (*points[ind])[i]->used) {
      rpoly.addPoint((*points[ind])[i]->p);

      (*points[ind])[i]->used = true;

      i += d;

      if (ind == 0) { FIX_IND(i, np1); }
      else          { FIX_IND(i, np2); }
    }

    if ((*points[ind])[i]->used || (*points[ind])[i]->type != NodePoint::INTERSECT)
      break;

    // add intersect point
    rpoly.addPoint((*points[ind])[i]->p);

    (*points[ind])[i]->used = true;

    if ((*points[ind])[i]->ip)
      (*points[ind])[i]->ip->used = true;

    // check for next outside point
    if (ind == 0) {
      int i1 = i + d                     ; FIX_IND(i1, np1);
      int j1 = (*points[0])[i]->ip->i - 1; FIX_IND(j1, np2);
      int j2 = (*points[0])[i]->ip->i + 1; FIX_IND(j2, np2);

      // switch sides if needed
      if ((*points[0])[i1]->type != NodePoint::OUTSIDE) {
        ind = 1 - ind;

        if      ((*points[1])[j1]->type == NodePoint::OUTSIDE && ! (*points[1])[j1]->used) {
          i = j1;
          d = -1;
        }
        else if ((*points[1])[j2]->type == NodePoint::OUTSIDE && ! (*points[1])[j2]->used) {
          i = j2;
          d = 1;
        }
        else
          break;
      }
    }
    else {
      int j1 = i + d                     ; FIX_IND(j1, np2);
      int i1 = (*points[1])[i]->ip->i - 1; FIX_IND(i1, np1);
      int i2 = (*points[1])[i]->ip->i + 1; FIX_IND(i2, np1);

      // switch sides if needed
      if ((*points[1])[j1]->type != NodePoint::OUTSIDE) {
        ind = 1 - ind;

        if      ((*points[0])[i1]->type == NodePoint::OUTSIDE && ! (*points[0])[i1]->used) {
          i = i1;
          d = -1;
        }
        else if ((*points[0])[i2]->type == NodePoint::OUTSIDE && ! (*points[0])[i2]->used) {
          i = i2;
          d = 1;
        }
        else
          break;
      }
    }
  }

  if (rpoly.getNumPoints() < 3)
    return false;

  return true;
}

static bool
CPolygonAndOp(NodePointArray &points1, NodePointArray &points2, CPolygon2D &rpoly)
{
  rpoly.clearPoints();

  // all inside and intersect
  NodePointArray *points[2];

  points[0] = &points1;
  points[1] = &points2;

  int np1 = points1.size();
  int np2 = points2.size();

  int ind = 0;

  //----

  // find first inside
  int i = 0;

  if (! CPolygonFindPoint(*points[ind], NodePoint::INSIDE, i)) {
    ind = 1 - ind;

    if (! CPolygonFindPoint(*points[ind], NodePoint::INSIDE, i)) {
      // if no inside points then result is all intersection points

      for (int i = 0; i < np1; ++i) {
        if ((*points[0])[i]->type != NodePoint::INTERSECT || (*points[0])[i]->used) continue;

        rpoly.addPoint((*points[0])[i]->p);

        (*points[0])[i]->used = true;
      }

      if (rpoly.getNumPoints() < 3)
        return false;

      return true;
    }
  }

  int d = 1;

  while (true) {
    // add inside points
    while ((*points[ind])[i]->type == NodePoint::INSIDE && ! (*points[ind])[i]->used) {
      rpoly.addPoint((*points[ind])[i]->p);

      (*points[ind])[i]->used = true;

      i += d;

      if (ind == 0) { FIX_IND(i, np1); }
      else          { FIX_IND(i, np2); }
    }

    if ((*points[ind])[i]->used || (*points[ind])[i]->type != NodePoint::INTERSECT)
      break;

    // add intersect point
    rpoly.addPoint((*points[ind])[i]->p);

    (*points[ind])[i]->used = true;

    if ((*points[ind])[i]->ip)
      (*points[ind])[i]->ip->used = true;

    // check for next outside point
    if (ind == 0) {
      int i1 = i + d                     ; FIX_IND(i1, np1);
      int j1 = (*points[0])[i]->ip->i - 1; FIX_IND(j1, np2);
      int j2 = (*points[0])[i]->ip->i + 1; FIX_IND(j2, np2);

      // switch sides if needed
      if ((*points[0])[i1]->type == NodePoint::OUTSIDE) {
        ind = 1 - ind;

        if      ((*points[1])[j1]->type != NodePoint::OUTSIDE && ! (*points[1])[j1]->used) {
          i = j1;
          d = -1;
        }
        else if ((*points[1])[j2]->type != NodePoint::OUTSIDE && ! (*points[1])[j2]->used) {
          i = j2;
          d = 1;
        }
        else
          break;
      }
    }
    else {
      int j1 = i + d                     ; FIX_IND(j1, np2);
      int i1 = (*points[1])[i]->ip->i - 1; FIX_IND(i1, np1);
      int i2 = (*points[1])[i]->ip->i + 1; FIX_IND(i2, np1);

      // switch sides if needed
      if ((*points[1])[j1]->type == NodePoint::OUTSIDE) {
        ind = 1 - ind;

        if      ((*points[0])[i1]->type != NodePoint::OUTSIDE && ! (*points[0])[i1]->used) {
          i = i1;
          d = -1;
        }
        else if ((*points[0])[i2]->type != NodePoint::OUTSIDE && ! (*points[0])[i2]->used) {
          i = i2;
          d = 1;
        }
        else
          break;
      }
    }
  }

  if (rpoly.getNumPoints() < 3)
    return false;

  return true;
}

static bool
CPolygonXorOp(NodePointArray &points1, NodePointArray &points2, CPolygon2D &rpoly)
{
  rpoly.clearPoints();

  // all outside first polygon
  NodePointArray *points[2];

  points[0] = &points1;
  points[1] = &points2;

  int np1 = points1.size();
  int np2 = points2.size();

  //----

  // find first outside
  int i = 0;

  if (! CPolygonFindPoint(*points[0], NodePoint::OUTSIDE, i))
    return false;

  if (i == 0) {
    int j = np1 - 1;

    while (j >= 0 && (*points[0])[j]->type == NodePoint::OUTSIDE && ! (*points[0])[i]->used)
      i = j--;
  }

  while (true) {
    // add outside points
    while ((*points[0])[i]->type == NodePoint::OUTSIDE && ! (*points[0])[i]->used) {
      rpoly.addPoint((*points[0])[i]->p);

      (*points[0])[i]->used = true;

      ++i; FIX_IND(i, np1);
    }

    if ((*points[0])[i]->used || (*points[0])[i]->type != NodePoint::INTERSECT)
      break;

    // add intersect point
    rpoly.addPoint((*points[0])[i]->p);

    (*points[0])[i]->used = true;

    if ((*points[0])[i]->ip)
      (*points[0])[i]->ip->used = true;

    // move to next intersect point
    int j1 = (*points[0])[i]->ip->i - 1; FIX_IND(j1, np2);
    int j2 = (*points[0])[i]->ip->i + 1; FIX_IND(j2, np2);

    if      ((*points[1])[j1]->type != NodePoint::OUTSIDE && ! (*points[1])[j1]->used) {
      int d = -1;

      while ((*points[1])[j1]->type == NodePoint::INSIDE && ! (*points[1])[j1]->used) {
        rpoly.addPoint((*points[1])[j1]->p);

        (*points[1])[j1]->used = true;

        j1 += d; FIX_IND(j1, np2);
      }

      if ((*points[1])[j1]->type != NodePoint::INTERSECT)
        return false;

      i = (*points[1])[j1]->ip->i;
    }
    else if ((*points[1])[j2]->type != NodePoint::OUTSIDE && ! (*points[1])[j2]->used) {
      int d = 1;

      while ((*points[1])[j2]->type == NodePoint::INSIDE && ! (*points[1])[j2]->used) {
        rpoly.addPoint((*points[1])[j2]->p);

        (*points[1])[j2]->used = true;

        j2 += d; FIX_IND(j2, np2);
      }

      if ((*points[1])[j2]->type != NodePoint::INTERSECT)
        return false;

      i = (*points[1])[j2]->ip->i;
    }
    else
      break;
  }

  if (rpoly.getNumPoints() < 3)
    return false;

  return true;
}

bool
CPolygonFindPoint(NodePointArray &points, uint type, int &i)
{
  int np = points.size();

  i = 0;

  while (i < np && (! (points[i]->type & type) || points[i]->used))
    ++i;

  return (i < np);
}
