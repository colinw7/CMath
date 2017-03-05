#include <CGrahamHullI.h>
#include <cstdio>
#include <cstdlib>

CGrahamHullI *hull = 0;

CGrahamHullI::
CGrahamHullI() :
 numDeleted_(0), debug_(false)
{
}

void
CGrahamHullI::
addPoint(int x, int y)
{
  Point p;

  p.v.x     = x;
  p.v.y     = y;
  p.vnum    = points_.size();
  p.removed = false;

  points_.push_back(p);
}

bool
CGrahamHullI::
calc()
{
  findLowest();

  sortPoints();

  squash();

  if (! graham())
    return false;

  if (debug_) {
    printf("Hull:\n");

    printStack();
  }

  return true;
}

/*---------------------------------------------------------------------
  findLowest finds the rightmost lowest point and swaps with 0-th.
  The lowest point has the min y-coord, and amongst those, the
  max x-coord: so it is rightmost among the lowest.
  ---------------------------------------------------------------------*/
void
CGrahamHullI::
findLowest()
{
  uint numPoints = points_.size();

  int m = 0; /* Index of lowest so far. */

  for (uint i = 1; i < numPoints; ++i) {
    if ((points_[i].v.y < points_[m].v.y) ||
        ((points_[i].v.y == points_[m].v.y) && (points_[i].v.x > points_[m].v.x)))
      m = i;
  }

  if (debug_)
    printf("Swapping %d with 0\n", m);

  std::swap(points_[0], points_[m]);

  if (debug_)
    printPoints();
}

void
CGrahamHullI::
sortPoints()
{
  hull = this;

  numDeleted_ = 0;

  int numPoints = points_.size();

  // reset removed
  for (int i = 0; i < numPoints; ++i)
    points_[i].removed = false;

  qsort(&points_[1], numPoints - 1, sizeof(Point), compare);

  if (debug_) {
    printf("After sorting, ndelete = %d:\n", numDeleted_);

    printPoints();
  }
}

/*---------------------------------------------------------------------
  squash removes all elements from points_ marked removed.
  ---------------------------------------------------------------------*/
void
CGrahamHullI::
squash()
{
  if (numDeleted_ == 0)
    return;

  int numPoints = points_.size();

  int j = 0;

  for (int i = 0; i < numPoints; ++i) {
    if (! points_[i].removed) { /* if not marked for deletion */
      points_[j].v = points_[i].v;

      j++;
    }
    else {
      /* else do nothing: removed by skipping. */
    }
  }

  while (numPoints > j) {
    points_.pop_back();

    --numPoints;
  }

  if (debug_) {
    printf("After Squash: numPoints=%d\n",numPoints);

    printPoints();
  }
}

/*---------------------------------------------------------------------
  Performs the Graham scan on an array of angularly sorted points points_.
  ---------------------------------------------------------------------*/
bool
CGrahamHullI::
graham()
{
  stack_.clear();

  /* Initialize stack. */
  int i = 0;

  int numPoints = points_.size();

  if (i < numPoints) push(&points_[i++]);
  if (i < numPoints) push(&points_[i++]);

  /* Bottom two elements will never be removed. */
  while (i < numPoints) {
    if (debug_) {
      printf("Stack at top of while loop, i=%d, vnum=%d:\n", i, points_[i].vnum);

      printStack();
    }

    uint num = stack_.size();

    if (num < 2) { printf("Error\n"); return false; }

    Point *p1 = stack_[num - 2];
    Point *p2 = stack_[num - 1];
    Point *p3 = &points_[i];

    if (pointLineLeft(p1->v, p2->v, p3->v)) {
      push(p3);

      ++i;
    }
    else
      pop();

    if (debug_) {
      printf("Stack at bot of while loop, i=%d, vnum=%d:\n", i, points_[i].vnum);

      printStack();

      putchar('\n');
    }
  }

  return true;
}

void
CGrahamHullI::
getHull(std::vector<Coord> &coords) const
{
  uint num = stack_.size();

  for (uint i = 0; i < num; ++i) {
    Point *p = stack_[num - i - 1];

    coords.push_back(p->v);
  }
}

/*---------------------------------------------------------------------
  compare: returns -1,0,+1 if p1 < p2, =, or > respectively;
  here "<" means smaller angle.  Follows the conventions of qsort.
  ---------------------------------------------------------------------*/
int
CGrahamHullI::
compare(const void *tpi, const void *tpj)
{
  Point *p0 = &hull->points_[0];
  Point *pi = (Point *) tpi;
  Point *pj = (Point *) tpj;

  int a = areaSign(hull->points_[0].v, pi->v, pj->v);

  if (a > 0) return -1;
  if (a < 0) return  1;

  /* Collinear with points_[0] */
  int x = abs(pi->v.x - p0->v.x) - abs(pj->v.x - p0->v.x);
  int y = abs(pi->v.y - p0->v.y) - abs(pj->v.y - p0->v.y);

  ++hull->numDeleted_;

  if ((x < 0) || (y < 0)) {
    pi->removed = true;
    return -1;
  }

  if ((x > 0) || (y > 0)) {
    pj->removed = true;
    return 1;
  }

  /* points are coincident */
  if (pi->vnum > pj->vnum)
    pj->removed = true;
  else
    pi->removed = true;

  return 0;
}

/*---------------------------------------------------------------------
  pops off top elment of stack
  ---------------------------------------------------------------------*/
CGrahamHullI::Point *
CGrahamHullI::
pop()
{
  Point *p = stack_.back();

  stack_.pop_back();

  return p;
}

/*---------------------------------------------------------------------
  push new point onto the stack.
  ---------------------------------------------------------------------*/
void
CGrahamHullI::
push(Point *p)
{
  stack_.push_back(p);
}

/*---------------------------------------------------------------------
  Returns true iff c is strictly to the left of the directed line through a to b.
  ---------------------------------------------------------------------*/
bool
CGrahamHullI::
pointLineLeft(const Coord &a, const Coord &b, const Coord &c)
{
  return areaSign(a, b, c) > 0;
}

int
CGrahamHullI::
areaSign(const Coord &a, const Coord &b, const Coord &c)
{
  double area2 = (b.x - a.x) * (double)(c.y - a.y) - (c.x - a.x) * (double)(b.y - a.y);

  //printf("areaSign: area2=%g\n", area2);

  /* The area should be an integer. */
  if      (area2 >  0.5) return  1;
  else if (area2 < -0.5) return -1;
  else                   return  0;
}

//-------------------------

void
CGrahamHullI::
printStack()
{
  uint num = stack_.size();

  if (! num) printf("Empty stack\n");

  for (uint i = 0; i < num; ++i) {
    Point *p = stack_[num - i - 1];

    printf("vnum=%d\tx=%d\ty=%d\n", p->vnum, p->v.x, p->v.y);
  }
}

void
CGrahamHullI::
printPoints()
{
  printf("Points:\n");

  uint numPoints = points_.size();

  for (uint i = 0; i < numPoints; i++)
    printf("vnum=%3d, x=%4d, y=%4d, removed=%d\n",
           points_[i].vnum, points_[i].v.x, points_[i].v.y, points_[i].removed);
}
