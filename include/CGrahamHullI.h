#ifndef CGRAHAM_HULL_I_H
#define CGRAHAM_HULL_I_H

#include <vector>
#include <sys/types.h>

class CGrahamHullI {
 public:
  struct Coord {
    int x, y;
  };

  struct Point {
    int   vnum;
    Coord v;
    bool  removed;
  };

  typedef std::vector<Coord> CoordArray;

 public:
  CGrahamHullI();

  void setDebug(bool debug) { debug_ = debug; }

  void addPoint(int x, int y);

  bool calc();

  void getHull(CoordArray &coords) const;

 private:
  void findLowest();
  void sortPoints();
  void squash();
  bool graham();

  Point *pop();
  void   push(Point *p);

  static int compare(const void *tp1, const void *tp2);

  static bool pointLineLeft(const Coord &a, const Coord &b, const Coord &c);
  static int  areaSign(const Coord &a, const Coord &b, const Coord &c);

  void printStack();
  void printPoints();

 private:
  typedef std::vector<Point>   PointArray;
  typedef std::vector<Point *> PointStack;

  PointArray points_;
  uint       numDeleted_;
  PointStack stack_;
  bool       debug_;
};

#endif
