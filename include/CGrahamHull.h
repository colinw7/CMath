#ifndef CGRAHAM_HULL_H
#define CGRAHAM_HULL_H

#include <CPoint2D.h>

#include <vector>
#include <set>

class CGrahamHull {
 public:
  typedef std::vector<CPoint2D> Points;

 public:
  CGrahamHull();

  void addPoint(const CPoint2D &point);

  bool calc();

  void getHull(Points &points) const;

 private:
  void sortLowestClockwise();
  void squash();
  void findLowest();
  bool doScan();

  static int sortLowestClockwiseCmp(const void *tp1, const void *tp2);

  static bool pointLineLeft(const CPoint2D &a, const CPoint2D &b, const CPoint2D &c);

  static int areaSign(const CPoint2D &a, const CPoint2D &b, const CPoint2D &c);

 private:
  typedef std::vector<int> IPoints;
  typedef std::set<int>    DelPoints;

  Points    points_;
  IPoints   ipoints_;
  DelPoints del_points_;
};

#endif
