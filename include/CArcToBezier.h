#ifndef CARC_TO_BEZIER_H
#define CARC_TO_BEZIER_H

#include <C3Bezier2D.h>
#include <vector>

class CArcToBezier {
 public:
  typedef std::vector<C3Bezier2D> BezierList;

 public:
  static void ArcToBeziers(double x, double y, double rx, double ry,
                           double angle1, double angle2, BezierList &beziers);

 public:
  CArcToBezier();

  virtual ~CArcToBezier();

  void calc (double x, double y, double xr, double yr, double angle1, double angle2);
  void calcN(double x, double y, double xr, double yr, double angle1, double angle2);

  virtual uint getCalcNumBeziers();

  uint getNumBeziers() const { return arc_beziers_.size(); }

  const C3Bezier2D &getBezier(uint i) const { return arc_beziers_[i]; }

 protected:
  double     x_          { 0 };
  double     y_          { 0 };
  double     xr_         { 1 };
  double     yr_         { 1 };
  double     angle1_     { 0 };
  double     angle2_     { 0 };
  double     angle_diff_ { 0 };
  BezierList arc_beziers_;
};

#endif
