#ifndef CLINEAR_GRADIENT_H
#define CLINEAR_GRADIENT_H

#include <CGenGradient.h>
#include <CPoint2D.h>

#include <cmath>

class CLinearGradient : public CGenGradient {
 public:
  CLinearGradient() :
   CGenGradient(), x1_(0), y1_(0), x2_(1), y2_(0), c_(0), s_(0), dmin_(0), dmax_(0) {
  }

  CLinearGradient(const CLinearGradient &lg) :
   CGenGradient(lg), x1_(lg.x1_), y1_(lg.y1_), x2_(lg.x2_), y2_(lg.y2_),
   c_(lg.c_), s_(lg.s_), dmin_(lg.dmin_), dmax_(lg.dmax_) {
  }

  CGenGradient *dup() const {
    return new CLinearGradient(*this);
  }

  double getX1() const { return x1_; }
  double getY1() const { return y1_; }
  double getX2() const { return x2_; }
  double getY2() const { return y2_; }

  void setX1(double x1) { x1_ = x1; }
  void setY1(double y1) { y1_ = y1; }
  void setX2(double x2) { x2_ = x2; }
  void setY2(double y2) { y2_ = y2; }

  CPoint2D getPoint1() const { return CPoint2D(x1_, y1_); }
  CPoint2D getPoint2() const { return CPoint2D(x2_, y2_); }

  void setPoint1(const CPoint2D &p) { x1_ = p.x; y1_ = p.y; }
  void setPoint2(const CPoint2D &p) { x2_ = p.x; y2_ = p.y; }

  double getAngle() const {
    return atan2(y2_ - y1_, x2_ - x1_);
  }

  void init(double width, double height) const {
    CLinearGradient *th = const_cast<CLinearGradient *>(this);

    double angle = getAngle();

    th->c_ = cos(angle);
    th->s_ = sin(angle);

    th->dmin_ = distance1(x1_*width, y1_*height);
    th->dmax_ = distance1(x2_*width, y2_*height);

    StopList::const_iterator p1 = beginStops();
    StopList::const_iterator p2 = endStops  ();

    for ( ; p1 != p2; ++p1) {
      CGradientStop *stop = const_cast<CGradientStop *>(&(*p1));

      stop->setOffset1(stop->getOffset()*(dmax_ - dmin_) + dmin_);
    }
  }

  void setAngle(double angle) {
    double c = cos((M_PI*angle)/180.0);
    double s = sin((M_PI*angle)/180.0);

    setLine(0, 0, fabs(c), fabs(s));
  }

  void setLine(double x1, double y1, double x2, double y2) {
    x1_ = x1; y1_ = y1;
    x2_ = x2; y2_ = y2;
  }

  CRGBA getColor(double x, double y) const {
    double d = distance(x, y);

    CRGBA rgba1, rgba2;

    double d1 = 0.0;
    double d2 = 0.0;

    StopList::const_iterator p1 = beginStops();
    StopList::const_iterator p2 = endStops  ();

    if (p1 != p2)
      rgba1 = (*p1).getColor();

    for ( ; p1 != p2; ++p1) {
      d2    = (*p1).getOffset1();
      rgba2 = (*p1).getColor();

      if (d >= d1 && d <= d2) {
        double d21 = d2 - d1;

        return rgba1*((d2 - d)/d21) + rgba2*((d - d1)/d21);
      }

      d1    = d2;
      rgba1 = rgba2;
    }

    d2    = 1.0;
    rgba2 = rgba1;

    double d21 = d2 - d1;

    return rgba1*((d2 - d)/d21) + rgba2*((d - d1)/d21);
  }

  double distance(double x, double y) const {
    double d = distance1(x, y);

    if      (spread_ == CGRADIENT_SPREAD_REPEAT) {
      double d1 = dmax_ - dmin_;

      while (d < dmin_)
        d += d1;

      while (d > dmax_)
        d -= d1;
    }
    else if (spread_ == CGRADIENT_SPREAD_REFLECT) {
      double d1 = dmax_ - dmin_;

      if (d < dmin_) {
        uint count = 0;

        while (d < dmin_) {
          d += d1;

          ++count;
        }

        if (count & 1)
          d = dmin_ + (dmax_ - d);
      }

      if (d > dmax_) {
        uint count = 0;

        while (d > dmax_) {
          d -= d1;

          ++count;
        }

        if (count & 1)
          d = dmin_ + (dmax_ - d);
      }
    }

    return d;
  }

 private:
  double distance1(double x, double y) const { return x*c_ + y*s_; }

  const CLinearGradient &operator=(const CLinearGradient &lg);

 private:
  double x1_, y1_, x2_, y2_;
  double c_, s_;
  double dmin_, dmax_;
};

#endif
