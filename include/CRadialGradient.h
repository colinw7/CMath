#ifndef CRADIAL_GRADIENT_H
#define CRADIAL_GRADIENT_H

#include <CGenGradient.h>
#include <CMathGeom2D.h>

class CRadialGradient : public CGenGradient {
 public:
  CRadialGradient() :
   CGenGradient(), cx_(0), cy_(0), r_(1), fx_(0), fy_(0),
   cx1_(0), cy1_(0), r1_(0), fx1_(0), fy1_(0) {
  }

  CRadialGradient(const CRadialGradient &rg) :
   CGenGradient(rg), cx_(rg.cx_), cy_(rg.cy_), r_(rg.r_), fx_(rg.fx_), fy_(rg.fy_),
   cx1_(rg.cx1_), cy1_(rg.cy1_), r1_(rg.r1_), fx1_(rg.fx1_), fy1_(rg.fy1_) {
  }

  CGenGradient *dup() const {
    return new CRadialGradient(*this);;
  }

  double getCenterX() const { return cx_; }
  double getCenterY() const { return cy_; }

  CPoint2D getCenter() const { return CPoint2D(cx_, cy_); }
  void setCenter(const CPoint2D &c) { cx_ = c.x; cy_ = c.y; }

  double getRadius() const { return r_; }

  double getFocusX() const { return fx_; }
  double getFocusY() const { return fy_; }

  CPoint2D getFocus() const { return CPoint2D(fx_, fy_); }
  void setFocus(const CPoint2D &f) { fx_ = f.x; fy_ = f.y; }

  void setCenter(double cx, double cy) { cx_ = cx; cy_ = cy; }

  void setRadius(double r) { r_ = r; }

  void setFocus(double fx, double fy) { fx_ = fx; fy_ = fy; }

  void init(double width, double height) const {
    CRadialGradient *th = const_cast<CRadialGradient *>(this);

    th->cx1_ = cx_*width ;
    th->cy1_ = cy_*height;
    th->r1_  = r_*std::max(width, height);
    th->fx1_ = fx_*width ;
    th->fy1_ = fy_*height;
  }

  double getOffset(double x, double y) const {
    uint   ni;
    double xi1, yi1, xi2, yi2;

    double dx = x - fx1_;
    double dy = y - fy1_;

    double dpf = dx*dx + dy*dy;

    if (! CMathGeom2D::CircleLineIntersect(cx1_, cy1_, r1_, fx1_, fy1_, x, y,
                                           &xi1, &yi1, &xi2, &yi2, &ni))
      return 0.0;

    double mu = 0.0;

    if (! realEq(xi1, fx1_))
      mu = (x - fx1_)/(xi1 - fx1_);
    else
      mu = (y - fy1_)/(yi1 - fy1_);

    if (mu >= 0) {
      dx = xi1 - fx1_;
      dy = yi1 - fy1_;
    }
    else {
      dx = xi2 - fx1_;
      dy = yi2 - fy1_;
    }

    double dif = dx*dx + dy*dy;

    double offset = sqrt(dpf/dif);

    if      (spread_ == CGRADIENT_SPREAD_REPEAT) {
      while (offset < 0.0)
        offset += 1;

      while (offset > 1.0)
        offset -= 1;
    }
    else if (spread_ == CGRADIENT_SPREAD_REFLECT) {
      if (offset < 0.0) {
        uint count = 0;

        while (offset < 0.0) {
          offset += 1;

          ++count;
        }

        if (count & 1)
          offset = 1.0 - offset;
      }

      if (offset > 1.0) {
        uint count = 0;

        while (offset > 1.0) {
          offset -= 1;

          ++count;
        }

        if (count & 1)
          offset = 1.0 - offset;
      }
    }

    return std::min(std::max(offset, 0.0), 1.0);
  }

  CRGBA getColor(double x, double y) const {
    double offset = getOffset(x, y);

    CRGBA  rgba1, rgba2;
    double offset1, offset2;

    StopList::const_iterator p1 = beginStops();
    StopList::const_iterator p2 = endStops  ();

    if (p1 != p2)
      rgba1 = (*p1).getColor();

    offset1 = 0.0;

    for ( ; p1 != p2; ++p1) {
      offset2 = (*p1).getOffset();
      rgba2   = (*p1).getColor();

      if (offset >= offset1 && offset <= offset2) {
        double o21 = offset2 - offset1;

        if (! realEq(offset1, offset2))
          return rgba1*((offset2 - offset)/o21) + rgba2*((offset - offset1)/o21);
        else
          return (rgba1 + rgba2)/2.0;
      }

      offset1 = offset2;
      rgba1   = rgba2;
    }

    offset2 = 1.0;
    rgba2   = rgba1;

    double o21 = offset2 - offset1;

    return rgba1*((offset2 - offset)/o21) + rgba2*((offset - offset1)/o21);
  }

 private:
  static bool realEq(double r1, double r2) {
    return (fabs(r1 - r2) < 1E-5);
  }

 private:
  const CRadialGradient &operator=(const CRadialGradient &rg);

 private:
  double cx_, cy_, r_;
  double fx_, fy_;
  double cx1_, cy1_, r1_;
  double fx1_, fy1_;
};

#endif
