#ifndef CGEN_GRADIENT_H
#define CGEN_GRADIENT_H

#include <vector>
#include <CRGBA.h>

enum CGradientSpreadType {
  CGRADIENT_SPREAD_PAD,
  CGRADIENT_SPREAD_REPEAT,
  CGRADIENT_SPREAD_REFLECT
};

//---

class CGradientStop {
 public:
  CGradientStop(const CRGBA &rgba, double offset) :
   rgba_(rgba), offset_(offset), offset1_(0.0) {
  }

  CGradientStop(const CGradientStop &s) :
   rgba_(s.rgba_), offset_(s.offset_), offset1_(s.offset1_) {
  }

  const CRGBA &getColor() const { return rgba_; }
  void setColor(const CRGBA &rgba) { rgba_ = rgba; }

  double getOffset() const { return offset_; }
  void setOffset(double offset) { offset_ = offset; }

  double getOffset1() const { return offset1_; }
  void setOffset1(double offset) { offset1_ = offset; }

 private:
  CRGBA  rgba_;
  double offset_;
  double offset1_;
};

//---

class CGenGradient {
 public:
  typedef std::vector<CGradientStop> StopList;

 public:
  CGenGradient() :
   stops_(), spread_(CGRADIENT_SPREAD_PAD) {
  }

  CGenGradient(const CGenGradient &g) :
   stops_(g.stops_), spread_(g.spread_) {
  }

  virtual ~CGenGradient() { }

  virtual CGenGradient *dup() const = 0;

  virtual void init(double width, double height) const = 0;

  void setSpread(CGradientSpreadType spread) { spread_ = spread; }

  CGradientSpreadType getSpread() const { return spread_; }

  const StopList &getStops() const { return stops_; }

  void setStops(const StopList &stops) { stops_ = stops; }

  void addStop(double offset, const CRGBA &rgba) {
    addStop(CGradientStop(rgba, offset));
  }

  void addStop(const CRGBA &rgba, double offset) {
    addStop(CGradientStop(rgba, offset));
  }

  void addStop(const CGradientStop &stop) {
    stops_.push_back(stop);
  }

  StopList::const_iterator beginStops() const {
    return stops_.begin();
  }

  StopList::const_iterator endStops() const {
    return stops_.end();
  }

  virtual CRGBA getColor(double x, double y) const = 0;

 private:
  const CGenGradient &operator=(const CGenGradient &g);

 protected:
  StopList            stops_;
  CGradientSpreadType spread_;
};

#endif
