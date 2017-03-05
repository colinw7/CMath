#ifndef CGEN_SHAPE_H
#define CGEN_SHAPE_H

#include <boost/shared_ptr.hpp>
#include <iostream>

// point class (integer values)
class IPoint {
 public:
  IPoint(int x=0, int y=0) :
   x_(x), y_(y) {
  }

  int x() const { return x_; }
  int y() const { return y_; }

  void setX(int x) { x_ = x; }
  void setY(int y) { y_ = y; }

  void incX(int dx) { x_ += dx; }
  void incY(int dy) { y_ += dy; }

  void decX(int dx) { x_ -= dx; }
  void decY(int dy) { y_ -= dy; }

  void setXY(int x, int y) { x_ = x; y_ = y; }

  void move(int dx, int dy) { x_ += dx; y_ += dy; }

  int diffX(const IPoint &p) const { return x_ - p.x_; }
  int diffY(const IPoint &p) const { return y_ - p.y_; }

  uint distX(const IPoint &p) const { return abs(diffX(p)); }
  uint distY(const IPoint &p) const { return abs(diffY(p)); }

  uint manhattenDist(const IPoint &p) const { return (distX(p) + distY(p)); }

  IPoint &operator+=(const IPoint &rhs) {
    x_ += rhs.x_; y_ += rhs.y_;

    return *this;
  }

  friend IPoint operator+(const IPoint &lhs, const IPoint &rhs) {
    return IPoint(lhs.x_ + rhs.x_, lhs.y_ + rhs.y_);
  }

  IPoint &operator-=(const IPoint &rhs) {
    x_ -= rhs.x_; y_ -= rhs.y_;

    return *this;
  }

  friend IPoint operator-(const IPoint &lhs, const IPoint &rhs) {
    return IPoint(lhs.x_ - rhs.x_, lhs.y_ - rhs.y_);
  }

  friend bool operator<(const IPoint &lhs, const IPoint &rhs) {
    return (lhs.x_ < rhs.x_ || (lhs.x_ == rhs.x_ && lhs.y_ < rhs.y_));
  }

  friend bool operator>(const IPoint &lhs, const IPoint &rhs) {
    return (lhs.x_ > rhs.x_ || (lhs.x_ == rhs.x_ && lhs.y_ > rhs.y_));
  }

  friend bool operator==(const IPoint &lhs, const IPoint &rhs) {
    return (lhs.x_ == rhs.x_ && lhs.y_ == rhs.y_);
  }

  friend bool operator!=(const IPoint &lhs, const IPoint &rhs) {
    return (lhs.x_ != rhs.x_ || lhs.y_ != rhs.y_);
  }

  void print(std::ostream &os) const {
    os << "(" << x_ << "," << y_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const IPoint &p) {
    p.print(os);

    return os;
  }

 private:
  int x_;
  int y_;
};

// point class (real values)
class RPoint {
 public:
  RPoint(double x=0.0, double y=0.0) :
   x_(x), y_(y) {
  }

  RPoint(const IPoint &ip) :
   x_(ip.x()), y_(ip.y()) {
  }

  double x() const { return x_; }
  double y() const { return y_; }

  void setX(double x) { x_ = x; }
  void setY(double y) { y_ = y; }

  void incX(double dx) { x_ += dx; }
  void incY(double dy) { y_ += dy; }

  void decX(double dx) { x_ -= dx; }
  void decY(double dy) { y_ -= dy; }

  void setXY(double x, double y) { x_ = x; y_ = y; }

  IPoint ipoint() const { return IPoint(round(x_), round(y_)); }

  RPoint &operator+=(const RPoint &rhs) {
    x_ += rhs.x_; y_ += rhs.y_;

    return *this;
  }

  friend RPoint operator+(const RPoint &lhs, const RPoint &rhs) {
    return RPoint(lhs.x_ + rhs.x_, lhs.y_ + rhs.y_);
  }

  RPoint &operator-=(const RPoint &rhs) {
    x_ -= rhs.x_; y_ -= rhs.y_;

    return *this;
  }

  friend RPoint operator-(const RPoint &lhs, const RPoint &rhs) {
    return RPoint(lhs.x_ - rhs.x_, lhs.y_ - rhs.y_);
  }

  friend RPoint operator*(const RPoint &lhs, double rhs) {
    return RPoint(lhs.x_ * rhs, lhs.y_ * rhs);
  }

  friend RPoint operator*(double lhs, const RPoint &rhs) {
    return RPoint(rhs.x_ * lhs, rhs.y_ * lhs);
  }

  RPoint &operator*=(double s) {
    x_ *= s; y_ *= s;

    return *this;
  }

  friend RPoint operator/(const RPoint &lhs, double rhs) {
    return RPoint(lhs.x_ / rhs, lhs.y_ / rhs);
  }

  RPoint &operator/=(double s) {
    x_ /= s; y_ /= s;

    return *this;
  }

  void print(std::ostream &os) const {
    os << "(" << x_ << "," << y_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const RPoint &p) {
    p.print(os);

    return os;
  }

 private:
  double x_;
  double y_;
};

// rectangle class (integer values)
class IRect {
 public:
  IRect() :
   set_(false), xmin_(0), ymin_(0), xmax_(0), ymax_(0) {
  }

  IRect(const IRect &rect) :
   set_(rect.set_), xmin_(rect.xmin_), ymin_(rect.ymin_), xmax_(rect.xmax_), ymax_(rect.ymax_) {
  }

  IRect(int xmin, int ymin, int xmax, int ymax) :
   set_(true), xmin_(xmin), ymin_(ymin), xmax_(xmax), ymax_(ymax) {
  }

  IRect(const IPoint &ll, const IPoint &ur) :
   set_(true), xmin_(ll.x()), ymin_(ll.y()), xmax_(ur.x()), ymax_(ur.y()) {
  }

  bool isSet() const { return set_; }

  void reset() { set_ = false; }

  int left  () const { return xmin_; }
  int right () const { return xmax_; }
  int bottom() const { return ymin_; }
  int top   () const { return ymax_; }

  void setLeft  (int l) { xmin_ = l; set_ = true; }
  void setRight (int r) { xmax_ = r; set_ = true; }
  void setBottom(int b) { ymin_ = b; set_ = true; }
  void setTop   (int t) { ymax_ = t; set_ = true; }

  void offsetLeft  (int dl) { assert(set_); xmin_ += dl; }
  void offsetRight (int dr) { assert(set_); xmax_ += dr; }
  void offsetBottom(int db) { assert(set_); ymin_ += db; }
  void offsetTop   (int dt) { assert(set_); ymax_ += dt; }

  int width () const { return abs(xmax_ - xmin_); }
  int height() const { return abs(ymax_ - ymin_); }

  int midX() const { return avg(xmin_, xmax_); }
  int midY() const { return avg(ymin_, ymax_); }

  IPoint ll() const { return IPoint(xmin_, ymin_); }
  IPoint lr() const { return IPoint(xmax_, ymin_); }
  IPoint ul() const { return IPoint(xmin_, ymax_); }
  IPoint ur() const { return IPoint(xmax_, ymax_); }

  IPoint center() const { IPoint c; c.setXY(midX(), midY()); return c; }

  void getCenter(IPoint &c) const { c.setXY(midX(), midY()); }

  void toPoints(std::vector<IPoint> &points) const {
    points.resize(4);

    points[0] = IPoint(left (), bottom());
    points[2] = IPoint(right(), top   ());
    points[1] = IPoint(points[2].x(), points[0].y());
    points[3] = IPoint(points[0].x(), points[2].y());
  }

  void moveBy(int dx, int dy) {
    assert(set_);

    xmin_ += dx; ymin_ += dy;
    xmax_ += dx; ymax_ += dy;
  }

  void expandBy(int dx, int dy) {
    assert(set_);

    xmin_ -= dx; ymin_ -= dy;
    xmax_ += dx; ymax_ += dy;
  }

  void expandBy(int dxl, int dyb, int dxr, int dyt) {
    assert(set_);

    xmin_ -= dxl; ymin_ -= dyb;
    xmax_ += dxr; ymax_ += dyt;
  }

  // is point (p) inside this rectangle
  bool contains(const IPoint &p) const {
    assert(set_);

    return (p.x() >= xmin_ && p.x() <= xmax_ && p.y() >= ymin_ && p.y() <= ymax_);
  }

  // is x value inside the bounds of this rectangle
  bool containsX(int x) const {
    assert(set_);

    return (x >= xmin_ && x <= xmax_);
  }

  // is y value inside the bounds of this rectangle
  bool containsY(int y) const {
    assert(set_);

    return (y >= ymin_ && y <= ymax_);
  }

  // is rectangle (ir) inside this rectangle
  bool contains(const IRect &ir) const {
    assert(set_ && ir.set_);

    return (ir.xmin_ >= xmin_ && ir.xmax_ <= xmax_ &&
            ir.ymin_ >= ymin_ && ir.ymax_ <= ymax_);
  }

  IRect &combine(const IPoint &p) {
    if (! set_) {
      xmin_ = p.x(); xmax_ = xmin_;
      ymin_ = p.y(); ymax_ = ymin_;
      set_  = true;
    }
    else {
      xmin_ = std::min(xmin_, p.x()); ymin_ = std::min(ymin_, p.y());
      xmax_ = std::max(xmax_, p.x()); ymax_ = std::max(ymax_, p.y());
    }

    return *this;
  }

  IRect &combine(const IRect &r) {
    if (! set_) {
      xmin_ = r.xmin_; xmax_ = r.xmax_;
      ymin_ = r.ymin_; ymax_ = r.ymax_;
      set_  = r.set_;
    }
    else {
      assert(r.set_);

      xmin_ = std::min(xmin_, r.xmin_); ymin_ = std::min(ymin_, r.ymin_);
      xmax_ = std::max(xmax_, r.xmax_); ymax_ = std::max(ymax_, r.ymax_);
    }

    return *this;
  }

  static IRect combine(const IRect &rect1, const IRect &rect2) {
    assert(rect1.set_);
    assert(rect2.set_);

    int l = std::min(rect1.xmin_, rect2.xmin_); int b = std::min(rect1.ymin_, rect2.ymin_);
    int r = std::max(rect1.xmax_, rect2.xmax_); int t = std::max(rect1.ymax_, rect2.ymax_);

    return IRect(l, b, r, t);
  }

  bool overlaps(const IRect &rect) const {
    return overlaps(*this, rect);
  }

  static bool overlaps(const IRect &rect1, const IRect &rect2) {
    assert(rect1.set_ && rect2.set_);

    if (rect1.xmax_ < rect2.xmin_ || rect1.xmin_ > rect2.xmax_ ||
        rect1.ymax_ < rect2.ymin_ || rect1.ymin_ > rect2.ymax_)
      return false;

    return true;
  }

  static bool intersect(const IRect &rect1, const IRect &rect2, IRect &irect) {
    assert(rect1.set_ && rect2.set_);

    if (rect1.xmax_ < rect2.xmin_ || rect1.xmin_ > rect2.xmax_ ||
        rect1.ymax_ < rect2.ymin_ || rect1.ymin_ > rect2.ymax_)
      return false;

    irect.set_  = true;
    irect.xmin_ = std::max(rect1.xmin_, rect2.xmin_);
    irect.ymin_ = std::max(rect1.ymin_, rect2.ymin_);
    irect.xmax_ = std::min(rect1.xmax_, rect2.xmax_);
    irect.ymax_ = std::min(rect1.ymax_, rect2.ymax_);

    return true;
  }

  static IRect intersect(const IRect &rect1, const IRect &rect2) {
    IRect irect;

    bool rc = intersect(rect1, rect2, irect);
    assert(rc);

    return irect;
  }

  friend bool operator<(const IRect &lhs, const IRect &rhs) {
    if (lhs.xmin_ < rhs.xmin_) return true ; if (lhs.xmin_ > rhs.xmin_) return false;
    if (lhs.ymin_ < rhs.ymin_) return true ; if (lhs.ymin_ > rhs.ymin_) return false;
    if (lhs.xmax_ < rhs.xmax_) return true ; if (lhs.xmax_ > rhs.xmax_) return false;
    if (lhs.ymax_ < rhs.ymax_) return true ; if (lhs.ymax_ > rhs.ymax_) return false;

    return false;
  }

  friend bool operator==(const IRect &lhs, const IRect &rhs) {
    return (lhs.xmin_ == rhs.xmin_ && lhs.ymin_ == rhs.ymin_ &&
            lhs.xmax_ == rhs.xmax_ && lhs.ymax_ == rhs.ymax_);
  }

  bool isValid() const { return (xmin_ <= xmax_ && ymin_ <= ymax_); }

  void print(std::ostream &os) const {
    os << "(" << xmin_ << "," << ymin_ << "," << xmax_ << "," << ymax_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const IRect &r) {
    r.print(os);

    return os;
  }

 private:
  void fixup() {
    if (xmin_ > xmax_) std::swap(xmin_, xmax_);
    if (ymin_ > ymax_) std::swap(ymin_, ymax_);
  }

 private:
  bool set_;
  int  xmin_, ymin_;
  int  xmax_, ymax_;
};

// rectangle class (real values)
class RRect {
 public:
  RRect() :
   set_(false), xmin_(0), ymin_(0), xmax_(0), ymax_(0) {
  }

  RRect(const RRect &rect) :
   set_(rect.set_), xmin_(rect.xmin_), ymin_(rect.ymin_), xmax_(rect.xmax_), ymax_(rect.ymax_) {
  }

  RRect(double xmin, double ymin, double xmax, double ymax) :
   set_(true), xmin_(xmin), ymin_(ymin), xmax_(xmax), ymax_(ymax) {
  }

  bool isSet() const { return set_; }

  void reset() { set_ = false; }

  double left  () const { return xmin_; }
  double right () const { return xmax_; }
  double bottom() const { return ymin_; }
  double top   () const { return ymax_; }

  double width () const { return abs(xmax_ - xmin_); }
  double height() const { return abs(ymax_ - ymin_); }

  double midX() const { return avg(xmin_, xmax_); }
  double midY() const { return avg(ymin_, ymax_); }

  void combine(const RPoint &p) {
    if (! set_) {
      xmin_ = p.x(); xmax_ = xmin_;
      ymin_ = p.y(); ymax_ = ymin_;

      set_ = true;
    }
    else {
      xmin_ = std::min(xmin_, p.x());
      ymin_ = std::min(ymin_, p.y());
      xmax_ = std::max(xmax_, p.x());
      ymax_ = std::max(ymax_, p.y());
    }
  }

  void combine(const RRect &r) {
    if (! set_) {
      xmin_ = r.xmin_; xmax_ = r.xmax_;
      ymin_ = r.ymin_; ymax_ = r.ymax_;

      set_ = true;
    }
    else {
      xmin_ = std::min(xmin_, r.xmin_);
      ymin_ = std::min(ymin_, r.ymin_);
      xmax_ = std::max(xmax_, r.xmax_);
      ymax_ = std::max(ymax_, r.ymax_);
    }
  }

  // is point (p) inside this rectangle
  bool contains(const RPoint &p) const {
    assert(set_);

    return (p.x() >= xmin_ && p.x() <= xmax_ &&
            p.y() >= ymin_ && p.y() <= ymax_);
  }

  IRect irect() const {
    return IRect(floor(xmin_), floor(ymin_), ceil(xmax_), ceil(ymax_));
  }

  bool isValid() const { return (xmin_ <= xmax_ && ymin_ <= ymax_); }

  void print(std::ostream &os) const {
    os << "(" << xmin_ << "," << ymin_ << xmax_ << "," << ymax_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const RRect &r) {
    r.print(os);

    return os;
  }

 private:
  bool   set_;
  double xmin_, ymin_;
  double xmax_, ymax_;
};

//-------------------

// class for set of integer rectangles
class IRectSet {
 public:
  typedef std::set<IRect> Rects;

 public:
  explicit IRectSet(const Rects &rects=Rects()) :
   rects_(rects), rect_() {
    Rects::const_iterator p1, p2;

    for (p1 = rects_.begin(), p2 = rects_.end(); p1 != p2; ++p1)
      rect_.combine(*p1);
  }

  void clear() {
    rects_.clear();

    rect_.reset();
  }

  void addRect(const IRect &rect) {
    if (! rect.isSet()) return;

    rects_.insert(rect);

    rect_.combine(rect);
  }

  void addRects(const IRectSet &rects) {
    if (! rects.isSet()) return;

    Rects::const_iterator p1, p2;

    for (p1 = rects.rects_.begin(), p2 = rects.rects_.end(); p1 != p2; ++p1)
      rects_.insert(*p1);

    rect_.combine(rects.rect_);
  }

  bool empty() const { return rects_.empty(); }

  uint size() const { return rects_.size(); }

  bool isSet() const { return rect_.isSet(); }

  const IRect &getRect() const {
    return rect_;
  }

  // is point (p) inside any of the rectangles
  bool contains(const IPoint &p) const {
    if (! rect_.contains(p))
      return false;

    Rects::const_iterator p1, p2;

    for (p1 = rects_.begin(), p2 = rects_.end(); p1 != p2; ++p1)
      if ((*p1).contains(p))
        return true;

      return false;
  }

  bool overlaps(const IRect &rect) const {
    if (! rect.overlaps(rect_))
      return false;

    Rects::const_iterator p1, p2;

    for (p1 = rects_.begin(), p2 = rects_.end(); p1 != p2; ++p1) {
      if (rect.overlaps(*p1))
        return true;
    }

    return false;
  }

  void expandBy(int dx, int dy) {
    rect_.expandBy(dx, dy);

    Rects::iterator p1, p2;

    for (p1 = rects_.begin(), p2 = rects_.end(); p1 != p2; ++p1) {
      IRect &rect = const_cast<IRect &>(*p1);

      rect.expandBy(dx, dy);
    }
  }

  void expandBy(int dxl, int dyb, int dxr, int dyt) {
    rect_.expandBy(dxl, dyb, dxr, dyt);

    Rects::iterator p1, p2;

    for (p1 = rects_.begin(), p2 = rects_.end(); p1 != p2; ++p1) {
      IRect &rect = const_cast<IRect &>(*p1);

      rect.expandBy(dxl, dyb, dxr, dyt);
    }
  }

  Rects::iterator begin() { return rects_.begin(); }
  Rects::iterator end  () { return rects_.end  (); }

  Rects::const_iterator begin() const { return rects_.begin(); }
  Rects::const_iterator end  () const { return rects_.end  (); }

 private:
  Rects rects_;
  IRect rect_;
};

//-------------------

// adaptor for integer rectangle to templated method interface
template<>
struct RectShimT<IRect> {
  static int width (const IRect &r) { return r.width (); }
  static int height(const IRect &r) { return r.height(); }

  static int midX(const IRect &r) { return r.midX(); }
  static int midY(const IRect &r) { return r.midY(); }

  static int left  (const IRect &r) { return r.left  (); }
  static int bottom(const IRect &r) { return r.bottom(); }
  static int right (const IRect &r) { return r.right (); }
  static int top   (const IRect &r) { return r.top   (); }

  static IPoint ll(const IRect &r) { return r.ll(); }
  static IPoint lr(const IRect &r) { return r.lr(); }
  static IPoint ul(const IRect &r) { return r.ul(); }
  static IPoint ur(const IRect &r) { return r.ur(); }

  static bool overlaps(const IRect &r1, const IRect &r2) {
    return rectsOverlap(r1, r2);
  }

  static IRect rotate45(const IRect &r) {
    return rotateRect45<IRect,IPoint>(r);
  }
};

// adaptor for integer poin to templated method interface
template<>
struct PointShimT<IPoint> {
  static int x(const IPoint &p) { return p.x(); }
  static int y(const IPoint &p) { return p.y(); }

  static IPoint rotate45(const IPoint &p) {
    return rotatePoint45(p);
  }
};

//-------------------

// polygon class (real values)
class Polygon {
 public:
  typedef std::vector<RPoint> Points;

 public:
  Polygon(const Points &points=Points()) :
   points_(points) {
  }

  Polygon(const IRect &rect) {
    points_.push_back(RPoint(rect.left (), rect.bottom()));
    points_.push_back(RPoint(rect.right(), rect.bottom()));
    points_.push_back(RPoint(rect.right(), rect.top   ()));
    points_.push_back(RPoint(rect.left (), rect.top   ()));
  }

  uint size() const { return points_.size(); }

  const RPoint &point(uint i) const { return points_[i]; }

  const Points &getPoints() const { return points_; }

  void addPoint(const RPoint &p) {
    points_.push_back(p);

    bbox_.reset();
  }

  void setPoints(const Points &points) {
    points_ = points;

    bbox_.reset();
  }

  const RRect &getRect() const {
    if (! bbox_.isSet()) {
      uint n = points_.size();

      for (uint i = 0; i < n; ++i)
        bbox_.combine(points_[i]);
    }

    return bbox_;
  }

  void getCentroid(RPoint &c) const {
    double xc = 0.0;
    double yc = 0.0;

    double area = 0.0;

    int n = points_.size();

    int i = n - 1;

    for (int j = 0; j < n; i = j++) {
      const RPoint &p1 = points_[i];
      const RPoint &p2 = points_[j];

      double xy = (p1.x()*p2.y()) - (p1.y()*p2.x());

      xc += (p1.x() + p2.x())*xy;
      yc += (p1.y() + p2.y())*xy;

      area += xy;
    }

    double f = 3.0*area;

    xc /= f;
    yc /= f;

    c = RPoint(xc, yc);
  }

  double getArea() const {
    double area = 0.0;

    int n = points_.size();

    int i = n - 1;

    for (int j = 0; j < n; i = j++) {
      const RPoint &p1 = points_[i];
      const RPoint &p2 = points_[j];

      double xy = (p1.x()*p2.y()) - (p1.y()*p2.x());

      area += xy;
    }

    return area;
  }

  void moveBy(double dx, double dy) {
    RPoint d(dx, dy);

    int n = points_.size();

    for (int i = 0; i < n; ++i)
      points_[i] += d;
  }

  void expandBy(double dx, double dy) {
    // get center
    RPoint c;

    getCentroid(c);

    int n = points_.size();

    for (int i = 0; i < n; ++i) {
      //move depending on whether left or right of center
      if      (points_[i].x() < c.x()) points_[i].decX(dx);
      else if (points_[i].x() > c.x()) points_[i].incX(dx);

      //move depending on whether below or above center
      if      (points_[i].y() < c.y()) points_[i].decY(dy);
      else if (points_[i].y() > c.y()) points_[i].incY(dy);
    }
  }

  // is point (p) inside this boundary
  bool inside(const RPoint &p) const {
    return inside(points_, p);
  }

  bool intersect(const IRect &rect, Polygon &ipoly) const {
    // TODO: quick bbox check
    return intersect(Polygon(rect), ipoly);
  }

  bool intersect(const Polygon &poly, Polygon &ipoly) const {
    Points ipoints;

    if (! intersect(points_, poly.points_, ipoints))
      return false;

    ipoly.setPoints(ipoints);

    return true;
  }

 private:
  // intersect two polygons (points1) and (points2) and return the
  // result in (ipoints) if successful
  static bool intersect(const Points &points1, const Points &points2, Points &ipoints) {
    static RPoint *f[2];
    static uint    num_f;

    ipoints.clear();

    uint num_points1 = points1.size();
    uint num_points2 = points2.size();

    // fail if polygons are degenerate
    if (num_points1 < 3 || num_points2 < 3)
      return false;

    int orient1 = PolygonOrientation(points1);
    int orient2 = PolygonOrientation(points2);

    if (! orient1 || ! orient2) return false;

    // max number of intersection
    uint ni = num_points1*num_points2;

    // make sure intersection buffer is large enough
    if (num_f < ni) {
      num_f = ni;

      delete [] f[0];
      delete [] f[1];

      f[0] = new RPoint [num_f];
      f[1] = new RPoint [num_f];
    }

    // store polygon one in start point array
    // Note: if orients don't match we invert the first polygon's point order
    int l1 = 0;

    ni = num_points1;

    if (orient1 == orient2) {
      for (uint i = 0; i < ni; ++i)
        f[l1][i] = points1[i];
    }
    else {
      for (uint i = 0, j = ni - 1; i < ni; ++i, --j)
        f[l1][i] = points1[j];
    }

    // intersect current set of points with each line (end1, end2)
    // of the second polygon (points2)
    RPoint end1 = points2[num_points2 - 1];

    for (uint i = 0; i < num_points2; ++i) {
      RPoint end2 = points2[i];

      // l2 is destination point index (inverse of current l1)
      int l2 = 1 - l1;

      // calc line coefficients
      double ca = end2.x() - end1.x(); // (x2 - x1), (y2 - y1)
      double cb = end1.y() - end2.y(); // (x2 - x1), (y2 - y1)
      double cc = -end1.x()*cb - end1.y()*ca; // -x1*(y2 - y1) - y1*(x2 - x1)

      // calc side of line for first point
      RPoint v1     = f[l1][ni - 1];
      double   fv1    = ca*v1.y() + cb*v1.x() + cc;
      double   absfv1 = fabs(fv1);

      int index1 = 0;

      if (absfv1 >= 1E-6)
        index1 = sign(fv1)*orient2;

      int ni1 = 0;

      for (uint j = 0; j < ni; ++j) {
        // calc side of line for second point
        RPoint v2     = f[l1][j];
        double   fv2    = ca*v2.y() + cb*v2.x() + cc;
        double   absfv2 = fabs(fv2);

        int index2 = 0;

        if (absfv2 >= 1E-6)
          index2 = sign(fv2)*orient2;

        // add start point
        if (index1 >= 0)
          f[l2][ni1++] = v1;

        // add intersection point (if changed sides)
        if (index1 != 0 && index1 != index2 && index2 != 0) {
          double delta = absfv1 + absfv2;

          double xi = (absfv2*v1.x() + absfv1*v2.x())/delta;
          double yi = (absfv2*v1.y() + absfv1*v2.y())/delta;

          f[l2][ni1++] = RPoint(xi, yi);
        }

        // move to next line
        v1     = v2;
        absfv1 = absfv2;
        index1 = index2;
      }

      // degenerate result so fail
      if (ni1 < 3)
        return false;

      l1   = l2;
      end1 = end2;
      ni   = ni1;
    }

    ipoints.resize(ni);

    for (uint i = 0; i < ni; ++i)
      ipoints[i] = f[l1][i];

    return true;
  }

  // is points inside the boundary specified by the array of points
  static bool inside(const Points &points, const RPoint &point) {
    uint num_points = points.size();

    int counter = 0;

    int           i2     = num_points - 1;
    const RPoint *point2 = &points[i2];

    // iterate through all lines of the polygon
    for (int i1 = 0; i1 < (int) num_points; ++i1) {
      const RPoint *point1 = &points[i1];

      // intersect current line with horizontal line at inside point
      if (point.y() > std::min(point1->y(), point2->y())) {
        if (point.y() <= std::max(point1->y(), point2->y())) {
          if (point.x() <= std::max(point1->x(), point2->x())) {
            if (point1->y() != point2->y()) {
              // if we have an intersection, increase count
              double xinters = (point . y() - point1->y())*(point2->x() - point1->x())/
                               (point2->y() - point1->y()) + point1->x();

              if (point1->x() == point2->x() || point.x() <= xinters)
                ++counter;
            }
          }
        }
      }

      // next line
      i2     = i1;
      point2 = point1;
    }

    // if odd then success
    return ((counter % 2) != 0);
  }

  //! Orientation on polygon - clockwise/anti-clockwise
  static int PolygonOrientation(const Points &points) {
    assert(points.size() >= 3);

    const RPoint &point1 = points[0];
    const RPoint &point2 = points[1];
    const RPoint &point3 = points[2];

    RPoint d1(point2.x() - point1.x(), point2.y() - point1.y());
    RPoint d2(point3.x() - point2.x(), point3.y() - point2.y());

    return sign(d1.x()*d2.y() - d1.y()*d2.x());
  }

 private:
  Points        points_;
  mutable RRect bbox_;
};

//---------

class Boundary;

typedef boost::shared_ptr<Boundary> BoundaryP;

// boundary base class
class Boundary {
 public:
  enum Type {
    RECT,
    RECT_LIST,
    POLY,
    POLY_LIST
  };

  typedef std::vector<IRect>   IRects;
  typedef std::vector<Polygon> Polygons;

 public:
  Boundary() { }

  virtual Type type() const = 0;

  bool isRect() const { return type() == RECT; }

  const IRect &getRect() const {
    if (! rect_.isSet())
      const_cast<Boundary *>(this)->calcRect();

    return rect_;
  }

  void setRect(const IRect &rect) { rect_ = rect; changed(); }

  void move(int dx, int dy) { moveImpl(dx, dy); changed(); }

  virtual void moveImpl(int dx, int dy) = 0;

  void expandBy(int dx, int dy) { expandByImpl(dx, dy); changed(); }

  virtual void expandByImpl(int dx, int dy) = 0;

  bool overlaps(BoundaryP boundary) const {
    if (type() < boundary->type())
      return boundary->overlapsImpl(this);
    else
      return overlapsImpl(boundary.get());
  }

  virtual bool overlapsImpl(const Boundary *boundary) const = 0;

  bool intersect(BoundaryP boundary, BoundaryP &iboundary) const {
    if (type() < boundary->type())
      return boundary->intersectImpl(this, iboundary);
    else
      return intersectImpl(boundary.get(), iboundary);
  }

  bool intersect(const IRect &rect, BoundaryP &iboundary) const;

  virtual bool intersectImpl(const Boundary *boundary, BoundaryP &iboundary) const = 0;

  virtual bool contains(const IPoint &p) const = 0;

  virtual void getCentroid(IPoint &c) const = 0;

  virtual double getArea() const = 0;

  virtual IPoint getInsidePoint() const = 0;

  virtual void getRects(std::vector<IRect> &rects) const = 0;

  const std::vector<IPoint> &getEnclosure() const {
    if (enclosure_.empty())
      const_cast<Boundary *>(this)->calcEnclosure();

    return enclosure_;
  }

  virtual void calcEnclosure() = 0;

  virtual void print(std::ostream &os) const = 0;

  friend std::ostream &operator<<(std::ostream &os, const BoundaryP &boundary) {
    boundary->print(os);

    return os;
  }

 protected:
  virtual void changed() { enclosure_.clear(); }

  virtual void calcRect() = 0;

 protected:
  IRect               rect_;
  std::vector<IPoint> enclosure_;
};

//------

// rectangle boundary class
class RectBoundary : public Boundary {
 public:
  RectBoundary(const IRect &rect) :
   Boundary() {
    rect_ = rect;
  }

  Type type() const { return RECT; }

  void moveImpl(int dx, int dy);

  void expandByImpl(int dx, int dy);

  bool overlapsImpl(const Boundary *boundary) const;

  bool intersectImpl(const Boundary *boundary, BoundaryP &iboundary) const;

  bool contains(const IPoint &p) const;

  void getCentroid(IPoint &c) const;

  double getArea() const;

  IPoint getInsidePoint() const { return rect_.ll(); }

  void getRects(std::vector<IRect> &rects) const;

  void calcEnclosure();

  void print(std::ostream &os) const;

  friend class RectListBoundary;
  friend class PolyBoundary;
  friend class PolyListBoundary;

 private:
  void calcRect() { }
};

typedef boost::shared_ptr<RectBoundary> RectBoundaryP;

//------

// rectangle list boundary class
class RectListBoundary : public Boundary {
 public:
  RectListBoundary(const IRects &rects=IRects()) :
   Boundary(), rects_(rects) {
  }

  Type type() const { return RECT_LIST; }

  uint size() const { return rects_.size(); }

  const IRect &rect(uint i) const { return rects_[i]; }

  void addRect(const IRect &r) {
    rects_.push_back(r);

    changed();
  }

  void moveImpl(int dx, int dy);

  void expandByImpl(int dx, int dy);

  bool overlapsImpl(const Boundary *boundary) const;

  bool intersectImpl(const Boundary *boundary, BoundaryP &iboundary) const;

  bool contains(const IPoint &p) const;

  void getCentroid(IPoint &c) const;

  double getArea() const;

  IPoint getInsidePoint() const { return rects_[0].ll(); }

  void getRects(std::vector<IRect> &rects) const;

  void calcEnclosure();

  void print(std::ostream &os) const;

  friend class RectBoundary;
  friend class PolyBoundary;
  friend class PolyListBoundary;

 private:
  void changed() { rect_.reset(); Boundary::changed(); }

  void calcRect() {
    uint n = rects_.size();

    for (uint i = 0; i < n; ++i)
      rect_.combine(rects_[i]);
  }

 private:
  IRects rects_;
};

typedef boost::shared_ptr<RectListBoundary> RectListBoundaryP;

//------

// polygon boundary class
class PolyBoundary : public Boundary {
 public:
  PolyBoundary(const Polygon &poly=Polygon()) :
   Boundary(), poly_(poly) {
  }

  Type type() const { return POLY; }

  const Polygon &getPolygon() const { return poly_; }

  void addPoint(const IPoint &p) {
    poly_.addPoint(RPoint(p.x(), p.y()));

    changed();
  }

  void moveImpl(int dx, int dy);

  void expandByImpl(int dx, int dy);

  bool overlapsImpl(const Boundary *boundary) const;

  bool intersectImpl(const Boundary *boundary, BoundaryP &iboundary) const;

  bool contains(const IPoint &p) const;

  void getCentroid(IPoint &c) const;

  double getArea() const;

  IPoint getInsidePoint() const { return poly_.point(0).ipoint(); }

  void getRects(std::vector<IRect> &rects) const;

  void calcEnclosure();

  void print(std::ostream &os) const;

  friend class RectBoundary;
  friend class RectListBoundary;
  friend class PolyListBoundary;

 private:
  void changed() { rect_.reset(); Boundary::changed(); }

  void calcRect() {
    const RRect &rect = poly_.getRect();

    rect_ = IRect(floor(rect.left ()), floor(rect.bottom()),
                  ceil (rect.right()), ceil (rect.top   ()));
  }

 private:
  Polygon poly_;
};

typedef boost::shared_ptr<PolyBoundary> PolyBoundaryP;

//------

// polygon list boundary class
class PolyListBoundary : public Boundary {
 public:
  PolyListBoundary(const Polygons &polys=Polygons()) :
   Boundary(), polys_(polys) {
  }

  Type type() const { return POLY_LIST; }

  uint size() const { return polys_.size(); }

  const Polygon &polygon(uint i) const { return polys_[i]; }

  void addPolygon(const Polygon &poly) {
    polys_.push_back(poly);

    changed();
  }

  void moveImpl(int dx, int dy);

  void expandByImpl(int dx, int dy);

  bool overlapsImpl(const Boundary *boundary) const;

  bool intersectImpl(const Boundary *boundary, BoundaryP &iboundary) const;

  bool contains(const IPoint &p) const;

  void getCentroid(IPoint &c) const;

  double getArea() const;

  IPoint getInsidePoint() const { return polys_[0].point(0).ipoint(); }

  void getRects(std::vector<IRect> &rects) const;

  void calcEnclosure();

  void print(std::ostream &os) const;

  friend class RectBoundary;
  friend class RectListBoundary;
  friend class PolyBoundary;

 private:
  void changed() { rect_.reset(); Boundary::changed(); }

  void calcRect() {
    RRect rect;

    uint n = polys_.size();

    for (uint i = 0; i < n; ++i)
      rect.combine(polys_[i].getRect());

    rect_ = IRect(floor(rect.left ()), floor(rect.bottom()),
                  ceil (rect.right()), ceil (rect.top   ()));
  }

 private:
  Polygons polys_;
};

typedef boost::shared_ptr<PolyListBoundary> PolyListBoundaryP;

//------

// boundary list class
class BoundaryList {
 public:
  typedef std::vector<BoundaryP> BoundaryArray;

 public:
  BoundaryList() :
   boundaries_() {
  }

  void clear() { boundaries_.clear(); }

  void addBoundary(BoundaryP boundary) {
    boundaries_.push_back(boundary);
  }

  bool empty() const { return boundaries_.empty(); }

  uint size() const { return boundaries_.size(); }

  BoundaryP getBoundary(int i) const {
    assert(i >= 0 && i < int(size()));

    return boundaries_[i];
  }

  bool contains(const IPoint &p) const {
    uint n = size();

    for (uint i = 0; i < n; ++i)
      if (boundaries_[i]->contains(p))
        return true;

    return false;
  }

  BoundaryArray::const_iterator begin() const { return boundaries_.begin(); }
  BoundaryArray::const_iterator end  () const { return boundaries_.end  (); }

 private:
  BoundaryArray boundaries_;
};

#endif
