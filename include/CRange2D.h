#ifndef CRANGE_2D_H
#define CRANGE_2D_H

template<typename T>
class CRange2DT {
 public:
  CRange2DT() :
   set_(false), x1_(0), y1_(0), x2_(0), y2_(0) {
  }

  CRange2DT(T x1, T y1, T x2, T y2) :
   set_(true), x1_(x1), y1_(y1), x2_(x2), y2_(y2) {
  }

  bool isSet() const { return set_; }

  void set(T x1, T y1, T x2, T y2) {
    set_ = true;

    x1_ = x1; y1_ = y1;
    x2_ = x2; y2_ = y2;
  }

  bool get(T *x1, T *y1, T *x2, T *y2) const {
    *x1 = x1_; *y1 = y1_;
    *x2 = x2_; *y2 = y2_;

    return set_;
  }

  T dx() const { assert(set_); return x2_ - x1_; }
  T dy() const { assert(set_); return y2_ - y1_; }

  T xmid() const { assert(set_); return (x2_ + x1_)/2; }
  T ymid() const { assert(set_); return (y2_ + y1_)/2; }

  T xmin() const { assert(set_); return std::min(x1_, x2_); }
  T ymin() const { assert(set_); return std::min(y1_, y2_); }
  T xmax() const { assert(set_); return std::max(x1_, x2_); }
  T ymax() const { assert(set_); return std::max(y1_, y2_); }

  T left  () const { assert(set_); return x1_; }
  T bottom() const { assert(set_); return y1_; }
  T right () const { assert(set_); return x2_; }
  T top   () const { assert(set_); return y2_; }

  void setLeft  (const T &t) { set_ = true; x1_ = t; }
  void setBottom(const T &t) { set_ = true; y1_ = t; }
  void setRight (const T &t) { set_ = true; x2_ = t; }
  void setTop   (const T &t) { set_ = true; y2_ = t; }

  T xsize() const { assert(set_); return fabs(x2_ - x1_); }
  T ysize() const { assert(set_); return fabs(y2_ - y1_); }

  void inc(T dx, T dy) {
    assert(set_);

    x1_ += dx; y1_ += dy;
    x2_ += dx; y2_ += dy;
  }

  void incX(T dx) { assert(set_); x1_ += dx; x2_ += dx; }
  void incY(T dy) { assert(set_); y1_ += dy; y2_ += dy; }

  CRange2DT &operator=(const CRange2DT &range) {
    set_ = range.set_;
    x1_  = range.x1_; y1_ = range.y1_;
    x2_  = range.x2_; y2_ = range.y2_;

    return *this;
  }

 private:
  bool set_;
  T    x1_, y1_, x2_, y2_;
};

typedef CRange2DT<double> CRange2D;
typedef CRange2DT<int>    CIRange2D;

#endif
