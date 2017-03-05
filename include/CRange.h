#ifndef CRANGE_H
#define CRANGE_H

class CRange {
 public:
  CRange() :
   set_(false), low_(0), high_(0) {
  }

  CRange(const double &low, const double &high) :
   set_(true), low_(low), high_(high) {
    fixup();
  }

  bool isSet() const { return set_; }

  void reset() { set_ = false; }

  void set(const double &low, const double &high) {
    set_ = true; low_ = low; high_ = high;

    fixup();
  }

  bool get(double *low, double *high) const {
    *low = low_; *high = high_;

    return set_;
  }

  double low () const { assert(set_); return low_ ; }
  double high() const { assert(set_); return high_; }

  double low (const double &d) const { return (set_ ? low_ : d); }
  double high(const double &d) const { return (set_ ? high_: d); }

  void setLow (const double &t) { set_ = true; low_  = t; fixup(); }
  void setHigh(const double &t) { set_ = true; high_ = t; fixup(); }

  double length() const { return (high_ - low_); }

  void move(const double &d) { assert(set_); low_ += d; high_ += d; }

  void expand(const double &d) { assert(set_); low_ -= d; high_ += d; }

  void expand(const double &d1, const double &d2) { assert(set_); low_ -= d1; high_ += d2; }

  CRange &operator=(const CRange &range) {
    set_ = range.set_; low_ = range.low_; high_ = range.high_;

    return *this;
  }

 private:
  void fixup() {
    if (low_ > high_)
      std::swap(low_, high_);
  }

 private:
  bool   set_ { false };
  double low_ { 0 }, high_ { 0 };
};

#endif
