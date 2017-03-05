#ifndef CIRANGE_H
#define CIRANGE_H

class CIRange {
 public:
  CIRange() :
   set_(false), low_(0), high_(0) {
  }

  CIRange(const int &low, const int &high) :
   set_(true), low_(low), high_(high) {
    fixup();
  }

  bool isSet() const { return set_; }

  void reset() { set_ = false; }

  void set(const int &low, const int &high) {
    set_ = true; low_ = low; high_ = high;

    fixup();
  }

  bool get(int *low, int *high) const {
    *low = low_; *high = high_;

    return set_;
  }

  int low () const { assert(set_); return low_ ; }
  int high() const { assert(set_); return high_; }

  int low (const int &d) const { return (set_ ? low_ : d); }
  int high(const int &d) const { return (set_ ? high_: d); }

  void setLow (const int &t) { set_ = true; low_  = t; fixup(); }
  void setHigh(const int &t) { set_ = true; high_ = t; fixup(); }

  int length() const { return (high_ - low_); }

  void move(const int &d) { assert(set_); low_ += d; high_ += d; }

  void expand(const int &d) { assert(set_); low_ -= d; high_ += d; }

  void expand(const int &d1, const int &d2) { assert(set_); low_ -= d1; high_ += d2; }

  CIRange &operator=(const CIRange &range) {
    set_ = range.set_; low_ = range.low_; high_ = range.high_;

    return *this;
  }

 private:
  void fixup() {
    if (low_ > high_)
      std::swap(low_, high_);
  }

 private:
  bool set_ { false };
  int  low_ { 0 }, high_ { 0 };
};

#endif
