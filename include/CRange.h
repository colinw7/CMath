#ifndef CRANGE_H
#define CRANGE_H

template<typename T>
class CRangeT {
 public:
  CRangeT() :
   set_(false), low_(0), high_(0) {
  }

  CRangeT(const T &low, const T &high) :
   set_(true), low_(low), high_(high) {
    fixup();
  }

  bool isSet() const { return set_; }

  void reset() { set_ = false; }

  void set(const T &low, const T &high) {
    set_ = true; low_ = low; high_ = high;

    fixup();
  }

  bool get(T *low, T *high) const {
    *low = low_; *high = high_;

    return set_;
  }

  T low () const { assert(set_); return low_ ; }
  T high() const { assert(set_); return high_; }

  T low (const T &d) const { return (set_ ? low_ : d); }
  T high(const T &d) const { return (set_ ? high_: d); }

  void setLow (const T &t) { set_ = true; low_  = t; fixup(); }
  void setHigh(const T &t) { set_ = true; high_ = t; fixup(); }

  T length() const { return (high_ - low_); }

  void move(const T &d) { assert(set_); low_ += d; high_ += d; }

  void expand(const T &d) { assert(set_); low_ -= d; high_ += d; }

  void expand(const T &d1, const T &d2) { assert(set_); low_ -= d1; high_ += d2; }

  CRangeT &operator=(const CRangeT &range) {
    set_ = range.set_; low_ = range.low_; high_ = range.high_;

    return *this;
  }

 private:
  void fixup() {
    if (low_ > high_)
      std::swap(low_, high_);
  }

 private:
  bool set_;
  T    low_, high_;
};

typedef CRangeT<double> CRange;
typedef CRangeT<int>    CIRange;

#endif
