#ifndef CMinMax_H
#define CMinMax_H

template<typename T>
class CMinMax {
 public:
  explicit CMinMax() = default;

  explicit CMinMax(const T &t) {
    add(t);
  }

  explicit CMinMax(const T &t1, const T &t2) {
    add(t1);
    add(t2);
  }

  //---

  void add(const T &t) {
    if (! set_) {
      min_ = t;
      max_ = t;
      set_ = true;
    }
    else {
      min_ = std::min(t, min_);
      max_ = std::max(t, max_);
    }
  }

  bool isSet() const { return set_; }

  void reset() { set_ = false; }

  const T &min() const { assert(set_); return min_; }
  const T &max() const { assert(set_); return max_; }

  T min(const T &t) const { return (set_ ? min_ : t); }
  T max(const T &t) const { return (set_ ? max_ : t); }

  //---

  T map(const T &r, const T &min, const T &max) const {
    return CMathUtil::map(r, min_, max_, min, max);
  }

  T normalize(const T &r) const {
    auto d = max() - min();

    if (d > T(0))
      return (r - min())/d;
    else
      return min();
  }

  //---
bool inside(const T &r) const {
    if (! set_) return false;

    return (r >= min_ && r <= max_);
  }

  bool insideHalfOpen(const T &r) const {
    if (! set_) return false;

    return (r >= min_ && r < max_);
  }

  bool overlaps(const CMinMax &r) const {
    if (! set_ || ! r.set_) return false;

    return (max_ >= r.min_ && min_ <= r.max_);
  }

  bool overlapsHalfOpen(const CMinMax &r) const {
    if (! set_ || ! r.set_) return false;

    return (max_ > r.min_ && min_ < r.max_);
  }

 private:
  T    min_ { };
  T    max_ { };
  bool set_ { false };
};

using CRMinMax = CMinMax<double>;
using CIMinMax = CMinMax<int>;

#endif
