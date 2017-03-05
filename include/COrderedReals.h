#ifndef CORDERED_REALS_H
#define CORDERED_REALS_H

#include <set>
#include <cmath>

class COrderedReals {
 private:
  struct Compare {
    bool operator()(const double &lhs, const double &rhs) const {
      if (fabs(lhs - rhs) < 1E-6) return false;

      return lhs < rhs;
    }
  };

 public:
  typedef typename std::set<double,Compare>::iterator       iterator;
  typedef typename std::set<double,Compare>::const_iterator const_iterator;

  COrderedReals() { }

  bool insert(double t) {
    std::pair<iterator,bool> p = reals_.insert(t);

    return p.second;
  }

  std::size_t size() {
    return reals_.size();
  }

  double min() { return *reals_.begin (); }
  double max() { return *reals_.rbegin(); }

  iterator begin() { return reals_.begin(); }
  iterator end  () { return reals_.end  (); }

  const_iterator begin() const { return reals_.begin(); }
  const_iterator end  () const { return reals_.end  (); }

 private:
  std::set<double,Compare> reals_;
};

#endif
