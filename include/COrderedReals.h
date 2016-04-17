#ifndef CORDERED_REALS_H
#define CORDERED_REALS_H

template<class T>
class COrderedRealsT {
 private:
  struct Compare {
    bool operator()(const T &lhs, const T &rhs) const {
      if (fabs(lhs - rhs) < 1E-6) return false;

      return lhs < rhs;
    }
  };

  set<T,Compare> reals_;

 public:
  typedef typename set<T,Compare>::iterator       iterator;
  typedef typename set<T,Compare>::const_iterator const_iterator;

  COrderedRealsT() { }

  bool insert(T t) {
    std::pair<iterator,bool> p = reals_.insert(t);

    return p.second;
  }

  uint size() {
    return reals_.size();
  }

  T min() { return *reals_.begin (); }
  T max() { return *reals_.rbegin(); }

  iterator begin() { return reals_.begin(); }
  iterator end  () { return reals_.end  (); }

  const_iterator begin() const { return reals_.begin(); }
  const_iterator end  () const { return reals_.end  (); }
};

typedef COrderedRealsT<double> COrderedReals;

#endif
