#ifndef CSQRT_TABLE_H
#define CSQRT_TABLE_H

#include <cmath>
#include <cassert>
#include <CMathMacros.h>

class CSqrtTable {
 public:
  CSqrtTable(double m1, double m2, uint size=500) :
   m1_(m1), m2_(m2), size_(size), lookup_(NULL) {
    assert(m1 >= 0.0);

    init();
  }

 ~CSqrtTable() {
    delete [] lookup_;
  }

  uint getSize() const { return size_; }

  void setSize(uint size) {
    if (size == size_) return;

    delete [] lookup_;

    size_   = size;
    lookup_ = NULL;

    init();
  }

  double sqrt(double x) {
    assert(x >= m1_ && x <= m2_);

    double x1 = (x - m1_)/(m2_ - m1_);

    double a  = size_*x1;
    int    ia = int(a);
    double ra = a - ia;

    double v1 = lookup_[ia    ];
    double v2 = lookup_[ia + 1];

    return (v1 + ra*(v2 - v1));
  }

 private:
  void init() {
    if (lookup_ == NULL)
      lookup_ = new double [size_ + 2];

    double di = (m2_ - m1_)/size_;

    for (uint i = 0; i <= size_; ++i)
      lookup_[i] = ::sqrt(m1_ + i*di);

    lookup_[size_ + 1] = lookup_[size_];
  }

 private:
  double  m1_, m2_;
  uint    size_;
  double *lookup_;
};

#endif
