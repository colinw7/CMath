#ifndef CEXP_TABLE_H
#define CEXP_TABLE_H

#include <cmath>
#include <cassert>

class CExpTable {
 public:
  CExpTable(double m1, double m2, unsigned int size=500) :
   m1_(m1), m2_(m2), size_(size) {
    init();
  }

 ~CExpTable() {
    delete [] lookup_;
  }

  unsigned int getSize() const { return size_; }

  void setSize(unsigned int size) {
    if (size == size_) return;

    delete [] lookup_;

    size_   = size;
    lookup_ = nullptr;

    init();
  }

  double exp(double x) {
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
    if (! lookup_)
      lookup_ = new double [size_ + 2];

    double di = (m2_ - m1_)/size_;

    for (unsigned int i = 0; i <= size_; ++i)
      lookup_[i] = ::exp(m1_ + i*di);

    lookup_[size_ + 1] = lookup_[size_];
  }

 private:
  double       m1_     { 1 };
  double       m2_     { 1 };
  unsigned int size_   { 500 };
  double*      lookup_ { nullptr };
};

#endif
