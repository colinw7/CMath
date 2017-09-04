#ifndef CCOS_TABLE_H
#define CCOS_TABLE_H

#include <cmath>

class CCosTable {
 public:
  CCosTable(uint size=500) :
   size_(size) {
    di_  = 2.0*M_PI/size_;
    idi_ = 1.0/di_;

    init();
  }

 ~CCosTable() {
    delete [] lookup_;
  }

  uint getSize() const { return size_; }

  void setSize(uint size) {
    if (size == size_) return;

    delete [] lookup_;

    size_   = size;
    di_     = 2.0*M_PI/size_;
    idi_    = 1.0/di_;
    lookup_ = nullptr;

    init();
  }

  double cos(double a) {
    while (a >= 2.0*M_PI) a -= 2.0*M_PI;
    while (a <       0.0) a += 2.0*M_PI;

    a *= idi_;

    int    ia = int(a);
    double ra = a - ia;

    double c1 = lookup_[ia    ];
    double c2 = lookup_[ia + 1];

    return (c1 + ra*(c2 - c1));
  }

 private:
  void init() {
    if (! lookup_)
      lookup_ = new double [size_ + 1];

    for (uint i = 0; i <= size_; ++i)
      lookup_[i] = ::cos(di_*i);
  }

 private:
  uint    size_   { 500 };
  double *lookup_ { nullptr };
  double  di_     { 1 };
  double  idi_    { 1 };
};

#endif
