#ifndef CDCOS_TABLE_H
#define CDCOS_TABLE_H

#include <cmath>
#include <CMathMacros.h>

class CDCosTable {
 public:
  CDCosTable(uint size=500) :
   size_(size) {
    di_  = 360.0/size_;
    idi_ = 1.0/di_;

    init();
  }

 ~CDCosTable() {
    delete [] lookup_;
  }

  uint getSize() const { return size_; }

  void setSize(uint size) {
    if (size == size_) return;

    delete [] lookup_;

    size_   = size;
    di_     = 360.0/size_;
    idi_    = 1.0/di_;
    lookup_ = nullptr;

    init();
  }

  double cos(double a) {
    while (a >= 360.0) a -= 360.0;
    while (a <    0.0) a += 360.0;

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
      lookup_[i] = ::cos(DEG_TO_RAD(di_*i));
  }

 private:
  uint    size_   { 500 };
  double *lookup_ { nullptr };
  double  di_     { 1 };
  double  idi_    { 1 };
};

#endif
