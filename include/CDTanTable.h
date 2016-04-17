#ifndef CDTAN_TABLE_H
#define CDTAN_TABLE_H

#include <cmath>
#include <CMathMacros.h>

class CTanTable {
 private:
  uint    size_;
  double *lookup_;
  double  di_, idi_;

 public:
  CTanTable(uint size=500) :
   size_(size), lookup_(NULL) {
    di_  = 360.0/size_;
    idi_ = 1.0/di_;

    init();
  }

 ~CTanTable() {
    delete [] lookup_;
  }

  uint getSize() const { return size_; }

  void setSize(uint size) {
    if (size == size_) return;

    delete [] lookup_;

    size_   = size;
    di_     = 360.0/size_;
    idi_    = 1.0/di_;
    lookup_ = NULL;

    init();
  }

  double tan(double a) {
    while (a >= 360.0) a -= 360.0;
    while (a <    0.0) a += 360.0;

    a *= idi_;

    int    ia = int(a);
    double ra = a - ia;

    double t1 = lookup_[ia    ];
    double t2 = lookup_[ia + 1];

    return (t1 + ra*(t2 - t1));
  }

 private:
  void init() {
    if (lookup_ == NULL)
      lookup_ = new double [size_ + 1];

    for (uint i = 0; i <= size_; ++i)
      lookup_[i] = ::tan(DEG_TO_RAD(di_*i));
  }
};

#endif
