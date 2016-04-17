#ifndef CDSIN_TABLE_H
#define CDSIN_TABLE_H

#include <cmath>
#include <CMathMacros.h>

class CDSinTable {
 private:
  uint    size_;
  double *lookup_;
  double  di_, idi_;

 public:
  CDSinTable(uint size=500) :
   size_(size), lookup_(NULL) {
    di_  = 360.0/size_;
    idi_ = 1.0/di_;

    init();
  }

 ~CDSinTable() {
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

  double sin(double a) {
    while (a >= 360.0) a -= 360.0;
    while (a <    0.0) a += 360.0;

    a *= idi_;

    int    ia = int(a);
    double ra = a - ia;

    double s1 = lookup_[ia    ];
    double s2 = lookup_[ia + 1];

    return (s1 + ra*(s2 - s1));
  }

 private:
  void init() {
    if (lookup_ == NULL)
      lookup_ = new double [size_ + 1];

    for (uint i = 0; i <= size_; ++i)
      lookup_[i] = ::sin(DEG_TO_RAD(di_*i));
  }
};

#endif
