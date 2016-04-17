#ifndef CSIN_TABLE_H
#define CSIN_TABLE_H

#include <cmath>

class CSinTable {
 private:
  uint    size_;
  double *lookup_;
  double  di_, idi_;

 public:
  CSinTable(uint size=500) :
   size_(size), lookup_(NULL) {
    di_  = 2.0*M_PI/size_;
    idi_ = 1.0/di_;

    init();
  }

 ~CSinTable() {
    delete [] lookup_;
  }

  uint getSize() const { return size_; }

  void setSize(uint size) {
    if (size == size_) return;

    delete [] lookup_;

    size_   = size;
    di_     = 2.0*M_PI/size_;
    idi_    = 1.0/di_;
    lookup_ = NULL;

    init();
  }

  double sin(double a) {
    while (a >= 2.0*M_PI) a -= 2.0*M_PI;
    while (a <       0.0) a += 2.0*M_PI;

    a *= idi_;

    int    ia = int(a);
    double ra = a - ia;

    double s1 = lookup_[ia    ];
    double s2 = lookup_[ia + 1];

    return (s1 + ra*(s2 - s1));
  }

 private:
  void init() {
    if (! lookup_)
      lookup_ = new double [size_ + 1];

    for (uint i = 0; i <= size_; ++i)
      lookup_[i] = ::sin(di_*i);
  }
};

#endif
