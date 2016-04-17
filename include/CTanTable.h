#ifndef CTAN_TABLE_H
#define CTAN_TABLE_H

#include <cmath>

class CTanTable {
 private:
  uint    size_;
  double *lookup_;
  double  di_, idi_;

 public:
  CTanTable(uint size=500) :
   size_(size), lookup_(NULL) {
    di_  = 2.0*M_PI/size_;
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
    di_     = 2.0*M_PI/size_;
    idi_    = 1.0/di_;
    lookup_ = NULL;

    init();
  }

  double tan(double a) {
    while (a >= 2.0*M_PI) a -= 2.0*M_PI;
    while (a <       0.0) a += 2.0*M_PI;

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
      lookup_[i] = ::tan(di_*i);
  }
};

#endif
