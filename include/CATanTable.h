#ifndef CATAN_TABLE_H
#define CATAN_TABLE_H

#include <cmath>
#include <CMathMacros.h>

class CATanTable {
 public:
  CATanTable(uint size=500) :
   size_(size) {
    init();
  }

 ~CATanTable() {
    delete [] lookup_;
  }

  uint getSize() const { return size_; }

  void setSize(uint size) {
    if (size == size_) return;

    delete [] lookup_;

    size_   = size;
    lookup_ = nullptr;

    init();
  }

  double atan(double y, double x) {
    if (x >= 0) {
      if (y < 0) return -atan_1(-y, x);
      else       return  atan_1( y, x);
    }
    else {
      if (y < 0) return -M_PI + atan_1(-y, -x);
      else       return  M_PI - atan_1( y, -x);
    }
  }

 private:
  double atan_1(double y, double x) {
    if (x < y) return M_PI/2.0 - atan_2(x, y);
    else       return atan_2(y, x);
  }

  double atan_2(double y, double x) {
    if (x == 0) return 0;

    double a  = size_*(y/x);
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

    double di = 1.0/size_;

    for (uint i = 0; i <= size_; ++i)
      lookup_[i] = ::atan(i*di);

    lookup_[size_ + 1] = lookup_[size_];
  }

 private:
  uint    size_   { 500 };
  double *lookup_ { nullptr };
};

#endif
