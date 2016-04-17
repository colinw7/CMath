#ifndef CRAND_H
#define CRAND_H

class CRand {
 private:
  unsigned long long seed_;
  unsigned long long mult_;
  unsigned long long llong_max_;
  float              float_max_;

 public:
  CRand(unsigned long long seed = 7564231ULL) :
   seed_(seed), mult_(6208991ULL), llong_max_(4294967295ULL),
   float_max_(4294967295.0) {
  }

  float operator()() {
    seed_ *= mult_;

    return float(seed_ % llong_max_) / float_max_;
  }
};

class CRandRange {
  private:
   double low_;
   double high_;
   double range_;
   CRand  r;

  public:
   CRandRange(double low=0.0, double high=1.0) :
    low_(low), high_(high) {
     range_ = high_ - low;
   }

   double operator()(double low, double high) {
     low_   = low;
     range_ = high - low;

     return r()*range_ + low_;
   }

   double operator()() {
     return r()*range_ + low_;
   }
};

#endif
