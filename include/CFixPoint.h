#ifndef CFIXED_POINT_H
#define CFIXED_POINT_H

// Unsigned fixed point float
// Need operator< ...

template<uint W, uint D>
class CUnsignedFixPoint {
 private:
  enum {
    WHOLE_MAX   = 1U<<W,
    DECIMAL_MAX = 1U<<D,
  };

  uint w_:W;
  uint d_:D;

 public:
  CUnsignedFixPoint(uint w, uint d = 0) :
   w_(w), d_(d) {
    if (w >= WHOLE_MAX || d >= DECIMAL_MAX)
      assert(false && "Overflow");
  }

  CUnsignedFixPoint(double r) :
   w_(0), d_(0) {
    if (r < 0)
      assert(false && "Negative");

    uint w = (uint)(r);
    uint d = (uint)((r - w)*DECIMAL_MAX + 0.5);

    if (w >= WHOLE_MAX || d >= DECIMAL_MAX)
      assert(false && "Overflow");

    w_ = w;
    d_ = d;
  }

  CUnsignedFixPoint operator+(const CUnsignedFixPoint &rhs) {
    uint w1 = w_ + rhs.w_;
    uint d1 = d_ + rhs.d_;

    uint w2 = d1 >> D;

    d1 -= w2 << D;

    return CUnsignedFixPoint(w1 + w2, d1);
  }

  CUnsignedFixPoint operator-(const CUnsignedFixPoint &rhs) {
    if (w_ < rhs.w_ || (w_ == rhs.w_ && d_ < rhs.d_)) {
      assert(false && "Negative");
      return CUnsignedFixPoint(0,0);
    }

    uint w1 = w_ - rhs.w_;
    int  d1 = d_ - rhs.d_;

    if (d1 < 0) {
      d1 += DECIMAL_MAX + 1;

      if (w1 <= 0) {
        assert(false && "Negative");
        return CUnsignedFixPoint(0,0);
      }

      --w1;
    }

    return CUnsignedFixPoint(w1, d1);
  }

  CUnsignedFixPoint operator*(const CUnsignedFixPoint &rhs) {
    uint w1 = w_*rhs.w_;

    uint d1 = w_*rhs.d_ + d_* rhs.w_ + ((d_ * rhs.d_) >> D);

    uint w2 = d1 >> D;

    d1 -= w2 << D;

    return CUnsignedFixPoint(w1 + w2, d1);
  }

  CUnsignedFixPoint operator/(const CUnsignedFixPoint &rhs) {
    uint t1 = (    w_ << D) +     d_;
    uint t2 = (rhs.w_ << D) + rhs.d_;

    uint w = t1 / t2;
    uint r = t1 - w*t2;

    while (r > WHOLE_MAX) {
      r  >>= 1;
      t2 >>= 1;
    }

    uint d = (r << D)/t2;

    return CUnsignedFixPoint(w, d);
  }

  void print(ostream &os) const {
    uint d1 = (uint)((1000000.0*d_)/DECIMAL_MAX + 0.5);

    os << w_ << ".";

    char buffer[6];

    sprintf(buffer, "%06u", d1);

    os << buffer;
  }

  friend ostream &operator<<(ostream &os, const CUnsignedFixPoint &fp) {
    fp.print(os);

    return os;
  }
};

template<uint W, uint D>
class CFixPoint {
 private:
  enum {
    WHOLE_MAX   = 1<<W,
    DECIMAL_MAX = 1<<D,
  };

  uint s_:1;
  uint w_:W;
  uint d_:D;

 public:
  CFixPoint(int w, uint d = 0) :
   s_(w < 0 ? 1 : 0), w_(abs(w)), d_(d) {
    if (abs(w) >= WHOLE_MAX || d >= DECIMAL_MAX)
      assert(false && "Overflow");
  }

  CFixPoint(double r) :
   s_(r < 0 ? 1 : 0), w_(0), d_(0) {
    uint w = (uint)(fabs(r));
    uint d = (uint)((fabs(r) - w)*DECIMAL_MAX + 0.5);

    if (w >= WHOLE_MAX || d >= DECIMAL_MAX)
      assert(false && "Overflow");

    w_ = w;
    d_ = d;
  }

  CFixPoint abs() const {
    return CFixPoint(0, w_, d_);
  }

  CFixPoint operator-() const {
    return CFixPoint(s_ ? 0 : 1, w_, d_);
  }

  CFixPoint operator+(const CFixPoint &rhs) const {
    if      (! s_ && ! rhs.s_) {
      uint w1 = w_ + rhs.w_;
      uint d1 = d_ + rhs.d_;

      uint w2 = d1 >> D;

      d1 -= w2 << D;

      return CFixPoint(0, w1 + w2, d1);
    }
    else if (! s_ &&   rhs.s_)
      return   *this - rhs.abs();
    else if (  s_ && ! rhs.s_)
      return -(abs() - rhs      );
    else
      return -(abs() + rhs.abs());
  }

  CFixPoint operator-(const CFixPoint &rhs) const {
    if      (! s_ && ! rhs.s_) {
      if (w_ < rhs.w_ || (w_ == rhs.w_ && d_ < rhs.d_))
        return -(rhs - *this);
      else {
        uint w1 = w_ - rhs.w_;
        int  d1 = d_ - rhs.d_;

        if (d1 < 0) {
          d1 += DECIMAL_MAX + 1;

          --w1;
        }

        return CFixPoint(0, w1, d1);
      }
    }
    else if (! s_ &&   rhs.s_)
      return   *this + rhs.abs();
    else if (  s_ && ! rhs.s_)
      return -(abs() + rhs      );
    else
      return -(abs() - rhs.abs());
  }

  CFixPoint operator*(const CFixPoint &rhs) const {
    uint s = (s_ != rhs.s_ ? 1 : 0);

    uint w1 = w_*rhs.w_;

    uint d1 = w_*rhs.d_ + d_* rhs.w_ + ((d_ * rhs.d_) >> D);

    uint w2 = d1 >> D;

    d1 -= w2 << D;

    return CFixPoint(s, w1 + w2, d1);
  }

  CFixPoint operator/(const CFixPoint &rhs) const {
    uint s = (s_ != rhs.s_ ? 1 : 0);

    uint t1 = (    w_ << D) +     d_;
    uint t2 = (rhs.w_ << D) + rhs.d_;

    uint w = t1 / t2;
    uint r = t1 - w*t2;

    while (r > WHOLE_MAX) {
      r  >>= 1;
      t2 >>= 1;
    }

    uint d = (r << D)/t2;

    return CFixPoint(s, w, d);
  }

  void print(ostream &os) const {
    if (s_)
      os << "-";

    uint d1 = (uint)((1000000.0*d_)/DECIMAL_MAX + 0.5);

    os << w_ << ".";

    char buffer[6];

    sprintf(buffer, "%06u", d1);

    os << buffer;
  }

  friend ostream &operator<<(ostream &os, const CFixPoint &fp) {
    fp.print(os);

    return os;
  }

 private:
  CFixPoint(uint s, uint w, uint d) :
   s_(s), w_(w), d_(d) {
  };
};

#endif
