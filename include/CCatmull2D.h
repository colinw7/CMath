#ifndef CCATMULL_2D_H
#define CCATMULL_2D_H

template<typename T>
class CCatmull2DT {
 private:
  T x1_, y1_;
  T x2_, y2_;
  T x3_, y3_;
  T x4_, y4_;

 public:
  CCatmull2DT(T x1, T y1, T x2, T y2, T x3, T y3, T x4, T y4) :
   x1_(x1), y1_(y1), x2_(x2), y2_(y2), x3_(x3), y3_(y3), x4_(x4), y4_(y4) {
  }

  void calc(T t, T *x, T *y) const {
    T tt  = t*t;
    T ttt = tt*t;

    T t1 = -0.5*ttt +  1.0*tt - 0.5*t      ;
    T t2 =  1.5*ttt -  2.5*tt         + 1.0;
    T t3 = -1.5*ttt +  2.0*tt + 0.5*t      ;
    T t4 =  0.5*ttt -  0.5*tt              ;

    *x = x1_*t1 + x2_*t2 + x3_*t3 + x4_*t4;
    *y = y1_*t1 + y2_*t2 + y3_*t3 + y4_*t4;
  }
};

typedef CCatmull2DT<double> CCatmull2D;

#endif
