#ifndef CCATMULL_2D_H
#define CCATMULL_2D_H

class CCatmull2D {
 public:
  CCatmull2D(double x1, double y1, double x2, double y2,
             double x3, double y3, double x4, double y4) :
   x1_(x1), y1_(y1), x2_(x2), y2_(y2), x3_(x3), y3_(y3), x4_(x4), y4_(y4) {
  }

  void calc(double t, double *x, double *y) const {
    double tt  = t*t;
    double ttt = tt*t;

    double t1 = -0.5*ttt +  1.0*tt - 0.5*t      ;
    double t2 =  1.5*ttt -  2.5*tt         + 1.0;
    double t3 = -1.5*ttt +  2.0*tt + 0.5*t      ;
    double t4 =  0.5*ttt -  0.5*tt              ;

    *x = x1_*t1 + x2_*t2 + x3_*t3 + x4_*t4;
    *y = y1_*t1 + y2_*t2 + y3_*t3 + y4_*t4;
  }

 private:
  double x1_ { 0.0 }, y1_ { 0.0 };
  double x2_ { 0.0 }, y2_ { 0.0 };
  double x3_ { 0.0 }, y3_ { 0.0 };
  double x4_ { 0.0 }, y4_ { 0.0 };
};

#endif
