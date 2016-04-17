#ifndef CRUNGE_KUTTA_H
#define CRUNGE_KUTTA_H

template<typename FnX>
class CRungeKuttaX {
 private:
  FnX    fnx_;
  double t_;
  double x_;
  double dt_, dt2_, dt6_;

 public:
  // f(t,x) step dt
  CRungeKuttaX(double t, double x, double dt=0.01) :
   t_(t), x_(x), dt_(dt) {
    dt2_ = dt_/2.0;
    dt6_ = dt_/6.0;
  }

  double getX() const { return x_; }

  void step() {
    double k1 = fnx_(t_       , x_         );
    double k2 = fnx_(t_ + dt2_, x_ + k1/2.0);
    double k3 = fnx_(t_ + dt2_, x_ + k2/2.0);
    double k4 = fnx_(t_ + dt_ , x_ + k3    );

    x_ += dt6_*(k1 + 2*k2 + 2*k3 + k4);
    t_ += dt_;
  }
};

template<typename FnX, typename FnY, typename FnZ>
class CRungeKuttaXYZ {
 private:
  FnX    fnx_;
  FnY    fny_;
  FnZ    fnz_;
  double t_;
  double x_, y_, z_;
  double dt_, dt2_, dt6_;

 public:
  // f1(t,x), f2(t,y), f3(t,z) step dt
  CRungeKuttaXYZ(double t, double x, double y, double z, double dt=0.01) :
   t_(t), x_(x), y_(y), z_(z), dt_(dt) {
    dt2_ = dt_/2.0;
    dt6_ = dt_/6.0;
  }

  double getX() const { return x_; }
  double getY() const { return y_; }
  double getZ() const { return z_; }

  void step() {
    double kx1 = fnx_(t_, x_, y_, z_);
    double ky1 = fny_(t_, x_, y_, z_);
    double kz1 = fnz_(t_, x_, y_, z_);

    double dt = t_ + dt2_;

    double xt = x_ + dt_*kx1/2.0;
    double yt = y_ + dt_*ky1/2.0;
    double zt = z_ + dt_*kz1/2.0;

    double kx2 = fnx_(dt, xt, yt, zt);
    double ky2 = fny_(dt, xt, yt, zt);
    double kz2 = fnz_(dt, xt, yt, zt);

    xt = x_ + dt_*kx2/2.0;
    yt = y_ + dt_*ky2/2.0;
    zt = z_ + dt_*kz2/2.0;

    double kx3 = fnx_(dt, xt, yt, zt);
    double ky3 = fny_(dt, xt, yt, zt);
    double kz3 = fnz_(dt, xt, yt, zt);

    xt = x_ + dt_*kx3;
    yt = y_ + dt_*ky3;
    zt = z_ + dt_*kz3;

    double kx4 = fnx_(dt, xt, yt, zt);
    double ky4 = fny_(dt, xt, yt, zt);
    double kz4 = fnz_(dt, xt, yt, zt);

    x_ += dt6_*(kx1 + 2.0*kx2 + 2.0*kx3 + kx4);
    y_ += dt6_*(ky1 + 2.0*ky2 + 2.0*ky3 + ky4);
    z_ += dt6_*(kz1 + 2.0*kz2 + 2.0*kz3 + kz4);
    t_ += dt_;
  }
};

#endif
