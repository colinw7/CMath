#ifndef CMATH_FUNCTOR_H
#define CMATH_FUNCTOR_H

class Functor {
 public:
  virtual ~Functor() { }

  double operator()(double x) const { return exec(x); }

  virtual double exec(double x) const = 0;

  virtual double gradient(double x) const = 0;
};

class CosFunctor : public Functor {
 public:
  double exec(double x) const override { return cos(x); }

  double gradient(double x) const override { return -sin(x); }
};

class SinFunctor : public Functor {
 public:
  double exec(double x) const override { return sin(x); }

  double gradient(double x) const override { return cos(x); }
};

class TanFunctor : public Functor {
 public:
  double exec(double x) const override { return tan(x); }

  double gradient(double x) const override { double t = tan(x); return 1 + t*t; }
};

class PolarFunctor {
 public:
  virtual ~PolarFunctor() { }

  double operator()(double theta) const { return exec(theta); }

  virtual double exec(double theta) const = 0;

  virtual double gradient(double theta) const = 0;
};

class ArchimedianSpiralFunctor : public PolarFunctor {
 public:
  ArchimedianSpiralFunctor(double a=0.0, double n=1.0) :
   a_(a), n_(n) {
    assert(a_ != 0.0 && n_ != 0.0);
  }

  virtual ~ArchimedianSpiralFunctor() { }

  double exec(double theta) const override {
    return a_*pow(theta, 1.0/n_);
  }

  double gradientAngle(double theta) const {
    return theta + atan(theta);
  }

  double gradient(double theta) const override {
    return tan(gradientAngle(theta));
  }

  double curvature(double theta) const {
    double theta2 = theta*theta;

    return (2 + theta2)/(a_*pow(1 + theta2, 1.5));
  }

 private:
  double a_, n_;
};

class ArchimedesSpiralFunctor : public PolarFunctor {
 public:
  ArchimedesSpiralFunctor(double a=0.0) :
   a_(a) {
    assert(a_ != 0.0);
  }

  virtual ~ArchimedesSpiralFunctor() { }

  double exec(double theta) const override {
    return a_*theta;
  }

  double gradientAngle(double theta) const {
    return theta + atan(theta);
  }

  double gradient(double theta) const override {
    return tan(gradientAngle(theta));
  }

  double curvature(double theta) const {
    double theta2 = theta*theta;

    return (2 + theta2)/(a_*pow(1 + theta2, 1.5));
  }

  ArchimedianSpiralFunctor inverseFunctor() {
    return ArchimedianSpiralFunctor(a_, -1);
  }

 private:
  double a_;
};

#endif
