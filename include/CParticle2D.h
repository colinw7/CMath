#ifndef CPARTICLE_2D_H
#define CPARTICLE_2D_H

#include <CVector2D.h>
#include <CRGBA.h>
#include <accessor.h>
#include <vector>

class CParticle2D;

class CParticleSystem2D {
 public:
  typedef std::vector<CParticle2D *> ParticleList;

 public:
  CParticleSystem2D() { }

  virtual ~CParticleSystem2D() { }

  void setGravity(double gravity) { gravity_ = CVector2D(0, gravity); }

  double getGravity() const { return gravity_.getY(); }

  const CVector2D &getGravityVector() const { return gravity_; }

  const ParticleList &getParticles() const { return particles_; }

  CParticle2D *addParticle();

  virtual CParticle2D *createParticle();

 private:
  CVector2D    gravity_;
  ParticleList particles_;
};

//---

class CParticle2D {
 public:
  enum State {
    ALIVE,
    DEAD
  };

 public:
  CParticle2D(const CParticleSystem2D &system) :
   system_(system), mass_(0.0), age_(0), state_(DEAD) {
  }

  virtual ~CParticle2D() { }

  void init() {
    mass_ = 0;
    age_  = 0;

    color_ = CRGBA(1,1,1);
    state_ = ALIVE;

    position_    .zero();
    velocity_    .zero();
    acceleration_.zero();
  }

  ACCESSOR(Mass        ,      double, mass)
  ACCESSOR(Position    , CVector2D, position)
  ACCESSOR(Velocity    , CVector2D, velocity)
  ACCESSOR(Acceleration, CVector2D, acceleration)
  ACCESSOR(Age         , uint  , age)
  ACCESSOR(Color       , CRGBA , color)
  ACCESSOR(State       , State , state)

  void setPosition(double x, double y) {
    position_ = CVector2D(x, y);
  }

  void setVelocity(double x, double y) {
    velocity_ = CVector2D(x, y);
  }

  void setAcceleration(double x, double y) {
    acceleration_ = CVector2D(x, y);
  }

  void incVelocity(const CVector2D &velocity) {
    velocity_ += velocity;
  }

  void incAcceleration(const CVector2D &acceleration) {
    acceleration_ += acceleration;
  }

  void incAge() {
    ++age_;
  }

  void setAlive() { state_ = ALIVE; }
  void setDead () { state_ = DEAD ; }

  bool isAlive() const { return state_ == ALIVE; }
  bool isDead () const { return state_ == DEAD ; }

  void step(double t) {
    CVector2D a = acceleration_ - system_.getGravityVector();

    CVector2D v = a*t;

    position_ += (velocity_ + 0.5*v)*t;

    velocity_ += v;
  }

 private:
  const CParticleSystem2D &system_;

  double    mass_;
  CVector2D position_;
  CVector2D velocity_;
  CVector2D acceleration_;
  uint      age_;
  CRGBA     color_;
  State     state_;
};

#endif
