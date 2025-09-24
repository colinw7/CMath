#ifndef CPARTICLE_2D_H
#define CPARTICLE_2D_H

#include <CVector2D.h>
#include <CRGBA.h>
#include <accessor.h>
#include <vector>

class CParticle2D;

class CParticleSystem2D {
 public:
  using ParticleList = std::vector<CParticle2D *>;

 public:
  CParticleSystem2D() { }

  virtual ~CParticleSystem2D() { }

  void setGravity(double gravity) { gravity_ = CVector2D(0, gravity); }
  double getGravity() const { return gravity_.getY(); }

  const CVector2D &getGravityVector() const { return gravity_; }

  const ParticleList &getParticles() const { return particles_; }

  size_t numParticles() const { return particles_.size(); }

  CParticle2D *particle(size_t i) const { return particles_[i]; }

  CParticle2D *addParticle();

  virtual CParticle2D *createParticle();

  void step(double dt);

 private:
  CVector2D    gravity_;
  ParticleList particles_;
};

//---

class CParticle2D {
 public:
  enum class State {
    ALIVE,
    DEAD
  };

 public:
  CParticle2D(const CParticleSystem2D &system) :
   system_(system) {
  }

  virtual ~CParticle2D() { }

  void init() {
    mass_ = 0.0;
    age_  = 0;

    color_ = CRGBA(1, 1, 1);
    state_ = State::ALIVE;

    position_    .zero();
    velocity_    .zero();
    acceleration_.zero();
  }

  ACCESSOR(Mass        , double   , mass)
  ACCESSOR(Position    , CVector2D, position)
  ACCESSOR(Velocity    , CVector2D, velocity)
  ACCESSOR(Acceleration, CVector2D, acceleration)
  ACCESSOR(Age         , uint     , age)
  ACCESSOR(Color       , CRGBA    , color)
  ACCESSOR(State       , State    , state)

  void setPosition(double x, double y) { position_ = CVector2D(x, y); }

  void setVelocity(double x, double y) { velocity_ = CVector2D(x, y); }

  void setAcceleration(double x, double y) { acceleration_ = CVector2D(x, y); }

  void incVelocity(const CVector2D &velocity) { velocity_ += velocity; }

  void incAcceleration(const CVector2D &acceleration) { acceleration_ += acceleration; }

  void incAge() { ++age_; }

  void setAlive() { state_ = State::ALIVE; }
  void setDead () { state_ = State::DEAD ; }

  bool isAlive() const { return state_ == State::ALIVE; }
  bool isDead () const { return state_ == State::DEAD ; }

  void step(double t) {
    auto a = acceleration_ - system_.getGravityVector();

    auto v = a*t;

    position_ += (velocity_ + 0.5*v)*t;

    velocity_ += v;
  }

 private:
  const CParticleSystem2D &system_;

  double    mass_         { 0.0 };
  CVector2D position_     { 0, 0 };
  CVector2D velocity_     { 0, 0 };
  CVector2D acceleration_ { 0, 0 };
  uint      age_          { 0 };
  CRGBA     color_        { 1, 1, 1 };
  State     state_        { State::DEAD };
  bool      fixed_        { false };
};

#endif
