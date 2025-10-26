#ifndef CPARTICLE_3D_H
#define CPARTICLE_3D_H

#include <CVector3D.h>
#include <CRGBA.h>

class CParticle3D;

class CParticleSystem3D {
 public:
  using ParticleList = std::vector<CParticle3D *>;

 public:
  CParticleSystem3D() { }

  virtual ~CParticleSystem3D() { }

  void setGravity(double gravity) { gravity_ = CVector3D(0, gravity, 0); }
  double getGravity() const { return gravity_.getY(); }

  const CVector3D &getGravityVector() const { return gravity_; }

  const ParticleList &getParticles() const { return particles_; }

  size_t numParticles() const { return particles_.size(); }

  CParticle3D *particle(size_t i) const { return particles_[i]; }

  CParticle3D *addParticle();

  virtual CParticle3D *createParticle();

  void step(double dt);
  void age();

  int maxParticles() const { return maxParticles_; }
  void setMaxParticles(int n) { maxParticles_ = n; }

  int maxAge() const { return maxAge_; }
  void setMaxAge(int n) { maxAge_ = n; }

 private:
  CVector3D    gravity_      { 0, 9.81, 0 };
  ParticleList particles_;
  int          maxParticles_ { -1 };
  int          maxAge_       { -1 };
};

//---

class CParticle3D {
 public:
  enum class State {
    ALIVE,
    DEAD
  };

 public:
  CParticle3D(const CParticleSystem3D &system) :
   system_(system) {
  }

  virtual ~CParticle3D() { }

  void init() {
    mass_ = 0.0;
    age_  = 0;

    color_ = CRGBA(1, 1, 1);
    state_ = State::ALIVE;

    position_    .zero();
    velocity_    .zero();
    acceleration_.zero();
  }

  //---

  double getMass() const { return mass_; }
  void setMass(double mass) { mass_ = mass; }

  //---

  const CVector3D &getPosition() const { return position_    ; }
  void setPosition(const CVector3D &position) { position_ = position; }
  void setPosition(double x, double y, double z) { setPosition(CVector3D(x, y, z)); }

  void incPosition(const CVector3D &position) { position_ += position; }

  //---

  const CVector3D &getVelocity() const { return velocity_    ; }
  void setVelocity(const CVector3D &velocity) { velocity_ = velocity; }
  void setVelocity(double x, double y, double z) { setVelocity(CVector3D(x, y, z)); }

  void incVelocity(const CVector3D &velocity) { velocity_ += velocity; }

  //---

  const CVector3D &getAcceleration() const { return acceleration_; }
  void setAcceleration(const CVector3D &acceleration) { acceleration_ = acceleration; }
  void setAcceleration(double x, double y, double z) { setAcceleration(CVector3D(x, y, z)); }

  void incAcceleration(const CVector3D &acceleration) { acceleration_ += acceleration; }

  //---

  void setForce(double x, double y, double z) { acceleration_ = CVector3D(x, y, z)/getMass(); }
  void setForce(const CVector3D &f) { acceleration_ = f/getMass(); }

  void incForce(const CVector3D &f) { acceleration_ += f/getMass(); }

  //---

  uint getAge() const { return age_; }
  void setAge(const uint &i) { age_ = i; }

  void incAge() { ++age_; }

  //---

  const CRGBA &getColor() const { return color_; }
  void setColor(const CRGBA &v) { color_ = v; }

  //---

  const State &getState() const { return state_; }
  void setState(const State &v) { state_ = v; }

  void setAlive() { state_ = State::ALIVE; }
  void setDead () { state_ = State::DEAD ; }

  bool isAlive() const { return state_ == State::ALIVE; }
  bool isDead () const { return state_ == State::DEAD ; }

  //---

  void setFixed(bool fixed=false) { fixed_ = fixed; }
  bool isFixed() const { return fixed_; }

  //---

  void step(double t) {
    auto a = acceleration_ - system_.getGravityVector();

    auto v = a*t;

    position_ += (velocity_ + 0.5*v)*t;

    velocity_ += v;
  }

 private:
  const CParticleSystem3D &system_;

  double    mass_         { 0.0 };
  CVector3D position_     { 0, 0, 0 };
  CVector3D velocity_     { 0, 0, 0 };
  CVector3D acceleration_ { 0, 0, 0 };
  uint      age_          { 0 };
  CRGBA     color_        { 1, 1, 1 };
  State     state_        { State::DEAD };
  bool      fixed_        { false };
};

#endif
