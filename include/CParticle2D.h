#ifndef CPARTICLE_2D_H
#define CPARTICLE_2D_H

#include <CVector2D.h>
#include <CRGBA.h>
#include <accessor.h>
#include <vector>

template<typename T>
class CParticle2DT;

template<typename T>
class CParticleSystem2DT {
 private:
  typedef CVector2DT<T> Vector;

 public:
  typedef CParticle2DT<T>         Particle;
  typedef std::vector<Particle *> ParticleList;

 private:
  Vector       gravity_;
  ParticleList particles_;

 public:
  CParticleSystem2DT() { }

  virtual ~CParticleSystem2DT() { }

  void setGravity(T gravity) { gravity_ = Vector(0, gravity); }

  T getGravity() const { return gravity_.getY(); }

  const Vector &getGravityVector() const { return gravity_; }

  const ParticleList &getParticles() const { return particles_; }

  Particle *addParticle() {
    typename ParticleList::iterator p1 = particles_.begin(), p2 = particles_.end();

    for ( ; p1 != p2; ++p1)
      if ((*p1)->isDead()) {
        (*p1)->init();

        return *p1;
      }

    Particle *particle = createParticle();

    particles_.push_back(particle);

    particle->init();

    return particle;
  }

  virtual Particle *createParticle() {
    return new Particle(*this);
  }
};

template<typename T>
class CParticle2DT {
 private:
  typedef CParticleSystem2DT<T> ParticleSystem;
  typedef CVector2DT<T>         Vector;

 public:
  enum State {
    ALIVE,
    DEAD
  };

 private:
  const ParticleSystem &system_;

  T      mass_;
  Vector position_;
  Vector velocity_;
  Vector acceleration_;
  uint   age_;
  CRGBA  color_;
  State  state_;

 public:
  CParticle2DT(const ParticleSystem &system) :
   system_(system), mass_(0.0), age_(0), state_(DEAD) {
  }

  virtual ~CParticle2DT() { }

  void init() {
    mass_ = 0;
    age_  = 0;

    color_ = CRGBA(1,1,1);
    state_ = ALIVE;

    position_    .zero();
    velocity_    .zero();
    acceleration_.zero();
  }

  ACCESSOR(Mass        ,      T, mass)
  ACCESSOR(Position    , Vector, position)
  ACCESSOR(Velocity    , Vector, velocity)
  ACCESSOR(Acceleration, Vector, acceleration)
  ACCESSOR(Age         , uint  , age)
  ACCESSOR(Color       , CRGBA , color)
  ACCESSOR(State       , State , state)

  void setPosition(T x, T y) {
    position_ = Vector(x, y);
  }

  void setVelocity(T x, T y) {
    velocity_ = Vector(x, y);
  }

  void setAcceleration(T x, T y) {
    acceleration_ = Vector(x, y);
  }

  void incVelocity(const Vector &velocity) {
    velocity_ += velocity;
  }

  void incAcceleration(const Vector &acceleration) {
    acceleration_ += acceleration;
  }

  void incAge() {
    ++age_;
  }

  void setAlive() { state_ = ALIVE; }
  void setDead () { state_ = DEAD ; }

  bool isAlive() const { return state_ == ALIVE; }
  bool isDead () const { return state_ == DEAD ; }

  void step(T t) {
    Vector a = acceleration_ - system_.getGravityVector();

    Vector v = a*t;

    position_ += (velocity_ + 0.5*v)*t;

    velocity_ += v;
  }
};

typedef CParticle2DT<double>       CParticle2D;
typedef CParticleSystem2DT<double> CParticleSystem2D;

#endif
