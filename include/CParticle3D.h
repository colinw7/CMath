#ifndef CPARTICLE_3D_H
#define CPARTICLE_3D_H

#include <CVector3D.h>

template<typename T>
class CParticle3DT {
 private:
  typedef CVector3DT<T> Vector;

 public:
  CParticle3DT(T m=0.0) :
   mass_(m), fixed_(false) {
  }

 ~CParticle3DT() { }

  void setPosition(T x, T y, T z) {
    position_ = Vector(x, y, z);
  }

  void setPosition(const Vector &position) {
    position_ = position;
  }

  void incPosition(const Vector &position) {
    position_ += position;
  }

  void setVelocity(T x, T y, T z) {
    velocity_ = Vector(x, y, z);
  }

  void setVelocity(const Vector &velocity) {
    velocity_ = velocity;
  }

  void incVelocity(const Vector &velocity) {
    velocity_ += velocity;
  }

  void setAcceleration(T x, T y, T z) {
    acceleration_ = Vector(x, y, z);
  }

  void setAcceleration(const Vector &acceleration) {
    acceleration_ = acceleration;
  }

  void incAcceleration(const Vector &acceleration) {
    acceleration_ += acceleration;
  }

  void setForce(T x, T y, T z) {
    acceleration_ = Vector(x, y, z)/mass_;
  }

  void setForce(const Vector &f) {
    acceleration_ = f/mass_;
  }

  void incForce(const Vector &f) {
    acceleration_ += f/mass_;
  }

  void setMass(T mass) {
    mass_ = mass;
  }

  void setFixed(bool fixed=false) {
    fixed_ = fixed;
  }

  const Vector &getPosition    () const { return position_    ; }
  const Vector &getVelocity    () const { return velocity_    ; }
  const Vector &getAcceleration() const { return acceleration_; }
  T             getMass        () const { return mass_        ; }
  bool          isFixed        () const { return fixed_       ; }

 private:
  Vector position_;
  Vector velocity_;
  Vector acceleration_;
  T      mass_;
  bool   fixed_;
};

typedef CParticle3DT<double> CParticle3D;

#endif
