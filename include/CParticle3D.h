#ifndef CPARTICLE_3D_H
#define CPARTICLE_3D_H

#include <CVector3D.h>

class CParticle3D {
 public:
  CParticle3D(double m=0.0) :
   mass_(m), fixed_(false) {
  }

 ~CParticle3D() { }

  void setPosition(double x, double y, double z) {
    position_ = CVector3D(x, y, z);
  }

  void setPosition(const CVector3D &position) {
    position_ = position;
  }

  void incPosition(const CVector3D &position) {
    position_ += position;
  }

  void setVelocity(double x, double y, double z) {
    velocity_ = CVector3D(x, y, z);
  }

  void setVelocity(const CVector3D &velocity) {
    velocity_ = velocity;
  }

  void incVelocity(const CVector3D &velocity) {
    velocity_ += velocity;
  }

  void setAcceleration(double x, double y, double z) {
    acceleration_ = CVector3D(x, y, z);
  }

  void setAcceleration(const CVector3D &acceleration) {
    acceleration_ = acceleration;
  }

  void incAcceleration(const CVector3D &acceleration) {
    acceleration_ += acceleration;
  }

  void setForce(double x, double y, double z) {
    acceleration_ = CVector3D(x, y, z)/mass_;
  }

  void setForce(const CVector3D &f) {
    acceleration_ = f/mass_;
  }

  void incForce(const CVector3D &f) {
    acceleration_ += f/mass_;
  }

  void setMass(double mass) {
    mass_ = mass;
  }

  void setFixed(bool fixed=false) {
    fixed_ = fixed;
  }

  const CVector3D &getPosition    () const { return position_    ; }
  const CVector3D &getVelocity    () const { return velocity_    ; }
  const CVector3D &getAcceleration() const { return acceleration_; }
  double           getMass        () const { return mass_        ; }
  bool             isFixed        () const { return fixed_       ; }

 private:
  CVector3D position_;
  CVector3D velocity_;
  CVector3D acceleration_;
  double    mass_;
  bool      fixed_;
};

#endif
