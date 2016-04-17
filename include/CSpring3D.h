#ifndef CSpring3D_H
#define CSpring3D_H

#include <CParticle3D.h>

class CSpring3D {
 private:
  typedef CParticle3D Particle;
  typedef CVector3D   Vector;

 public:
  CSpring3D(const Particle &a, const Particle &b, double ks, double d, double r) :
   a_(a), b_(b) {
    springConstant_ = ks;
    damping_        = d;
    restLength_     = r;
    on_             = true;
  }

  void turnOff() {
    on_ = false;
  }

  void turnOn() {
    on_ = true;
  }

  bool isOn() const {
    return on_;
  }

  bool isOff() const {
    return ! on_;
  }

  const Particle &getOneEnd() const {
    return a_;
  }

  const Particle &getTheOtherEnd() const {
    return b_;
  }

  double currentLength() const {
    return a_.getPosition().getDistance(b_.getPosition());
  }

  double restLength() const {
    return restLength_;
  }

  double strength() const {
    return springConstant_;
  }

  void setStrength(double ks) {
    springConstant_ = ks;
  }

  double damping() const {
    return damping_;
  }

  void setDamping(double d) {
    damping_ = d;
  }

  void setRestLength(double l) {
    restLength_ = l;
  }

  void apply() {
    if (on_ && (! a_.isFixed() || ! b_.isFixed())) {
      double a2bX = a_.getPosition().x() - b_.getPosition().x();
      double a2bY = a_.getPosition().y() - b_.getPosition().y();
      double a2bZ = a_.getPosition().z() - b_.getPosition().z();

      double a2bDistance = sqrt(a2bX*a2bX + a2bY*a2bY + a2bZ*a2bZ);

      if (a2bDistance != 0) {
        a2bX /= a2bDistance;
        a2bY /= a2bDistance;
        a2bZ /= a2bDistance;
      }
      else {
        a2bX = 0;
        a2bY = 0;
        a2bZ = 0;
      }

      // spring force is proportional to how much it stretched
      double springForce = -(a2bDistance - restLength_)*springConstant_;

      // want velocity along line b/w a & b, damping force is proportional to this
      double Va2bX = a_.getVelocity().x() - b_.getVelocity().x();
      double Va2bY = a_.getVelocity().y() - b_.getVelocity().y();
      double Va2bZ = a_.getVelocity().z() - b_.getVelocity().z();

      double dampingForce = -damping_*(a2bX*Va2bX + a2bY*Va2bY + a2bZ*Va2bZ);

      // forceB is same as forceA in opposite direction
      double r = springForce + dampingForce;

      a2bX *= r;
      a2bY *= r;
      a2bZ *= r;

      if (! a_.isFixed()) a_.incForce(Vector( a2bX,  a2bY,  a2bZ));
      if (! b_.isFixed()) b_.incForce(Vector(-a2bX, -a2bY, -a2bZ));
    }
  }

 protected:
  void setA(const Particle &p) { a_ = p; }
  void setB(const Particle &p) { b_ = p; }

 private:
  double   springConstant_;
  double   damping_;
  double   restLength_;
  Particle a_;
  Particle b_;
  bool     on_;
};

#endif
