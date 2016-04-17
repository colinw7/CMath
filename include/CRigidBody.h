#ifndef CRIGID_BODY_H
#define CRIGID_BODY_H

#include <CCoordFrame3D.h>
#include <CSymmetricMatrix3D.h>

class CRigidBodyState : public CCoordFrame3D {
 protected:
  CSymmetricMatrix3D i_;     // inertia tensor with respect to
                             // parent frame, kg*m^2
  CSymmetricMatrix3D i_inv_; // inverse inertia tensor
  CVector3D          v_;     // linear velocity, meters/sec
  CVector3D          w_;     // angular velocity, radians/sec

 public:
  CRigidBodyState() { }

  CRigidBodyState(const CVector3D &v, const CVector3D &w) : v_(v), w_(w) { }

  void reset() {
    CCoordFrame3D::reset();
  }

  const CVector3D &linearVelocity () { return v_; }
  const CVector3D &angularVelocity() { return w_; }

  void linearVelocity (const CVector3D &v) { v_ = v; }
  void angularVelocity(const CVector3D &w) { w_ = w; }

  // calculate the inertia tensor and its
  // inverse from the current orientation
  // and the principal moments of inertia
  void calculateInertiaTensor(const CVector3D &ip);
};

//
// Describe a rigid body
//
class CRigidBody : public CRigidBodyState {
 protected:
  // previous dynamic state
  CRigidBodyState prev_;

  // mass and inertia
  double          mass_;  // mass, kg
  CVector3D       ip_;    // principal moments of inertia, kg * m^2

  // bounding box dimensions
  CVector3D       dim_;

 public:
  CRigidBody() : mass_(0) { }

  CRigidBody(const double m, const CVector3D &d,
             const CVector3D ip = CVector3D(0,0,0)) :
   mass_(m), ip_(ip), dim_(d) {
    reset();
  }

  void reset() {
    CRigidBodyState::reset();

    // if no moments were passed in, give the moments for a box
    if (ip_.getX() == 0 && ip_.getY() == 0 && ip_.getZ() == 0) {
      double x2 = dim_.getX()*dim_.getX();
      double y2 = dim_.getY()*dim_.getY();
      double z2 = dim_.getZ()*dim_.getZ();

      ip_ = (mass_/12) * CVector3D(y2 + z2, z2 + x2, x2 + y2);
    }

    // before any physics are done
    calculateInertiaTensor(ip_);
  }

  // physics
  void mass(const double m) { mass_ = m; }

  double mass() const { return mass_; }

  void principalInertia(const CVector3D &ip) { ip_ = ip; }

  const CVector3D &principalInertia() const { return ip_; }

  void dim(const CVector3D &d) { dim_ = d; }

  const CVector3D &dim() const { return dim_; }

  const CRigidBodyState currentState() const { return *this; }
  const CRigidBodyState previousState() const { return prev_; }

  // time derivative of angular velocity
  const CVector3D dwdt(const CVector3D &T) const {
    CVector3D v1;

    i_.multiplyVector(w_, v1);

    CVector3D v2;

    i_inv_.multiplyVector(T - w_.crossProduct(v1), v2);

    return v2;
  }

  void rotate(const CVector3D &v) {
    CRigidBodyState::rotate(v);

    calculateInertiaTensor(ip_);
  }
};

#endif
