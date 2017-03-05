#include <CRigidBody.h>

//////////////////////
// CRigidBodyState //
//////////////////////

void
CRigidBodyState::
calculateInertiaTensor(const CVector3D &p)
{
  const CVector3D ip = CVector3D(1/p.getX(), 1/p.getY(), 1/p.getZ());

  double a, b, c, d, e, f, g, h, i;

  getMatrix().getValues(&a, &b, &c, &d, &e, &f, &g, &h, &i);

  //since it's a symmetric matrix, there are only 6 moments and products
  i_.setValues(p.getX()*a*a + p.getY()*b*b + p.getZ()*c*c,
               p.getX()*a*d + p.getY()*b*e + p.getZ()*c*f,
               p.getX()*a*g + p.getY()*b*h + p.getZ()*c*i,
               p.getX()*d*d + p.getY()*e*e + p.getZ()*f*f,
               p.getX()*d*g + p.getY()*e*h + p.getZ()*f*i,
               p.getX()*g*g + p.getY()*h*h + p.getZ()*i*i);

  i_inv_.setValues(ip.getX()*a*a + ip.getY()*b*b + ip.getZ()*c*c,
                   ip.getX()*a*d + ip.getY()*b*e + ip.getZ()*c*f,
                   ip.getX()*a*g + ip.getY()*b*h + ip.getZ()*c*i,
                   ip.getX()*d*d + ip.getY()*e*e + ip.getZ()*f*f,
                   ip.getX()*d*g + ip.getY()*e*h + ip.getZ()*f*i,
                   ip.getX()*g*g + ip.getY()*h*h + ip.getZ()*i*i);
}
