#include <CParticle2D.h>

CParticle2D *
CParticleSystem2D::
addParticle()
{
  typename ParticleList::iterator p1 = particles_.begin(), p2 = particles_.end();

  for ( ; p1 != p2; ++p1)
    if ((*p1)->isDead()) {
      (*p1)->init();

    return *p1;
  }

  CParticle2D *particle = createParticle();

  particles_.push_back(particle);

  particle->init();

  return particle;
}

CParticle2D *
CParticleSystem2D::
createParticle()
{
  return new CParticle2D(*this);
}
