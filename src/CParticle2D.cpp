#include <CParticle2D.h>

CParticle2D *
CParticleSystem2D::
addParticle()
{
  for (auto p1 = particles_.begin(); p1 != particles_.end(); ++p1)
    if ((*p1)->isDead()) {
      (*p1)->init();

    return *p1;
  }

  auto *particle = createParticle();

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

void
CParticleSystem2D::
step(double dt)
{
  for (auto *particle : getParticles())
    particle->step(dt);
}
