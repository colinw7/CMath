#include <CParticle2D.h>

CParticle2D *
CParticleSystem2D::
addParticle()
{
  // reuse dead particle
  for (auto *particle : particles_) {
    if (particle->isDead()) {
      particle->init();
      return particle;
    }
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
