#include <CParticle3D.h>
#include <map>

CParticle3D *
CParticleSystem3D::
addParticle()
{
  CParticle3D *particle = nullptr;

  int numAlive = 0;

  using Particles     = std::vector<CParticle3D *>;
  using AgedParticles = std::map<int, Particles>;

  AgedParticles agedParticles;

  for (auto *particle1 : particles_) {
    if (particle1->isDead()) {
      // reuse first dead particle
      particle1->init();
      return particle1;
    }

    agedParticles[-particle1->getAge()].push_back(particle1);

    ++numAlive;
  }

  // reuse oldest particle
  if (maxParticles() > 0 && numAlive >= maxParticles()) {
    auto *particle2 = agedParticles.begin()->second.back();
    particle2->init();
    return particle2;
  }

  // create new particle
  particle = createParticle();

  particles_.push_back(particle);

  particle->init();

  return particle;
}

CParticle3D *
CParticleSystem3D::
createParticle()
{
  return new CParticle3D(*this);
}

void
CParticleSystem3D::
step(double dt)
{
  for (auto *particle : getParticles())
    particle->step(dt);
}

void
CParticleSystem3D::
age()
{
  for (auto *particle : getParticles()) {
    particle->incAge();
  }
}
