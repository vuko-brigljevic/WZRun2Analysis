#include "Particles.h"
#include "WZEvent.h"
#include "Constants.h"


using namespace std;


// Initialize static members
WZEvent* Particle::fWZTree = 0;


void Particle::SetWZEvent(WZEvent* wzt)
{
  fWZTree = wzt;
}


Particle::Particle(unsigned int index, double pt, double eta, double phi)
{
  SetPtEtaPhiM(pt, eta, phi, 0.);
  fIndex = index;
}


GenParticle::GenParticle(unsigned int index, double pt, double eta, double phi)
  : Particle(index, pt, eta, phi)
{
  if (fWZTree == 0)  cout << "WZEvent pointer is ZERO!!!! \n";

  fPdgId = fWZTree->mcPID->at(fIndex);
}

