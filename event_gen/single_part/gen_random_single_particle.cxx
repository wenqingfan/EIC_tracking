#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include <iostream>
#include <stdlib.h> 
#include <random>
#include <cmath>
#include <math.h>
#include <TMath.h>

using namespace HepMC3;

const int verbosity = 0;

int get_seed()
{
    int fSeed;

    ifstream devrandom;
    devrandom.open("/dev/urandom",ios::binary);
    devrandom.read((char*)&fSeed,sizeof(fSeed));
    devrandom.close();
    if ( fSeed != -1 )
    {
      cout << "Got seed from /dev/urandom" << endl;
      fSeed = abs(fSeed)%900000000;
    }
    else fSeed = 0;
    cout << "Seed is " << fSeed << endl;
    return fSeed;
}

/** Generate single particle event with random three momentum **/
void gen_random_single_particle(int n_events = 1000, const int pdgid = 13, float p_min = 0.2, const float p_max = 2, float eta_min = -3.5, const float eta_max = 3.5, const char* out_fname = "single_particle.hepmc")
{
  WriterAscii hepmc_output(out_fname);
  int events_parsed = 0;
  GenEvent evt(Units::GEV, Units::MM);
  if (abs(pdgid)==211)  std::cout << "Generating single pion events" << std::endl;
  else if (abs(pdgid)==321)  std::cout << "Generating single kaon events" << std::endl;
  else if (abs(pdgid)==13)   std::cout << "Generating single muon events" << std::endl;
  else if (abs(pdgid)==11)   std::cout << "Generating single electon events" << std::endl;
  else if (abs(pdgid)==2212) std::cout << "Generating single proton events" << std::endl;
  else
  {
    cout << "Particle with PDG ID " << pdgid << "not implemented. Abort... " << endl;
    exit(0);
  }
  cout << "-- p range from " << p_min << " to " << p_max << " GeV" << endl;
  cout << "-- eta range from " << eta_min << " to " << eta_max << endl;

  int seed = get_seed();
  gRandom->SetSeed(seed);

  TF1 flat("flat","1.0",0.0,1.0);
      
  for (events_parsed = 0; events_parsed < n_events; events_parsed++) {
    // FourVector(px,py,pz,e,pdgid,status)
    // type 4 is beam
    // pdgid 2212 - proton
    GenParticlePtr p1 =
        std::make_shared<GenParticle>(FourVector(0.0, 0.0, 10.0, 10.0), 11, 4);
    GenParticlePtr p2 = std::make_shared<GenParticle>(
        FourVector(0.0, 0.0, 0.0, 0.938), 2212, 4);

    Double_t p  = flat.GetRandom()*(p_max-p_min)+p_min; // flat in [pt_min, pt_max]GeV
    Double_t phi = 2*TMath::Pi()*flat.GetRandom();
    Double_t eta = flat.GetRandom()*(eta_max-eta_min)+eta_min;

    // Define momentum
    Double_t th    = 2 * atan(exp(-eta));
    Double_t px    = p * cos(phi) * sin(th);
    Double_t py    = p * sin(phi) * sin(th);
    Double_t pz    = p * cos(th);

    if (events_parsed%100==0) cout << "Generated particle with p " << p << " and eta " << eta << " => momentum is (" << px << "," << py << "," << pz << ")" << endl;

    // type 1 is final state
    Double_t mass = 0;
    if (abs(pdgid)==13) mass = 0.10565837; // charged muon
    if (abs(pdgid)==211) mass = 0.13957039; // charged pion
    if (abs(pdgid)==11) mass = 0.000511; // electron/positron
    if (abs(pdgid)==321) mass = 0.493677; // charged kaon
    if (abs(pdgid)==2212) mass = 0.983272; // proton/anti-proton

    GenParticlePtr p3 = make_shared<GenParticle>(
        FourVector(px, py, pz, sqrt(p*p + (mass * mass))), pdgid, 1); // 1 means stable particle

    GenVertexPtr v1 = std::make_shared<GenVertex>();
    v1->add_particle_in(p1);
    v1->add_particle_in(p2);

    v1->add_particle_out(p3);
    evt.add_vertex(v1);

    if (events_parsed == 0) {
      std::cout << "First event: " << std::endl;
      Print::listing(evt);
    }

    hepmc_output.write_event(evt);
    if (events_parsed % 1000 == 0) {
      std::cout << "Event: " << events_parsed << std::endl;
    }
    evt.clear();
  }
  hepmc_output.close();
  std::cout << "Events parsed and written: " << events_parsed << std::endl;
}
