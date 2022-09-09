#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include <iostream>
#include <random>
#include <cmath>
#include <math.h>
#include <TMath.h>

using namespace HepMC3;

/** Generate single particle event with fixed three momentum **/
void gen_fixed_single_particle(int n_events = 100, const int pdgid = 13, const double p = 10, const double phi = 0, const double eta = 0, const char* out_fname = "single_particle.hepmc")
{
  WriterAscii hepmc_output(out_fname);
  int events_parsed = 0;
  GenEvent evt(Units::GEV, Units::MM);
  std::cout << "generating single particle events with p = " << p << std::endl;

  for (events_parsed = 0; events_parsed < n_events; events_parsed++) {
    // FourVector(px,py,pz,e,pdgid,status)
    // type 4 is beam
    // pdgid 2212 - proton
    GenParticlePtr p1 =
        std::make_shared<GenParticle>(FourVector(0.0, 0.0, 10.0, 10.0), 11, 4);
    GenParticlePtr p2 = std::make_shared<GenParticle>(
        FourVector(0.0, 0.0, 0.0, 0.938), 2212, 4);

    // Define momentum
    Double_t th    = 2*std::atan(exp(-eta));
    Double_t px    = p * std::cos(phi) * std::sin(th);
    Double_t py    = p * std::sin(phi) * std::sin(th);
    Double_t pz    = p * std::cos(th);

    std::cout << "Generated particle momentum is (" << px << "," << py << "," << pz << ")" << std::endl;

    // type 1 is final state
    Double_t mass = 0;
    if (abs(pdgid)==13) mass = 0.10565837; // charged muon
    else if (abs(pdgid)==211) mass = 0.13957039; // charged pion
    else if (abs(pdgid)==11) mass = 0.000511; // electron/positron
    else if (abs(pdgid)==321) mass = 0.493677; // charged kaon
    else if (abs(pdgid)==2212) mass = 0.983272; // proton/anti-proton
    else
    {
      std::cout << "particle with PDG ID " << pdgid << "not implemented. Abort... " << std::endl;
      exit(0);
    }
    GenParticlePtr p3 = std::make_shared<GenParticle>(
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
    if (events_parsed % 10000 == 0) {
      std::cout << "Event: " << events_parsed << std::endl;
    }
    evt.clear();
  }
  hepmc_output.close();
  std::cout << "Events parsed and written: " << events_parsed << std::endl;
}
