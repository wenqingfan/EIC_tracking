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

/** Generate single muon event with fixed three momentum **/
void gen_fixed_single_muon(int n_events = 100, const double p = 10, const double phi = 0, const double eta = 0, const char* out_fname = "single_muon.hepmc")
{
  WriterAscii hepmc_output(out_fname);
  int events_parsed = 0;
  GenEvent evt(Units::GEV, Units::MM);
  std::cout << "generating single muon events with p = " << p << std::endl;

  for (events_parsed = 0; events_parsed < n_events; events_parsed++) {
    // FourVector(px,py,pz,e,pdgid,status)
    // type 4 is beam
    // pdgid 13 - muon (neg charge)
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

    std::cout << "Generated muon momentum is (" << px << "," << py << "," << pz << ")" << std::endl;

    // type 1 is final state
    // pdgid 13 - muon 105.65837 MeV/c^2
    GenParticlePtr p3 = std::make_shared<GenParticle>(
        FourVector(
            px, py, pz,
            sqrt(p*p + (0.10565837 * 0.10565837))),
        13, 1);

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
