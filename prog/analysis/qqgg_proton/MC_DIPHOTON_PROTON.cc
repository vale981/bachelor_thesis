// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {

/// @brief Generate some simple histograms of the diphoton process. TODO: more
class MC_DIPHOTON_PROTON : public Analysis {
public:
  /// Constructor
  DEFAULT_RIVET_ANALYSIS_CTOR(MC_DIPHOTON_PROTON);

  /// @name Analysis methods
  //@{

  /// Book histograms and initialise projections before the run
  void init() {
    // Initialise and register projections

    // The basic final-state projection:
    // all final-state particles within
    // the given eta acceptance
    // We hardcode this TODO: make this configurable
    FinalState fs;
    declare(fs, "FS");

    // cut has been made in sherpa
    IdentifiedFinalState ifs {};
    ifs.acceptId(PID::PHOTON);
    declare(ifs, "IFS");

    auto energy = info().energies()[0].first;
    book(_h_pT, "pT", 50, 1000, energy);
    book(_h_eta, "eta", 50, -1, 1);
    book(_h_cos_theta, "cos_theta", 50, -.986, .986);
  }

  /// Perform the per-event analysis
  void analyze(const Event &event) {
    const double weight = 1.0;

    const Particles &photons = apply<IdentifiedFinalState>(event, "IFS").particles();

    // they are both the same, so we take the first
    const auto &photon = photons.front();

    _h_pT->fill(photon.pT(), weight);
    _h_eta->fill(photon.eta(), weight);
    _h_cos_theta->fill(cos(photon.theta()), weight);
  }

  //@}

  void finalize() {
    normalize(_h_pT);
    normalize(_h_eta);
    normalize(_h_cos_theta);
  }

  /// @name Histograms
  //@{
  Histo1DPtr _h_pT;
  Histo1DPtr _h_eta;
  Histo1DPtr _h_cos_theta;
  //@}
};

DECLARE_RIVET_PLUGIN(MC_DIPHOTON_PROTON);

} // namespace Rivet
