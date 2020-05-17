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

  /// Create a logarithmicly spaced vector
  template <typename T = double>
  std::vector<T> logspace(T min, T max, size_t count, bool endpoint = false) {
    std::vector<T> result(count);
    min = std::log10(min);
    max = std::log10(max);

    T step = (max - min) / (count - (endpoint ? 1 : 0));
    std::generate(result.begin(), result.end(), [&, n = 0] () mutable {
      return std::pow(10, min + (n ++) * step);
    });

    return result;
  }

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
    double min_pT = 20;
    double eta = 2.5;

    book(_h_pT, "pT", logspace(min_pT, energy, 51, true));
    book(_h_eta, "eta", 50, -eta, eta);
    book(_h_cos_theta, "cos_theta", 50, -1, 1);
    book(_h_inv_m, "inv_m", logspace(2 * min_pT, 2 * energy, 51, true));
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
    const auto &moms = photons.moms();
    const auto total_momentum = std::accumulate(moms.begin(), moms.end(), FourMomentum(0, 0, 0, 0));

    _h_inv_m->fill(total_momentum.mass());
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
  Histo1DPtr _h_inv_m;
  //@}
};

DECLARE_RIVET_PLUGIN(MC_DIPHOTON_PROTON);

} // namespace Rivet
