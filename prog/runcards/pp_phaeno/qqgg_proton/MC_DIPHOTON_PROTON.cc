// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include <fstream>
#include <iostream>

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
    // Calorimeter particles for photon isolation
    VisibleFinalState visFS;
    VetoedFinalState calo_fs(visFS);
    calo_fs.addVetoPairId(PID::MUON); // we also want photon
    declare(calo_fs, "calo_fs");

    // we chain in a prompt final state just to be save
    PromptFinalState prompt{};
    IdentifiedFinalState ifs(prompt);
    ifs.acceptId(PID::PHOTON);
    declare(ifs, "Photons");

    auto energy = info().energies()[0].first;
    double min_pT = 20;
    double eta = 2.5;

    _observables = {"pT",      "eta",        "cos_theta", "inv_m",
                    "o_angle", "o_angle_cs", "total_pT"};

    book(_histos["pT"], "pT", logspace(50, min_pT, energy, true));
    book(_histos["eta"], "eta", 50, -eta, eta);
    book(_histos["cos_theta"], "cos_theta", 50, -1, 1);
    book(_histos["inv_m"], "inv_m", logspace(50, 2 * min_pT, 2 * energy, true));
    book(_histos["o_angle"], "o_angle", 50, 0, 1);
    book(_histos["o_angle_cs"], "o_angle_cs", 50, 0, 1);
    book(_histos["total_pT"], "total_pT", logspace(50, .01, 2 * energy));
  }

  /// Perform the per-event analysis
  void analyze(const Event &event) {
    Particles photons =
        apply<IdentifiedFinalState>(event, "Photons").particlesByPt();

    // make sure that there are only two photons
    if (photons.size() < 2)
      vetoEvent;

    photons.resize(2);

    // Require the two photons to be separated in dR
    if (deltaR(photons[0], photons[1]) < 0.45)
      vetoEvent;

    const Particles fs = apply<VetoedFinalState>(event, "calo_fs").particles();
    // Loop over photons and require isolation
    for (const Particle &photon : photons) {
      // Compute calo isolation via particles within an R=0.4 cone of the photon

      FourMomentum mom_in_EtCone;
      for (const Particle &p : fs) {
        // Reject if not in cone
        if (deltaR(photon.momentum(), p.momentum()) > 0.4)
          continue;
        // Sum momentum
        mom_in_EtCone += p.momentum();
      }

      // subtract core photon
      mom_in_EtCone -= photon.momentum();

      // Use photon if energy in isolation cone is low enough
      if (mom_in_EtCone.Et() > 0.045 * photon.momentum().pT() + 6.0 * GeV) {
        vetoEvent;
      }
    }

    const auto &photon = photons.front();

    std::map<string, double> obs;
    obs["pT"] = photon.pT();
    obs["eta"] = photon.eta();
    obs["cos_theta"] = cos(photon.theta());

    const auto &moms = photons.moms();
    const auto total_momentum =
        std::accumulate(moms.begin(), moms.end(), FourMomentum(0, 0, 0, 0));

    obs["inv_m"] = total_momentum.mass();
    obs["o_angle"] = std::abs(tanh((photons[1].eta() - photons[0].eta()) / 2));

    // std::abs(photons[0].theta() + photons[1].theta());

    obs["o_angle_cs"] = std::abs(
        sinh((photons[0].eta() - photons[1].eta())) * 2.0 * photons[0].pT() *
        photons[1].pT() / sqrt(sqr(obs["inv_m"]) + sqr(total_momentum.pT())) /
        obs["inv_m"]);

    obs["total_pT"] = total_momentum.pT();

    for (const auto &name : _observables) {
      _histos[name]->fill(obs[name]);
    }
  }

  //@}

  void finalize() {
    const double sf = crossSection() / (picobarn * sumOfWeights());

    for (auto name : _observables) {
      scale(_histos[name], sf);
    }
  }

  /// @name Histograms
  //@{
  std::map<string, Histo1DPtr> _histos;
  std::vector<string> _observables;
  //@}
};

DECLARE_RIVET_PLUGIN(MC_DIPHOTON_PROTON);

} // namespace Rivet
