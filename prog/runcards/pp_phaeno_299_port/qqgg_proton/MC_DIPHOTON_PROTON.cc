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

struct Observable {
  double _min, _max;
  int _bins = 50;
  bool _log = false;
  double _value = 0;
  Histo1DPtr _hist = nullptr;

  // Observable(double min, double max, int bins = 50, bool log = false)
  //     : _min{min}, _max{max}, _bins{bins}, _log{log} {};

  /// @brief Wraps Histo1D::fill.
  template<typename... Params>
  void fill(double value, Params&&... args) {
    _value = value;
    if (_hist) {
      _hist->fill(value, std::forward<Params>(args)...);
    }
  }
};

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

    double min_pT = 20;
    double eta = 2.5;

    _observables = {{"pT", {min_pT, 1100, 50, true}},
                    {"pT_subl", {min_pT, 1100, 50, true}},
                    {"eta", {-eta, eta}},
                    {"cos_theta", {-1, 1}},
                    {"inv_m", {2 * min_pT - 1, 2000, 50, true}},
                    {"o_angle", {0, 1}},
                    {"o_angle_cs", {0, 1}},
                    {"total_pT", {1e-5, 1100, 50, true}},
                    {"azimuthal_angle", {1e-5, PI, 50, true}}};

    for (auto &[name, observable] : _observables) {
      if (observable._log) {
        book(observable._hist, name,
             logspace(observable._bins, observable._min, observable._max));
        continue;
      }

      book(observable._hist, name, observable._bins, observable._min,
           observable._max);
      book(_h_XS, "xs");
    }
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


    // require an isolation cone around the photons
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

    const auto &photon = photons[0];
    const auto &photon_subl = photons[1];

    _observables.at("pT").fill(photon.pT());
    _observables.at("pT_subl").fill(photon_subl.pT());
    _observables.at("eta").fill(photon.eta());
    _observables.at("cos_theta").fill(cos(photon.theta()));
    _observables.at("azimuthal_angle")
      .fill(PI - mapAngle0ToPi(photon.phi() - photon_subl.phi()));

    const auto &moms = photons.moms();
    const auto total_momentum = moms[0] + moms[1];
    _observables.at("total_pT").fill(total_momentum.pT());
    _observables.at("inv_m").fill(total_momentum.mass());

    _observables.at("o_angle").fill(
        std::abs(tanh((photons[1].eta() - photons[0].eta()) / 2)));

    _observables.at("o_angle_cs").fill(std::abs(
        sinh((photons[0].eta() - photons[1].eta())) * 2.0 * photons[0].pT() *
        photons[1].pT() /
        sqrt(sqr(_observables.at("inv_m")._value) + sqr(total_momentum.pT())) /
        _observables.at("inv_m")._value));
  }

  //@}

  void finalize() {
    _h_XS->addPoint(0, crossSection() / picobarn, .5, crossSectionError() / picobarn);

    const double sf = crossSection() / (picobarn * sumOfWeights());
    for (auto const &[_, observable] : _observables) {
      scale(observable._hist, sf);
    }
  }

  /// @name Histograms
  //@{
  std::map<const std::string, Observable> _observables;
  Scatter2DPtr _h_XS;
  //@}
};

DECLARE_RIVET_PLUGIN(MC_DIPHOTON_PROTON);

} // namespace Rivet
