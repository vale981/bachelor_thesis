// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {
  /// @brief MC validation analysis for isolated di-photon events
  class MC_DIPHOTON : public Analysis {
  public:

    /// Constructor
    MC_DIPHOTON()
      : Analysis {"MC_DIPHOTON"}
    {    }


    /// @name Analysis methods
    //@{

    void init() {
      FinalState fs;
      declare(fs, "FS");

      IdentifiedFinalState ifs(Cuts::abseta < 2 && Cuts::pT > 20*GeV);
      ifs.acceptId(PID::PHOTON);
      declare(ifs, "IFS");

      book(_h_m_PP ,"m_PP", logspace(50, 1.0, 0.25*(sqrtS()>0.?sqrtS():14000.)));
      book(_h_pT_PP ,"pT_PP", logspace(50, 1.0, 0.25*(sqrtS()>0.?sqrtS():14000.)));
      book(_h_pT_P1 ,"pT_P1", 50, 0.0, 70.0);
      book(_h_pT_P2 ,"pT_P2", 50, 0.0, 70.0);
      book(_h_dphi_PP ,"dphi_PP", 20, 0.0, M_PI);
    }


    void analyze(const Event& event) {
      const double weight = 1.0;

      Particles photons = apply<IdentifiedFinalState>(event, "IFS").particles();

      if (photons.size() < 2) {
        vetoEvent;
      }

      // Isolate photons with ET_sum in cone
      Particles isolated_photons;
      Particles fs = apply<FinalState>(event, "FS").particlesByPt();
      for (const Particle& photon : photons) {
        FourMomentum mom_in_cone;
        double eta_P = photon.eta();
        double phi_P = photon.phi();
        for (const Particle& p : fs) {
          if (deltaR(eta_P, phi_P, p.eta(), p.phi()) < 0.4) {
            mom_in_cone += p.momentum();
          }
        }
        if (mom_in_cone.Et()-photon.Et() < 4.0*GeV) {
          isolated_photons.push_back(photon);
        }
      }

      if (isolated_photons.size() != 2) {
        vetoEvent;
      }

      _h_pT_P1->fill(isolated_photons[0].pT(), weight);
      _h_pT_P2->fill(isolated_photons[1].pT(), weight);
      FourMomentum mom_PP = isolated_photons[0].momentum() + isolated_photons[1].momentum();
      _h_m_PP->fill(mom_PP.mass(), weight);
      _h_pT_PP->fill(mom_PP.pT(), weight);
      _h_dphi_PP->fill(deltaPhi(isolated_photons[0].phi(),
                                isolated_photons[1].phi()), weight);
    }


    void finalize() {
      scale(_h_m_PP, crossSection()/sumOfWeights());
      scale(_h_pT_PP, crossSection()/sumOfWeights());
      scale(_h_pT_P1, crossSection()/sumOfWeights());
      scale(_h_pT_P2, crossSection()/sumOfWeights());
      scale(_h_dphi_PP, crossSection()/sumOfWeights());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_m_PP;
    Histo1DPtr _h_pT_PP;
    Histo1DPtr _h_pT_P1;
    Histo1DPtr _h_pT_P2;
    Histo1DPtr _h_dphi_PP;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_DIPHOTON);

}
