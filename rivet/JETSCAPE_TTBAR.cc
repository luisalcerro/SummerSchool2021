#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// @brief Add a short analysis description here
  class JETSCAPE_TTBAR : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(JETSCAPE_TTBAR);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // A FinalState is used to select particles within |eta| < 4.2 and with pT
      // > 30 GeV, out of which the ChargedLeptons projection picks only the muons, to be accessed later as "LFS".
      ChargedLeptons lfs(FinalState(Cuts::abspid == PID::MUON && Cuts::abseta < 4.2 && Cuts::pT > 30*GeV));
      declare(lfs, "LFS");

      // A second FinalState is used to select all particles in |eta| < 4.2,
      // with no pT cut. This is used to construct jets and measure missing
      // transverse energy.
      VetoedFinalState fs(FinalState(Cuts::abseta < 4.2));
      fs.addVetoOnThisFinalState(lfs);
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "Jets");
      declare(MissingMomentum(fs), "MissingET");

      //booking histograms for leptons
      book(_h["nLep"], "lep_mult", 11, -0.5, 10.5);
      book(_h["lepPt"], "lep_pT", 200, 0.0, 700.0);
      book(_h["lepEta"], "lep_eta", 50, -5.0, 5.0);
      book(_h["lepPhi"], "lep_phi", 20, -0.5, 7.0);

      //book(_h["lep_1_Pt"], "lep1_pT", 200, 0.0, 700.0);
      //book(_h["lep_2_Pt"], "lep2_pT", 200, 0.0, 700.0);
      book(_h["miss_E"], "miss_E", 200, 0.0, 700.0);
      book(_h["W_mass"], "W_mass", 75, 30, 180);
      //booking histograms for Jets
      book(_h["njets"], "jet_mult", 11, -0.5, 10.5);
      book(_h["jetPt"], "jet_pT", 100, 0.0, 500.);
      book(_h["jetEta"], "jet_eta", 50, -5.0, 5.0);
      book(_h["jetPhi"], "jet_phi", 20, -0.5, 7.0);
      book(_h["nBjets"], "bjet_mult", 11, -0.5, 10.5);
      book(_h["nLjets"], "ljet_mult", 11, -0.5, 10.5);
    }




    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Use the "LFS" projection to require at least one hard charged
      // lepton. This is an experimental signature for the leptonically decaying
      // W. This helps to reduce pure QCD backgrounds.
      const ChargedLeptons& lfs = apply<ChargedLeptons>(event, "LFS");
      MSG_DEBUG("Charged lepton multiplicity = " << lfs.chargedLeptons().size());
      for (const Particle& lepton : lfs.chargedLeptons()) {
        MSG_DEBUG("Lepton pT = " << lepton.pT());
      }

      size_t nLeps = lfs.chargedLeptons().size();
      if (nLeps < 1){
        MSG_DEBUG("Event failed lepton multiplicity cut");
        vetoEvent;
      }


      // Use a missing ET cut to bias toward events with a hard neutrino from
      // the leptonically decaying W. This helps to reduce pure QCD backgrounds.
      // not applied in all-hadronic mode
      const Vector3& met = apply<MissingMomentum>(event, "MissingET").vectorMissingPt();
      MSG_DEBUG("Vector pT = " << met.mod() << " GeV");
      if (met.mod() < 30*GeV) {
        MSG_DEBUG("Event failed missing ET cut");
        vetoEvent;
      }

      // Use the "Jets" projection to check how many jets with pT > 30 GeV there are
      // remove jets overlapping with any lepton (dR < 0.3)
      // cut on jet multiplicity depending on ttbar decay mode
      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      const Jets jets = discardIfAnyDeltaRLess(jetpro.jetsByPt(30*GeV), lfs.chargedLeptons(), 0.3);

      if ( jets.size() < 4)  vetoEvent;

      // Fill lepton histograms
      _h["nLep"]->fill(lfs.chargedLeptons().size());
      for (const Particle& lepton : lfs.chargedLeptons()) {
        _h["lepPt"]->fill(lepton.pT());
        _h["lepEta"]->fill(lepton.eta());
        _h["lepPhi"]->fill(lepton.phi());

      }

      // Fill MET histogram
      _h["miss_E"]->fill(met.mod());



      // Fill All jets histograms
      _h["njets"]->fill(jets.size());
      for (const Jet& j : jets) {
        _h["jetPt"]->fill( j.pT());
        _h["jetEta"]->fill(j.eta());
        _h["jetPhi"]->fill(j.phi());
      }

      // distinguish heavy and light flavor jets
      Jets bjets, ljets;
      for (const Jet& jet : jets) {
        if (jet.bTagged())  bjets += jet;
        else                ljets += jet;
      }

      _h["nBjets"]->fill(bjets.size());
      _h["nLjets"]->fill(ljets.size());


      // Construct the hadronically decaying W momentum 4-vector from pairs of
      // non-b-tagged jets. The pair which best matches the W mass is used. We start
      // with an always terrible 4-vector estimate which should always be "beaten" by
      // a real jet pair.
      FourMomentum W(10*(sqrtS()>0.?sqrtS():14000.), 0, 0, 0);
      for (size_t i = 0; i < ljets.size()-1; ++i) {
        for (size_t j = i + 1; j < ljets.size(); ++j) {
          const FourMomentum Wcand = ljets[i].momentum() + ljets[j].momentum();
          MSG_TRACE(i << "," << j << ": candidate W mass = " << Wcand.mass()/GeV
                    << " GeV, vs. incumbent candidate with " << W.mass()/GeV << " GeV");
          if (fabs(Wcand.mass() - 80.4*GeV) < fabs(W.mass() - 80.4*GeV)) {
            W = Wcand;
          }
        }
      }
      MSG_DEBUG("Candidate W mass = " << W.mass() << " GeV");
      _h["W_mass"]->fill(W.mass());

      }

    /// Normalise histograms etc., after the run
    void finalize() {

   //   const double sf = crossSection() / sumOfWeights();
   //   scale(_h, sf);

    }

    //@}

    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    size_t _mode;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(JETSCAPE_TTBAR);
}
