#ifdef RIVET_DECLARE_PLUGIN

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"

namespace Rivet {


  /// @brief Semivisible jets t-channel search
  class ATLAS_2023_I2663256 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2023_I2663256);


  private:

    /// @name Analysis methods
    /// @{

    /// Initialise and register projections
    void init() {

      //Leptons
      ChargedLeptons lfs(FinalState(Cuts::abseta < 2.5 && Cuts::pT > 7*GeV));
      declare(lfs, "LFS");

      FinalState fs(Cuts::abseta < 4.5);
      VetoedFinalState vfs(fs);
      vfs.vetoNeutrinos();
      vfs.addVetoPairId(PID::MUON);

      FastJets j04(vfs, JetAlg::ANTIKT, 0.4);
      declare(j04, "Jets");

      SmearedJets sj04(j04, JET_SMEAR_ATLAS_RUN2, [](const Jet& j){
        return j.bTagged() ? 0.7*(1 - exp(-j.pT()/(10*GeV))) : 0.01;
      });
      declare(sj04, "SJets");

      // Book histograms

      book(_h["ninebin"], 1, 1, 8);
      book(_h["ht"],      2, 1, 8);
      book(_h["met"],     3, 1, 8);
      book(_h["ptbal"],   4, 1, 8);
      book(_h["difphi"],  5, 1, 8);

    }



    // Perform the per-event analysis
    void analyze(const Event& event) {

      const ChargedLeptons& lfs = apply<ChargedLeptons>(event, "LFS");
      if (lfs.chargedLeptons().size() > 0) vetoEvent;

      const Jets j04 = apply<JetFinder>(event, "Jets").jetsByPt(Cuts::abseta < 2.8 && Cuts::pT > 30*GeV);
      const Jets sj04 = apply<JetFinder>(event, "SJets").jetsByPt(Cuts::pT > 30*GeV);
      if(sj04.size() < 2) vetoEvent;
      const Jets sj04b = select(sj04, [&](const Jet& j) { return j.bTagged(); });
      if(sj04b.size() > 1) vetoEvent;
      if( sj04[0].pT() < 250) vetoEvent;

      FourMomentum pTmiss;
      for (const Jet& j : sj04) {
        pTmiss -= j.momentum();
      }


      double met = pTmiss.perp();
      if(met < 600) vetoEvent;


      double sumpt = 0;
      double minphi = 99;
      double maxphi = -99;

      Jet svj;
      Jet antisvj;


      for (const Jet& jet : sj04) {
        sumpt += jet.pT();
        if (deltaPhi(jet, pTmiss) < minphi) {
          minphi = deltaPhi(jet, pTmiss);
          svj=jet;
        }

        if (deltaPhi(jet, pTmiss) > maxphi) {
          maxphi = deltaPhi(jet, pTmiss);
          antisvj=jet;
        }
      }


      if (sumpt < 600) vetoEvent;
      if (minphi > 2) vetoEvent;
      double difphi = deltaPhi(svj, antisvj);
      double delta_jets = (svj.momentum() + antisvj.momentum()).pT();
      double total_pt = svj.momentum().pT() + antisvj.momentum().pT();
      double delta_jets_n = delta_jets/ total_pt;

      _h["met"]->fill(met);
      _h["ht"]->fill(sumpt);
      _h["difphi"]->fill(difphi);
      _h["ptbal"]->fill(delta_jets_n);

      int i=0;
      int j=0;
      if(inRange(delta_jets_n, 0.0, 0.6)) i=0;
      if(inRange(delta_jets_n, 0.6, 0.9)) i=1;
      if(inRange(delta_jets_n, 0.9, 1.0)) i=2;
      if(inRange(difphi, 0.0, 2.0)) j=0;
      if(inRange(difphi, 2.0, 2.7)) j=1;
      if(inRange(difphi, 2.7, 3.2)) j=2;
      int binindex = 3*i + j + 1;
      _h["ninebin"]->fill(binindex);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      for (auto& hist : _h) {
        for (size_t i=0; i < hist.second->numBins(); ++i) {
          double bW = hist.second->bin(i).xWidth();
          hist.second->bin(i).scaleW(bW);
        }
      }

      double norm = 139*crossSection()/femtobarn/sumOfWeights();
      scale(_h, norm);
    }

    /// @}


  private:

    /// Histograms
    map<string, Histo1DPtr> _h;

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2023_I2663256);

}

#endif
