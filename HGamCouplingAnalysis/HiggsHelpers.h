#ifndef HIGGSHELPERS_H
#define HIGGSHELPERS_H

// other includes
#include <vector>

namespace HIGGS {
  
  //  Helper function to reweight Powheg+Pythia ggH pT ( NJets < 2, Jet pT > 25 GeV )
  // ------------------------------------------------------------------------------------------
  double ggF_01jet_hpT_weight(double H_pT) {
    if ( H_pT <  20 )  return 1.11;
    if ( H_pT <  45 )  return 1.11 - (H_pT-20)/25*0.2;  // -> 0.91
    if ( H_pT < 135 )  return 0.91 - (H_pT-45)/90*0.36; // -> 0.55
    return 0.55;
  }

  
  //  QCD Scale Uncertainties for ggF 
  // ------------------------------------------------------------------------------------------

  // enum for QCD scale uncertainty source 
  enum ggF_qcdUncSource { yield=1, res=2, cut01=3, cut12=4 };

  // Event weight for propagation of QCD scale uncertainty
  // Input: Number of truth (particle) jets with pT > 30, built excluding the Higgs decay
  //        Number of signma varaiation (+1 for "up", -1 for "down")
  double getJetBinUncertaintyWeight( ggF_qcdUncSource source, int Njets30, double Nsig=+1.0) {

    // Cross sections in the =0, =1, and >=2 jets of Powheg ggH after rewighting scaled to sigma (N3LO)
    static vector<double> yieldUnc({  1.12,  0.66,  0.42 });
    static vector<double>   resUnc({  0.03,  0.57,  0.42 });
    static vector<double> cut01Unc({ -1.22,  1.00,  0.21 });
    static vector<double> cut12Unc({  0.00, -0.86,  0.86 });

    // account for missing EW + quark mass effects by scaling BLPTW total cross section to sigma (N3LO)
    double SF = 48.52/47.4;

    int jetBin = (Njets30 > 1 ? 2 : Njets30);
    if ( source == yield ) return 1.0 + Nsig * yieldUnc[jetBin] / sig[jetBin] * SF;
    if ( source == res   ) return 1.0 + Nsig *   resUnc[jetBin] / sig[jetBin] * SF;
    if ( source == cut01 ) return 1.0 + Nsig * cut01Unc[jetBin] / sig[jetBin] * SF;
    if ( source == cut12 ) return 1.0 + Nsig * cut12Unc[jetBin] / sig[jetBin] * SF;

    // Should never get here, add warning
    return -999.;

  }

}

#endif // HIGGSHELPERS_H

