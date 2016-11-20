#include "HGamCouplingAnalysis/CouplingAnalysis.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamAnalysisFramework/HGamVariables.h"

// Additional Helpers
#include "HGamCouplingAnalysis/STXSHelpers.h"
#include "HGamCouplingAnalysis/HiggsHelpers.h"


// this is needed to distribute the algorithm to the workers
ClassImp(CouplingAnalysis)



CouplingAnalysis::CouplingAnalysis(const char *name)
: HgammaAnalysis(name)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



CouplingAnalysis::~CouplingAnalysis()
{
  // Here you delete any memory you allocated during your analysis.
}



EL::StatusCode CouplingAnalysis::createOutput()
{
  // Only apply the reweighting to ggH MC
  m_reweightHiggsPt = config()->getBool("STXS.ReweightHiggsPt", false);

  if (isMC()) m_reweightHiggsPt &= getMCSampleName(eventInfo()->mcChannelNumber()).Contains("ggH125");
  else m_reweightHiggsPt = false;

  if (m_reweightHiggsPt) std::cout << "*** !!! REWEIGHTING HIGGS PT !!! ***" << std::endl;

  // Create Histograms
  histoStore()->createTH2F( "h2_catICHEP_tbin", 40, -0.5, 39.5, 13, 0.5, 13.5 );
  histoStore()->createTH2F( "h2_catSTXS_tbin",  40, -0.5, 39.5, 30, 0.5, 30.5 );
  histoStore()->createTH2F( "h2_catSTXS_HCombBin", 53, -0.5, 52.5, 30, 0.5, 30.5 );
  histoStore()->createTH1F( "h_truthAcc_weightMC", 40, -0.5, 39.5 );
  histoStore()->createTH1F( "h_truthAcc_weight",   40, -0.5, 39.5 );

  for ( int icat(1); icat < 30; icat++ ) {
    TString histName = TString::Format("h_cat%d_myy",icat);
    histoStore()->createTH1F(histName, 55, 105, 160, ";m_{#gamma#gamma} [GeV];Events / GeV");
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode CouplingAnalysis::execute()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  // Important to keep this, so that internal tools / event variables
  // are filled properly.
  HgammaAnalysis::execute();

  // Blind the data
  if (isData() && var::m_yy() > 120*HG::GeV && var::m_yy() < 130*HG::GeV) return EL::StatusCode::SUCCESS;
  
  // Apply Higgs pT reweighting
  double w_pT = 1.0;
  if ( isMC() && m_reweightHiggsPt && (var::N_j.truth() < 2) ) {
    w_pT = HIGGS::ggF_01jet_hpT_weight( var::pT_h1.truth()*HG::invGeV );
  }

  double wMC = (isData()) ? 1.0 : w_pT * eventHandler()->mcWeight();
  double w   = (isData()) ? 1.0 : w_pT * weight() * lumiXsecWeight();

  
  // Save Histogram for Truth Acceptance
  int stage1(0), prodMode(0), errorCode(0);
  if (isMC() && eventInfo()->isAvailable<int>("HTXS_Stage1_Category_pTjet30")) {
    stage1   = eventInfo()->auxdata<int>("HTXS_Stage1_Category_pTjet30");
    prodMode = eventInfo()->auxdata<int>("HTXS_prodMode");

    errorCode = eventInfo()->auxdata<int>("HTXS_errorCode");
    if (errorCode != 0) stage1 = -1;

    // Fix a bug in the STXS tool
    if (stage1 == 101 || stage1 == 102) {
      if (eventInfo()->auxdata<float>("HTXS_Higgs_pt") > 200 ) {
        int stxsNJ30 = eventInfo()->auxdata<int>("HTXS_Njets_pTjet30");
        if (stxsNJ30 == 1) stage1 = 107;
        if (stxsNJ30 == 2) stage1 = 111;
      }
    }
  }

  int STXSbin = STXS::stage1_to_index( stage1 );
  
  // Create Histograms for Truth Acceptances ( requires unskimmed samples )
  histoStore()->fillTH1F( "h_truthAcc_weightMC", STXSbin, wMC );
  histoStore()->fillTH1F( "h_truthAcc_weight",   STXSbin, w   );
  
  //int mcID = (isData()) ? 0 : eventInfo()->mcChannelNumber();
  //bool isTWH = (341997 <= mcID && mcID <= 341999);
  //int stage1_Hcomb = HTXSstage1_to_Hcomb( stage1, prodMode, isTWH );
  //int hbin = HComb_to_index( stage1_Hcomb );

  if (w == 0.) return EL::StatusCode::SUCCESS;



  return EL::StatusCode::SUCCESS;
}
