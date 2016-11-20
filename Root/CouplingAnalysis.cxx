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
  // Loop Over Systematics? 
  m_useSystematics = config()->getBool("HGamCoupling.UseSystematics", false);

  // Only apply the reweighting to ggH MC
  m_reweightHiggsPt = config()->getBool("HGamCoupling.ReweightHiggsPt", false);

  if (isMC()) m_reweightHiggsPt &= getMCSampleName(eventInfo()->mcChannelNumber()).Contains("ggH125");
  else m_reweightHiggsPt = false;

  if (m_reweightHiggsPt) std::cout << "*** !!! REWEIGHTING HIGGS PT !!! ***" << std::endl;
  else               std::cout << "*** !!! NOT REWEIGHTING HIGGS PT !!! ***" << std::endl;

  // Create Histograms
  int nCats(30), nBins(40);

  histoStore()->createTH1F( "h_truthAcc_weightMC", nBins, -0.5, nBins-0.5 );
  histoStore()->createTH1F( "h_truthAcc_weight",   nBins, -0.5, nBins-0.5 );

  TString suffix = ""; 
  for (auto sys: getSystematics()) {
    if (sys.name() != "") {
      if (not m_useSystematics) break;
      if (isData()) break;

      TString sysName = sys.name();
      if (sysName.Contains("Trig")) continue;
      if (sysName.Contains("_CorrUncertainty")) continue;

      suffix = "_" + sys.name();
      suffix.ReplaceAll(" ","_");
    }

    histoStore()->createTH2F( "h2_catSTXS"+suffix,  nCats, 0.5, nCats+0.5, nBins, -0.5, nBins-0.5 );
  }
    
  histoStore()->createTH2F( "h2_catSTXS_QCDyield", nCats, 0.5, nCats+0.5, nBins, -0.5, nBins-0.5 );
  histoStore()->createTH2F( "h2_catSTXS_QCDres",   nCats, 0.5, nCats+0.5, nBins, -0.5, nBins-0.5 );
  histoStore()->createTH2F( "h2_catSTXS_QCDcut01", nCats, 0.5, nCats+0.5, nBins, -0.5, nBins-0.5 );
  histoStore()->createTH2F( "h2_catSTXS_QCDcut12", nCats, 0.5, nCats+0.5, nBins, -0.5, nBins-0.5 );

  for ( int icat(1); icat < nCats; icat++ ) {
    TString histName = TString::Format("h_myy_cat%d",icat);
    histoStore()->createTH1F(histName, 55, 105, 160, ";m_{#gamma#gamma} [GeV];Events / GeV");
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode CouplingAnalysis::execute()
{
  // Important to keep this, so that internal tools / event variables are filled properly.
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

    // Fix a bug in the STXS tool, remove this for new samples
    int   stxsNJ30 = eventInfo()->auxdata<int>("HTXS_Njets_pTjet30");
    double stxsHPT = eventInfo()->auxdata<float>("HTXS_Higgs_pt");
    stage1 = STXS::fixBins_stage1( stage1, stxsNJ30, stxsHPT );
  }


  // Create Histograms for Truth Acceptances ( requires unskimmed samples )
  int STXSbin = STXS::stage1_to_index( stage1 );
  histoStore()->fillTH1F( "h_truthAcc_weightMC", STXSbin, wMC );
  histoStore()->fillTH1F( "h_truthAcc_weight",   STXSbin, w   );
  

  // Loop over systematic variations
  TString suffix = "";
  for (auto sys: getSystematics()) {

    bool nominal = (sys.name() == "");
    if (not nominal) {
      if (not m_useSystematics) break;
      if (isData()) break;

      TString sysName = sys.name();
      if (sysName.Contains("Trig")) continue;
      if (sysName.Contains("_CorrUncertainty")) continue;

      suffix = "_" + sys.name();
      suffix.ReplaceAll(" ","_");
      applySystematicVariation(sys);
    }

    if (not var::isPassed()) return EL::StatusCode::SUCCESS;

    w = (isData()) ? 1.0 : w_pT * weight() * lumiXsecWeight();
    if (w == 0.) return EL::StatusCode::SUCCESS;

    int category = var::catCoup_dev();
    histoStore()->fillTH2F( "h2_catSTXS"+suffix, category, STXSbin, w );

    if (nominal) {
      TString histName = TString::Format("h_myy_cat%d", category);
      histoStore()->fillTH1F( histName, var::m_yy()*HG::invGeV, w );

      int Njets30 = var::N_j_30.truth();

      double wQCDyield = w * HIGGS::getJetBinUncertaintyWeight( HIGGS::yield, Njets30, +1.0);
      double wQCDres   = w * HIGGS::getJetBinUncertaintyWeight( HIGGS::res,   Njets30, +1.0);
      double wQCDcut01 = w * HIGGS::getJetBinUncertaintyWeight( HIGGS::cut01, Njets30, +1.0);
      double wQCDcut12 = w * HIGGS::getJetBinUncertaintyWeight( HIGGS::cut12, Njets30, +1.0);

      histoStore()->fillTH2F( "h2_catSTXS_QCDyield", category, STXSbin, wQCDyield );
      histoStore()->fillTH2F( "h2_catSTXS_QCDres",   category, STXSbin, wQCDres   );
      histoStore()->fillTH2F( "h2_catSTXS_QCDcut01", category, STXSbin, wQCDcut01 );
      histoStore()->fillTH2F( "h2_catSTXS_QCDcut12", category, STXSbin, wQCDcut12 );
    }
    
  }

  return EL::StatusCode::SUCCESS;
}
