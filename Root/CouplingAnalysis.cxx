#include "HGamCouplingAnalysis/CouplingAnalysis.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamAnalysisFramework/HGamVariables.h"

// Additional Helpers
#include "HGamCouplingAnalysis/STXSHelpers.h"
#include "HGamCouplingAnalysis/HiggsHelpers.h"
#include <cmath>


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
  
  int mcid = (isData()) ? -999 : eventInfo()->mcChannelNumber();

  m_isGGH = false;
  m_isTWH = false;
  if (isMC()) {
    m_isGGH = getMCSampleName(mcid).Contains("ggH125");
    m_isTWH = getMCSampleName(mcid).Contains("tWH125");
    m_reweightHiggsPt &= m_isGGH;
  } else m_reweightHiggsPt = false;
    

  if (m_reweightHiggsPt) std::cout << "*** !!! REWEIGHTING HIGGS PT !!! ***" << std::endl;
  else               std::cout << "*** !!! NOT REWEIGHTING HIGGS PT !!! ***" << std::endl;
      
  m_usePDFUncerts = false;
  if (isMC()) {
    StrV higgsTypes = config()->getStrV("EventHandler.HiggsWeights.Types", {});
    for (TString sample: higgsTypes) {
      StrV dsids = config()->getStrV("EventHandler.HiggsWeights."+sample, {""});
      for (TString dsid: dsids) {
        if (mcid == std::atoi(dsid.Data())) {
          m_usePDFUncerts = true;
        }
      }
    }
  }
  m_usePDFUncerts &= config()->getBool("HGamCoupling.UsePDFSystematics", false);

  // Create Data TTree
  if (isData()) {
    m_tree = new TTree("physics","physics");
    wk()->addOutput(m_tree);
    m_tree->Branch( "category", &m_category );
    m_tree->Branch( "m_yy",     &m_myy );
  }

  // Create Histograms
  int nCats(33), nBins(42), nIndex(53);

  histoStore()->createTH1F( "h_truthAcc_fineIndex_weightMC", nIndex, -0.5, nIndex-0.5 );
  histoStore()->createTH1F( "h_truthAcc_fineIndex_weight", nIndex, -0.5, nIndex-0.5 );

  histoStore()->createTH1F( "h_truthAcc_weightMC", nBins, -0.5, nBins-0.5 );
  histoStore()->createTH1F( "h_truthAcc_weight", nBins, -0.5, nBins-0.5 );

  histoStore()->createTH1F( "h_initialEvents_weighted",      1, -10, 10 );
  histoStore()->createTH1F( "h_initialEvents_weighted_pTRW", 1, -10, 10 ); 

  histoStore()->createTH1F( "h_initialEvents_noDalitz_weighted",      1, -10, 10 );
  histoStore()->createTH1F( "h_initialEvents_noDalitz_weighted_pTRW", 1, -10, 10 ); 

  HG::SystematicList allSys = getSystematics();
  int Nsysts = (isData()) ? 0 : allSys.size()-1;
  if (Nsysts > 0) std::cout << "Running over " << Nsysts << " systematics:" << std::endl;
  if (Nsysts > 0) {
    for (int i(0); i <= Nsysts; i++ )
      std::cout << "\t" << i << " " << allSys[i].name() << std::endl;
  }
  std::cout << std::endl;

  static int sysIndex = config()->getInt("HGamCoupling.SystematicIndex", -1);
  if (sysIndex > Nsysts) 
    std::cout << "Canceling job for syst index " << sysIndex 
              << ", there are only " << Nsysts << " variations." << std::endl;

  if (sysIndex < 0) sysList = allSys;
  else sysList.push_back( allSys[sysIndex] );

  if (Nsysts > 0) std::cout << "Running over " << Nsysts << " systematics:" << std::endl;

  TString suffix = ""; 
  for (auto sys: sysList) {
    if (sys.name() != "") {
      if (not m_useSystematics) break;
      if (isData()) break;

      TString sysName = sys.name();
      //if (sysName.Contains("Trig")) continue;
      //if (sysName.Contains("_CorrUncertainty")) continue;
      
      if (Nsysts > 0) std::cout << "\t" << sysName << std::endl;

      suffix = "_" + sys.name();
      suffix.ReplaceAll(" ","_");
    }

    histoStore()->createTH1F( "h_inclusive"+suffix,  1, -10, 10 );
    histoStore()->createTH1F(   "h_pT_y1", 50, 0, 250 );
    histoStore()->createTH1F(   "h_pT_y2", 50, 0, 250 );
    histoStore()->createTH1F(  "h_eta_y1", 52, -2.6, 2.6 );
    histoStore()->createTH1F(  "h_eta_y2", 52, -2.6, 2.6 );
 
    histoStore()->createTH1F(  "h_catSTXS"+suffix,  nCats, 0.5, nCats+0.5 );
    histoStore()->createTH2F( "h2_catSTXS"+suffix,  nCats, 1, nCats+1, nBins, 0, nBins );
    histoStore()->createTH2F( "h2_fineIndex"+suffix,  nCats, 1, nCats+1, nIndex, 0, nIndex );
  }

  // Histograms below only for nominal samples
  if (sysIndex > 0) return EL::StatusCode::SUCCESS;
    
  if (m_isGGH) {
    histoStore()->createTH2F( "h2_catSTXS_QCDyield", nCats, 0.5, nCats+0.5, nBins, -0.5, nBins-0.5 );
    histoStore()->createTH2F( "h2_catSTXS_QCDres",   nCats, 0.5, nCats+0.5, nBins, -0.5, nBins-0.5 );
    histoStore()->createTH2F( "h2_catSTXS_QCDcut01", nCats, 0.5, nCats+0.5, nBins, -0.5, nBins-0.5 );
    histoStore()->createTH2F( "h2_catSTXS_QCDcut12", nCats, 0.5, nCats+0.5, nBins, -0.5, nBins-0.5 );
    
    histoStore()->createTH2F( "h2_fineIndex_QCDyield", nCats, 0.5, nCats+0.5, nIndex, -0.5, nIndex-0.5 );
    histoStore()->createTH2F( "h2_fineIndex_QCDres",   nCats, 0.5, nCats+0.5, nIndex, -0.5, nIndex-0.5 );
    histoStore()->createTH2F( "h2_fineIndex_QCDcut01", nCats, 0.5, nCats+0.5, nIndex, -0.5, nIndex-0.5 );
    histoStore()->createTH2F( "h2_fineIndex_QCDcut12", nCats, 0.5, nCats+0.5, nIndex, -0.5, nIndex-0.5 );
    
    histoStore()->createTH2F( "h2_catSTXS_WG1_QCDmu",    nCats, 0.5, nCats+0.5, nBins, -0.5, nBins-0.5 );
    histoStore()->createTH2F( "h2_catSTXS_WG1_QCDqm",    nCats, 0.5, nCats+0.5, nBins, -0.5, nBins-0.5 );
    histoStore()->createTH2F( "h2_catSTXS_WG1_QCDres",   nCats, 0.5, nCats+0.5, nBins, -0.5, nBins-0.5 );
    histoStore()->createTH2F( "h2_catSTXS_WG1_QCDpTH",   nCats, 0.5, nCats+0.5, nBins, -0.5, nBins-0.5 );
    histoStore()->createTH2F( "h2_catSTXS_WG1_QCDmig01", nCats, 0.5, nCats+0.5, nBins, -0.5, nBins-0.5 );
    histoStore()->createTH2F( "h2_catSTXS_WG1_QCDmig12", nCats, 0.5, nCats+0.5, nBins, -0.5, nBins-0.5 );
    
    histoStore()->createTH2F( "h2_fineIndex_WG1_QCDmu",    nCats, 0.5, nCats+0.5, nIndex, -0.5, nIndex-0.5 );
    histoStore()->createTH2F( "h2_fineIndex_WG1_QCDqm",    nCats, 0.5, nCats+0.5, nIndex, -0.5, nIndex-0.5 );
    histoStore()->createTH2F( "h2_fineIndex_WG1_QCDres",   nCats, 0.5, nCats+0.5, nIndex, -0.5, nIndex-0.5 );
    histoStore()->createTH2F( "h2_fineIndex_WG1_QCDpTH",   nCats, 0.5, nCats+0.5, nIndex, -0.5, nIndex-0.5 );
    histoStore()->createTH2F( "h2_fineIndex_WG1_QCDmig01", nCats, 0.5, nCats+0.5, nIndex, -0.5, nIndex-0.5 );
    histoStore()->createTH2F( "h2_fineIndex_WG1_QCDmig12", nCats, 0.5, nCats+0.5, nIndex, -0.5, nIndex-0.5 );
  }

  if (m_usePDFUncerts) {
    for (TString suffix: {"_alphaS_up","_alphaS_dn"}) {
      histoStore()->createTH1F( "h_catSTXS"+suffix,  nCats, 0.5, nCats+0.5 );
      histoStore()->createTH2F( "h2_catSTXS"+suffix,  nCats, 1, nCats+1, nBins, 0, nBins );
      histoStore()->createTH2F( "h2_fineIndex"+suffix,  nCats, 1, nCats+1, nIndex, 0, nIndex );
    }

    for (int ipdf(0); ipdf < 30; ipdf++) {
      TString suffixPDF = TString::Format("_PDF%d",ipdf);
      histoStore()->createTH1F( "h_catSTXS"+suffixPDF,  nCats, 0.5, nCats+0.5 );
      histoStore()->createTH2F( "h2_catSTXS"+suffixPDF,  nCats, 1, nCats+1, nBins, 0, nBins );
      histoStore()->createTH2F( "h2_fineIndex"+suffixPDF,  nCats, 1, nCats+1, nIndex, 0, nIndex );
    }
  }

  for ( int icat(1); icat < nCats+1; icat++ ) {
    TString histName = TString::Format("h_myy_cat%d",icat);
    histoStore()->createTH1F( histName, 110, 105, 160, ";m_{#gamma#gamma} [GeV];Events / GeV");
  }

  return EL::StatusCode::SUCCESS;
}




EL::StatusCode CouplingAnalysis::execute()
{
  // Important to keep this, so that internal tools / event variables are filled properly.
  HgammaAnalysis::execute();

  // Blind the data
  if (isData() && var::m_yy() >= 120*HG::GeV && var::m_yy() < 130*HG::GeV) return EL::StatusCode::SUCCESS;

  
  // Apply Higgs pT reweighting
  double w_pT = 1.0; 
  if ( isMC() && m_reweightHiggsPt && (eventInfo()->auxdata<int>("HTXS_Njets_pTjet25") < 2) ) {
    w_pT = HIGGS::ggF_01jet_hpT_weight( var::pT_h1.truth()*HG::invGeV );
  }

  // get correction factor for N_init with pT reweighting
  double corrDenom = ( isMC() && m_reweightHiggsPt ) ? 1.0 : 1.003;
  
  double wMC = (isData()) ? 1.0 : w_pT * eventHandler()->mcWeight();
  double w  = (isData()) ? 1.0 : w_pT * weightCatCoup_Moriond2017BDT() * lumiXsecWeight();
  //double w  = (isData()) ? 1.0 : w_pT * weightCatCoup_Moriond2017() * lumiXsecWeight();

  // Save Histogram for Truth Acceptance
  int stage1(0), errorCode(0);
  if (isMC() && eventInfo()->isAvailable<int>("HTXS_Stage1_Category_pTjet30")) {
    stage1   = eventInfo()->auxdata<int>("HTXS_Stage1_Category_pTjet30");
    errorCode = eventInfo()->auxdata<int>("HTXS_errorCode");
    if (errorCode != 0) stage1 = -1;
  }
  
  int fineIndex(0);
  if (isMC() && eventInfo()->isAvailable<int>("HTXS_Stage1_FineIndex_pTjet30")) {
    fineIndex = eventInfo()->auxdata<int>("HTXS_Stage1_FineIndex_pTjet30");
    if (stage1 / 100 == 8) {
      fineIndex = 49 + (stage1%100);
      if (m_isTWH) fineIndex += 2;
    }
  }

  // Create Histograms for Truth Acceptances ( requires unskimmed samples )
  int STXSbin = STXS::stage1_to_index( stage1 );
  if (m_isTWH && STXSbin > 0) STXSbin += 2;

  histoStore()->fillTH1F( "h_truthAcc_fineIndex_weightMC", fineIndex, wMC );
  histoStore()->fillTH1F( "h_truthAcc_fineIndex_weight", fineIndex, w );

  histoStore()->fillTH1F( "h_truthAcc_weightMC", STXSbin, wMC );
  histoStore()->fillTH1F( "h_truthAcc_weight", STXSbin, w );

  // Create histos for weighting
  if ( isMC() ) {
    histoStore()->fillTH1F( "h_initialEvents_weighted_pTRW", 1.0, weightInitial() * w_pT );
    histoStore()->fillTH1F( "h_initialEvents_weighted",      1.0, weightInitial() );
  }

  if ( isMC() && not var::isDalitzEvent() ) {
    histoStore()->fillTH1F( "h_initialEvents_noDalitz_weighted_pTRW", 1.0, weightInitial() * w_pT );
    histoStore()->fillTH1F( "h_initialEvents_noDalitz_weighted",      1.0, weightInitial() );
  }

  // Loop over systematic variations
  TString suffix = "";
  for (auto sys: sysList) {

    bool nominal = (sys.name() == "");
    if (not nominal) {
      if (not m_useSystematics) break;
      if (isData()) break;

      TString sysName = sys.name();

      suffix = "_" + sys.name();
      suffix.ReplaceAll(" ","_");
      applySystematicVariation(sys);
    }

    if (not var::isPassed()) return EL::StatusCode::SUCCESS;

    //w = (isData()) ? 1.0 : w_pT * corrDenom * weightCatCoup_Moriond2017() * lumiXsecWeight();
    w = (isData()) ? 1.0 : w_pT * corrDenom * weightCatCoup_Moriond2017BDT() * lumiXsecWeight();
    if (w == 0.) return EL::StatusCode::SUCCESS;
    
    histoStore()->fillTH1F( "h_inclusive"+suffix, 0, w );

    //m_category = var::catCoup_Moriond2017();
    m_category = var::catCoup_Moriond2017BDT();

    histoStore()->fillTH1F(  "h_catSTXS"+suffix,   m_category, w );
    histoStore()->fillTH2F( "h2_catSTXS"+suffix,   m_category, STXSbin, w );
    histoStore()->fillTH2F( "h2_fineIndex"+suffix, m_category, fineIndex, w );

    if (nominal) {
      TString histName = TString::Format("h_myy_cat%d", m_category);
      histoStore()->fillTH1F( histName, var::m_yy()*HG::invGeV, w );

      xAOD::PhotonContainer all_photons = photonHandler()->getCorrectedContainer();
      xAOD::PhotonContainer selPhotons = photonHandler()->applySelection(all_photons);
      histoStore()->fillTH1F( "h_pT_y1", var::pT_y1()*HG::invGeV, w );
      histoStore()->fillTH1F( "h_pT_y2", var::pT_y2()*HG::invGeV, w );
      histoStore()->fillTH1F( "h_eta_y1", selPhotons[0]->eta(), w );
      histoStore()->fillTH1F( "h_eta_y2", selPhotons[1]->eta(), w );

      m_myy = var::m_yy()*HG::invGeV;
      if (isData()) m_tree->Fill();

      if (m_isGGH) {
        int Njets30 = eventInfo()->auxdata<int>("HTXS_Njets_pTjet30");
        //int Njets30 = var::N_j_30.truth();
        
        {
          double wQCDyield = w * HIGGS::getJetBinUncertaintyWeight( HIGGS::yield, Njets30, +1.0);
          double wQCDres   = w * HIGGS::getJetBinUncertaintyWeight( HIGGS::res,   Njets30, +1.0);
          double wQCDcut01 = w * HIGGS::getJetBinUncertaintyWeight( HIGGS::cut01, Njets30, +1.0);
          double wQCDcut12 = w * HIGGS::getJetBinUncertaintyWeight( HIGGS::cut12, Njets30, +1.0);
          
          histoStore()->fillTH2F( "h2_catSTXS_QCDyield", m_category, STXSbin, wQCDyield );
          histoStore()->fillTH2F( "h2_catSTXS_QCDres",   m_category, STXSbin, wQCDres   );
          histoStore()->fillTH2F( "h2_catSTXS_QCDcut01", m_category, STXSbin, wQCDcut01 );
          histoStore()->fillTH2F( "h2_catSTXS_QCDcut12", m_category, STXSbin, wQCDcut12 );
          
          histoStore()->fillTH2F( "h2_fineIndex_QCDyield", m_category, fineIndex, wQCDyield );
          histoStore()->fillTH2F( "h2_fineIndex_QCDres",   m_category, fineIndex, wQCDres   );
          histoStore()->fillTH2F( "h2_fineIndex_QCDcut01", m_category, fineIndex, wQCDcut01 );
          histoStore()->fillTH2F( "h2_fineIndex_QCDcut12", m_category, fineIndex, wQCDcut12 );
        }
        
        {
          xAOD::HiggsWeights higgsWeights = eventHandler()->higgsWeights();

          double wQCDmu    = w * higgsWeights.qcd_wg1_mu    / higgsWeights.nominal;
          double wQCDqm    = w * higgsWeights.qcd_wg1_qm    / higgsWeights.nominal;
          double wQCDres   = w * higgsWeights.qcd_wg1_res   / higgsWeights.nominal;
          double wQCDpTH   = w * higgsWeights.qcd_wg1_pTH   / higgsWeights.nominal;
          double wQCDmig01 = w * higgsWeights.qcd_wg1_mig01 / higgsWeights.nominal;
          double wQCDmig12 = w * higgsWeights.qcd_wg1_mig12 / higgsWeights.nominal;
          
          histoStore()->fillTH2F( "h2_catSTXS_WG1_QCDmu",    m_category, STXSbin, wQCDmu    );
          histoStore()->fillTH2F( "h2_catSTXS_WG1_QCDqm",    m_category, STXSbin, wQCDqm    );
          histoStore()->fillTH2F( "h2_catSTXS_WG1_QCDres",   m_category, STXSbin, wQCDres   );
          histoStore()->fillTH2F( "h2_catSTXS_WG1_QCDpTH",   m_category, STXSbin, wQCDpTH   );
          histoStore()->fillTH2F( "h2_catSTXS_WG1_QCDmig01", m_category, STXSbin, wQCDmig01 );
          histoStore()->fillTH2F( "h2_catSTXS_WG1_QCDmig12", m_category, STXSbin, wQCDmig12 );
          
          histoStore()->fillTH2F( "h2_fineIndex_WG1_QCDmu",    m_category, fineIndex, wQCDmu    );
          histoStore()->fillTH2F( "h2_fineIndex_WG1_QCDqm",    m_category, fineIndex, wQCDqm    );
          histoStore()->fillTH2F( "h2_fineIndex_WG1_QCDres",   m_category, fineIndex, wQCDres   );
          histoStore()->fillTH2F( "h2_fineIndex_WG1_QCDpTH",   m_category, fineIndex, wQCDpTH   );
          histoStore()->fillTH2F( "h2_fineIndex_WG1_QCDmig01", m_category, fineIndex, wQCDmig01 );
          histoStore()->fillTH2F( "h2_fineIndex_WG1_QCDmig12", m_category, fineIndex, wQCDmig12 );
        }
    
      }

      // PDF and alpha_S uncertainties
      if (isMC() && m_usePDFUncerts) {
        xAOD::HiggsWeights higgsWeights = eventHandler()->higgsWeights();

        double wASHI = w * higgsWeights.alphaS_up / higgsWeights.nominal;
        if (!std::isfinite(wASHI)) wASHI = w;
        histoStore()->fillTH1F(  "h_catSTXS_alphaS_up",   m_category, wASHI );
        histoStore()->fillTH2F( "h2_catSTXS_alphaS_up",   m_category, STXSbin, wASHI );
        histoStore()->fillTH2F( "h2_fineIndex_alphaS_up", m_category, fineIndex, wASHI );

        double wASLO = w * higgsWeights.alphaS_dn / higgsWeights.nominal;
        if (!std::isfinite(wASLO)) wASLO = w;
        histoStore()->fillTH1F(  "h_catSTXS_alphaS_dn",   m_category, wASHI );
        histoStore()->fillTH2F( "h2_catSTXS_alphaS_dn",   m_category, STXSbin, wASHI );
        histoStore()->fillTH2F( "h2_fineIndex_alphaS_dn", m_category, fineIndex, wASHI );

        for (int ipdf(0); ipdf < 30; ipdf++)  {
          double wPDF = w * higgsWeights.pdf4lhc_unc[ipdf] / higgsWeights.nominal;
          if (!std::isfinite(wPDF)) wPDF = w;
          TString suffixPDF = TString::Format("_PDF%d",ipdf);
          histoStore()->fillTH1F(  "h_catSTXS"+suffixPDF,   m_category, wPDF );
          histoStore()->fillTH2F( "h2_catSTXS"+suffixPDF,   m_category, STXSbin, wPDF );
          histoStore()->fillTH2F( "h2_fineIndex"+suffixPDF, m_category, fineIndex, wPDF );
        }
      }

    }
    
  }

  return EL::StatusCode::SUCCESS;
}
