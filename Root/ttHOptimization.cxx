#include "HGamCouplingAnalysis/ttHOptimization.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamAnalysisFramework/HGamVariables.h"

// Additional Helpers
//#include "HGamCouplingAnalysis/STXSHelpers.h"
//#include "HGamCouplingAnalysis/HiggsHelpers.h"

// this is needed to distribute the algorithm to the workers
ClassImp(ttHOptimization)



ttHOptimization::ttHOptimization(const char *name)
: HgammaAnalysis(name)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



ttHOptimization::~ttHOptimization()
{
  // Here you delete any memory you allocated during your analysis.
}



EL::StatusCode ttHOptimization::createOutput()
{
  // Here you setup the histograms needed for you analysis. This method
  // gets called after the Handlers are initialized, so that the systematic
  // registry is already filled.

  // Only apply the reweighting to ggH MC
  m_reweightHiggsPt = config()->getBool("HGamCoupling.ReweightHiggsPt", false);

  m_isGGH = false;
  m_isTWH = false;
  if (isMC()) {
    m_isGGH = getMCSampleName(eventInfo()->mcChannelNumber()).Contains("ggH125");
    m_isTWH = getMCSampleName(eventInfo()->mcChannelNumber()).Contains("tWH125");
    m_reweightHiggsPt &= m_isGGH;
  } else m_reweightHiggsPt = false;
    

  if (m_reweightHiggsPt) std::cout << "*** !!! REWEIGHTING HIGGS PT !!! ***" << std::endl;
  else               std::cout << "*** !!! NOT REWEIGHTING HIGGS PT !!! ***" << std::endl;
  
  // Create Histograms
  int nCats(10);
  histoStore()->createTH1F(  "h_catMxAOD",  nCats, 0.5, nCats+0.5 );
  histoStore()->createTH1F(  "h_catOpt",    nCats, 0.5, nCats+0.5 );


  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ttHOptimization::execute()
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
  if ( isMC() && m_reweightHiggsPt && (eventInfo()->auxdata<int>("HTXS_Njets_pTjet25") < 2) ) {
    w_pT = ggF_01jet_hpT_weight( var::pT_h1.truth()*HG::invGeV );
  }

  // get correction factor for N_init with pT reweighting
  double corrDenom = ( isMC() && m_reweightHiggsPt ) ? 1.0 : 1.003;
  
  //double wMC = (isData()) ? 1.0 : w_pT * eventHandler()->mcWeight();
  double w  = (isData()) ? 1.0 : w_pT * weightCatCoup_Moriond2017() * lumiXsecWeight();

  if (isMC() && not var::isPassed()) return EL::StatusCode::SUCCESS;
  else if ( isData() ) {
    bool failTI = ( var::cutFlow() == 11 || var::cutFlow() == 12 );
    if (!failTI) return EL::StatusCode::SUCCESS;

    // Retrieve photons
    xAOD::PhotonContainer all_photons  = photonHandler() -> getCorrectedContainer();
    xAOD::PhotonContainer loosePhotons = photonHandler() -> applyPreSelection(all_photons);

    // Relative pT cuts
    if ( !passRelativePtCuts(loosePhotons) ) return EL::StatusCode::SUCCESS;

    // Mass window cuts
    if ( !passMyyWindowCut(loosePhotons) ) return EL::StatusCode::SUCCESS;
  }

  w = (isData()) ? 1.0 : w_pT * corrDenom * weightCatCoup_Moriond2017() * lumiXsecWeight();
  if (w == 0.) return EL::StatusCode::SUCCESS;

  int catMxAOD = var::catCoup_Moriond2017() - 24;
  if (catMxAOD < 0) catMxAOD = 0;
  
  xAOD::ElectronContainer electrons = electronHandler()->getCorrectedContainer();
  xAOD::MuonContainer muons = muonHandler()->getCorrectedContainer();
  xAOD::JetContainer jets   = jetHandler()->getCorrectedContainer();
  
  // do my own categorization  
  int ttHCat = 0;
  double wCat = 1.0;
  int n_electrons = electrons.size();
  int n_muons = muons.size();
  int n_leptons = (n_electrons + n_muons);

  int n_jetsCen25(0), n_jetsCen35(0), n_jetsFwd25(0);
  int n_tags25_60(0), n_tags25_70(0), n_tags25_77(0);
  int n_tags35_60(0), n_tags35_70(0);

  double jetWeight25(1.0), jetWeight35(1.0); // tag and JVT SFs
  double bweight25_60(1.0), bweight25_70(1.0), bweight25_77(1.0);
  double bweight35_60(1.0), bweight35_70(1.0);

  // count only central jets s.t. |eta| < 2.5
  for ( auto jet: jets ) {
    if (fabs(jet->eta()) > 2.5) { n_jetsFwd25++; continue; }

    n_jetsCen25++;

    jetWeight25  *= jet->auxdata<float>("SF_jvt");
    bweight25_60 *= jet->auxdata<float>("SF_MV2c10_FixedCutBEff_60");
    bweight25_70 *= jet->auxdata<float>("SF_MV2c10_FixedCutBEff_70");
    bweight25_77 *= jet->auxdata<float>("SF_MV2c10_FixedCutBEff_77");

    if (jet->auxdata<char>("MV2c10_FixedCutBEff_60")) n_tags25_60++;
    if (jet->auxdata<char>("MV2c10_FixedCutBEff_70")) n_tags25_70++;
    if (jet->auxdata<char>("MV2c10_FixedCutBEff_77")) n_tags25_77++;

    // Count 35 GeV Jets and Tags
    if (jet->pt() > 35.e3) {
      n_jetsCen35++;

      jetWeight35  *= jet->auxdata<float>("SF_jvt");
      bweight35_60 *= jet->auxdata<float>("SF_MV2c10_FixedCutBEff_60");
      bweight35_70 *= jet->auxdata<float>("SF_MV2c10_FixedCutBEff_70");

      if (jet->auxdata<char>("MV2c10_FixedCutBEff_60")) n_tags35_60++;
      if (jet->auxdata<char>("MV2c10_FixedCutBEff_70")) n_tags35_70++;
    }

  }

  // force 70% WP everywhere
  n_tags25_60 = n_tags25_70;
  n_tags25_77 = n_tags25_70;
  n_tags35_60 = n_tags35_70;

  // Leptonic Channel Selection
  if (n_leptons >= 1) {

    wCat *= (jetWeight25 * bweight25_77); // jvt & b-tag SFs
    for( auto el : electrons ) wCat *= el->auxdata< float >("scaleFactor"); // lepton SFs
    for( auto mu : muons )     wCat *= mu->auxdata< float >("scaleFactor");

    bool passPreselTH = ( (n_leptons == 1) && (n_tags25_77 >= 1) );
    if      ( passPreselTH && (n_jetsCen25 <= 3) && (n_jetsFwd25 == 0) )  ttHCat = 9;
    else if ( passPreselTH && (n_jetsCen25 <= 4) && (n_jetsFwd25 >= 1) )  ttHCat = 8;
    else if (                 (n_jetsCen25 >= 2) && (n_tags25_77 >= 1) )  ttHCat = 7;

    double mll = -999.; // require OS dileptons in future?
    if ( n_muons >= 2 ) mll = ( muons[0]->p4() + muons[1]->p4() ).M() * HG::invGeV;
    if ( n_electrons >= 2 ) mll = ( electrons[0]->p4() + electrons[1]->p4() ).M() * HG::invGeV;
    if (fabs(mll-91) < 10) ttHCat = 0; // fail selection

  // Hadronic Channel Selection
  } else {

    if        ((n_jetsCen35 == 4) && (n_tags35_60 == 1)) {
      wCat *= (jetWeight35 * bweight35_60); ttHCat = 6;

    } else if ((n_jetsCen35 == 4) && (n_tags35_70 >= 2)) {
      wCat *= (jetWeight35 * bweight35_70); ttHCat = 5;

    } else if ((n_jetsCen25 == 5) && (n_tags25_60 == 1)) {
      wCat *= (jetWeight25 * bweight25_60); ttHCat = 4;

    } else if ((n_jetsCen25 == 5) && (n_tags25_70 >= 2)) {
      wCat *= (jetWeight25 * bweight25_70); ttHCat = 3;

    } else if ((n_jetsCen25 >= 6) && (n_tags25_60 == 1)) {
      wCat *= (jetWeight25 * bweight25_60); ttHCat = 2;

    } else if ((n_jetsCen25 >= 6) && (n_tags25_70 >= 2)) {
      wCat *= (jetWeight25 * bweight25_70); ttHCat = 1;
    }
  }

  //wCat = (isData()) ? 1.0 : w_pT * corrDenom * wCat * lumiXsecWeight();
  wCat = (isData()) ? 1.0 : w_pT * corrDenom * weight() * lumiXsecWeight();

  histoStore()->fillTH1F(  "h_catMxAOD",  catMxAOD, wCat );
  histoStore()->fillTH1F(  "h_catOpt",    ttHCat,   wCat );
  
  // do weights and get njets and ntags needed.
  return EL::StatusCode::SUCCESS;
}

double ttHOptimization::ggF_01jet_hpT_weight(double H_pT) {
  if ( H_pT <  20 )  return 1.11;
  if ( H_pT <  45 )  return 1.11 - (H_pT-20)/25*0.2;  // -> 0.91
  if ( H_pT < 135 )  return 0.91 - (H_pT-45)/90*0.36; // -> 0.55
  return 0.55;
}

