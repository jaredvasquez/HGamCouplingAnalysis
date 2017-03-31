#ifndef HGamCouplingAnalysis_CouplingAnalysis_H
#define HGamCouplingAnalysis_CouplingAnalysis_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"

class CouplingAnalysis : public HgammaAnalysis
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)

private:
  bool m_useSystematics;
  bool m_reweightHiggsPt;
  bool m_usePDFUncerts;
  bool m_isGGH, m_isTWH;
  HG::SystematicList sysList;
  
  TTree *m_tree; //!
  int m_category; 
  double m_myy;
  
  TString qcdNames[9] = { "mu", "res", "mig01", "mig12", "vbf2j", "vbf3j", "pTH60", "pTH120", "qm_t" }; //!

public:
  // this is a standard constructor
  CouplingAnalysis() { }
  CouplingAnalysis(const char *name);
  virtual ~CouplingAnalysis();

  // these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode execute();


  // this is needed to distribute the algorithm to the workers
  ClassDef(CouplingAnalysis, 1);
};

#endif // HGamCouplingAnalysis_CouplingAnalysis_H
