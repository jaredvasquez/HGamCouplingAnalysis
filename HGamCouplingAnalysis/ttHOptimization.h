#ifndef HGamCouplingAnalysis_ttHOptimization_H
#define HGamCouplingAnalysis_ttHOptimization_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"

class ttHOptimization : public HgammaAnalysis
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
private:
  // Tree *myTree; //!
  // TH1 *myHist; //!
  bool m_reweightHiggsPt;
  bool m_isGGH, m_isTWH;
  double ggF_01jet_hpT_weight(double);



public:
  // this is a standard constructor
  ttHOptimization() { }
  ttHOptimization(const char *name);
  virtual ~ttHOptimization();

  // these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode execute();


  // this is needed to distribute the algorithm to the workers
  ClassDef(ttHOptimization, 1);
};

#endif // HGamCouplingAnalysis_ttHOptimization_H
