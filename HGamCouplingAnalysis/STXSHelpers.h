#ifndef STXSHELPERS_H
#define STXSHELPERS_H

// other includes
#include <vector>

namespace STXS {
  
  // ********************
  //   Stage 0 Helpers 
  // ********************
  
  int stage0_to_index( int stage0 ) {
  // Helper function to define histogram binnings using stage-0 STXS code
  // ------------------------------------------------------------------------------------------
    int P = (int)(stage0 / 100);
    int F = (int)(stage0 % 100);
    static std::vector<int> offset({0,1,3,7,9,11,13,15,16,19});
    return ( F + offset[P] );
  }
  

  // ********************
  //   Stage 1 Helpers 
  // ********************
  

  int fixBins_stage1( int stage1, int NTruthJets=-1, double pTHiggs=-999. ) {
  // Helper function to correct issues in binning in versions before TruthRivetTools-00-00-07
  // ------------------------------------------------------------------------------------------
    if (stage1/100 >= 6) stage1--;
    else if ((pTHiggs > 200) && (stage1 == 101 || stage1 == 102)) {
      if (NTruthJets == 1) stage1 = 107;
      if (NTruthJets == 2) stage1 = 111;
    }

    return stage1;
  }


  int stage1_to_index( int stage1 ) {
  // ------------------------------------------------------------------------------------------
  // Helper function to define histogram binnings using stage-1 STXS code
    int P = (int)(stage1 / 100);
    int F = (int)(stage1 % 100);
    static std::vector<int> offset({0,1,13,19,24,29,33,35,37,39});
    if (stage1 < 0) return -1;
    return ( F + offset[P] );
  }


  int define_stage1HComb( int stage1, int prodMode ) { 
  // Helper function to redefine STXS codes needed for HComb
  // ------------------------------------------------------------------------------------------
    if (stage1/100 == 2)  return stage1 + (prodMode-2)*10;
    return stage1;
  }
  
  int stage1HComb_to_index( int stage1 ) {
  // Helper function to define histogram binnings using stage-1 HComb STXS code
  // ------------------------------------------------------------------------------------------
    int P = (int)(stage1 / 100);
    int G = (int)(stage1 / 10 );
    int F = (int)(stage1 % 100);
    static std::vector<int> offset({0,1,13,31,36,41,45,47,49,51});
    if (P == 2) return ( F + G*6 + offset[P] );
    return ( F + offset[P] );
  }
  

}

#endif // STXSHELPERS_H

