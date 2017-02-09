#include "HGamCouplingAnalysis/ttHOptimization.h"
#include "HGamAnalysisFramework/RunUtils.h"

int main(int argc, char *argv[])
{
  // Set up the job for xAOD access
  xAOD::Init().ignore();

  // Create our algorithm
  ttHOptimization *alg = new ttHOptimization("ttHOptimization");

  // Use helper to start the job
  HG::runJob(alg, argc, argv);

  return 0;
}
