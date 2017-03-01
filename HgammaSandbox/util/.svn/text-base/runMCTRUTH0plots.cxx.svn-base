#include "HgammaSandbox/MCTRUTH0plots.h"
#include "HGamAnalysisFramework/RunUtils.h"

int main(int argc, char *argv[])
{
  // Set up the job for xAOD access
  xAOD::Init().ignore();

  // Create our algorithm
  MCTRUTH0plots *alg = new MCTRUTH0plots("MCTRUTH0plots");

  // Use helper to start the job
  HG::runJob(alg, argc, argv);

  return 0;
}
