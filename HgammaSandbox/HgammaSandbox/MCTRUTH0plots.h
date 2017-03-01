#ifndef HgammaSandbox_MCTRUTH0plots_H
#define HgammaSandbox_MCTRUTH0plots_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include "HgammaSandbox/HGamTruthCategoryTool.h"

class MCTRUTH0plots : public HgammaAnalysis
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

  xAOD::TEvent *m_event; //!
  // Couplings categories
  #ifndef __CINT__
  HG::HGamTruthCategoryTool *m_catTool; //!
  #endif // __CINT__

public:
  // this is a standard constructor
  MCTRUTH0plots() { }
  MCTRUTH0plots(const char *name);
  virtual ~MCTRUTH0plots();

  // these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode execute();


  // this is needed to distribute the algorithm to the workers
  ClassDef(MCTRUTH0plots, 1);
};

#endif // HgammaSandbox_MCTRUTH0plots_H
