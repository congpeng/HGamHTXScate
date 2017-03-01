#ifndef HGamTruthCategoryTool_H
#define HGamTruthCategoryTool_H

#include <EventLoop/Worker.h>
#include "HGamAnalysisFramework/EventHandler.h"
#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"

namespace TMVA {
  class Reader;
}

namespace HG {
  class HGamTruthCategoryTool {

  private:
    xAOD::TEvent *m_event; //!
    xAOD::TStore *m_store; //!

    TMVA::Reader *m_readerVBF_low; //!
    TMVA::Reader *m_readerVBF_high; //!
    TMVA::Reader *m_readerVH_had_Moriond2017; //!

    int t_Njets30;
    float t_pTt_yyGeV;
    float t_m_jjGeV;
    float t_pTt_yy;
    float t_m_jj;
    float t_dEta_jj;
    float t_dPhi_yy_jj;
    float t_Zepp;
    float t_Dy_yy_jj;
    float t_Drmin_y_j2;
    float t_cosTS_yy_jj;
    float t_cosTS_yy;

    double m_weight;

    float getMVAWeight(TMVA::Reader *XReader);
    
    void prepVars( xAOD::TruthParticleContainer &gams, xAOD::TruthParticleContainer &electrons, xAOD::TruthParticleContainer &muons, xAOD::JetContainer &jets, xAOD::JetContainer &bjets, xAOD::MissingETContainer &met); 
    
    int ttH_split( xAOD::TruthParticleContainer &gams, xAOD::TruthParticleContainer &electrons, xAOD::TruthParticleContainer &muons, xAOD::JetContainer &jets, xAOD::JetContainer &bjets, xAOD::MissingETContainer &met); 
    int VH2lep_split( xAOD::TruthParticleContainer &gams, xAOD::TruthParticleContainer &electrons, xAOD::TruthParticleContainer &muons, xAOD::JetContainer &jets, xAOD::JetContainer &bjets, xAOD::MissingETContainer &met);
    int VHlep_split( xAOD::TruthParticleContainer &gams, xAOD::TruthParticleContainer &electrons, xAOD::TruthParticleContainer &muons, xAOD::JetContainer &jets, xAOD::JetContainer &bjets, xAOD::MissingETContainer &met);
    int VHMET_split( xAOD::TruthParticleContainer &gams, xAOD::TruthParticleContainer &electrons, xAOD::TruthParticleContainer &muons, xAOD::JetContainer &jets, xAOD::JetContainer &bjets, xAOD::MissingETContainer &met);
    bool Passes_qqH_BSM(xAOD::JetContainer &jets);
    bool pass_VHhad( xAOD::JetContainer &jets );
    int Passes_VBF( xAOD::TruthParticleContainer &gams, xAOD::JetContainer &jets);
    int ggH_split( xAOD::TruthParticleContainer &selphotons, xAOD::JetContainer &jets);

  public:
    HGamTruthCategoryTool(xAOD::TEvent *event, xAOD::TStore *store);
    virtual ~HGamTruthCategoryTool();

    virtual EL::StatusCode initialize(Config &config);

    // Use below for splitting ggH categories
    enum HGamM17Index { 
      M17_FailDiphoton=-1, M17_ggH_0J_Cen=1, M17_ggH_0J_Fwd=2, M17_ggH_1J_LOW=3, M17_ggH_1J_MED=4, M17_ggH_1J_HIGH=5, 
      M17_ggH_1J_BSM=6, M17_ggH_2J_LOW=7, M17_ggH_2J_MED=8, M17_ggH_2J_HIGH=9, M17_ggH_2J_BSM=10, M17_VBF_HjjLOW_loose=11, 
      M17_VBF_HjjLOW_tight=12, M17_VBF_HjjHIGH_loose=13, M17_VBF_HjjHIGH_tight=14, M17_VHhad_loose=15, M17_VHhad_tight=16, 
      M17_qqH_BSM=17, M17_VHMET_LOW=18, M17_VHMET_HIGH=19, M17_VHMET_BSM=20, M17_VHlep_LOW=21, M17_VHlep_HIGH=22, 
      M17_VHdilep_LOW=23, M17_VHdilep_HIGH=24, M17_ttH=25 
    };

    int getEventCategory( xAOD::TruthParticleContainer &gams, xAOD::TruthParticleContainer &electrons, xAOD::TruthParticleContainer &muons, xAOD::JetContainer &jets, xAOD::JetContainer &bjets, xAOD::MissingETContainer &met); 

  };

}
#endif // HGamTruthCategoryTool_H
