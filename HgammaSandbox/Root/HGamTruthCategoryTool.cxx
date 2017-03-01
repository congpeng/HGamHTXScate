#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

#include <HGamAnalysisFramework/HgammaAnalysis.h>
#include <HGamAnalysisFramework/HgammaUtils.h>
#include <HgammaSandbox/HGamTruthCategoryTool.h>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


namespace HG {

  //______________________________________________________________________________
  HGamTruthCategoryTool::HGamTruthCategoryTool(xAOD::TEvent *event, xAOD::TStore *store)
  : m_event(event)
  , m_store(store)
  , m_readerVBF_low(nullptr)
  , m_readerVBF_high(nullptr)
  , m_readerVH_had_Moriond2017(nullptr)
  { }


  //______________________________________________________________________________
  HGamTruthCategoryTool::~HGamTruthCategoryTool()
  {
    SafeDelete(m_readerVBF_high);
    SafeDelete(m_readerVBF_low);
    SafeDelete(m_readerVH_had_Moriond2017);
  }


  //______________________________________________________________________________
  EL::StatusCode HGamTruthCategoryTool::initialize(Config &config)
  {
    // Setup MVA reader

     // VBF_high
    m_readerVBF_high = new TMVA::Reader( "!Color:!Silent" );
    m_readerVBF_high->AddVariable("pTt_yy",             &t_pTt_yy);
    m_readerVBF_high->AddVariable("m_jj",               &t_m_jj);
    m_readerVBF_high->AddVariable("jj_DeltaEta",        &t_dEta_jj);
    m_readerVBF_high->AddVariable("Dphi_yy_jj",     &t_dPhi_yy_jj);
    m_readerVBF_high->AddVariable("abs(Zepp)",&t_Zepp);
    m_readerVBF_high->AddVariable("DRmin_y_j",     &t_Drmin_y_j2);
    TString readerPathVBF_high = PathResolverFindCalibFile("HGamAnalysisFramework/TMVA_VBF_high.weights.xml");
    m_readerVBF_high->BookMVA( "BDTG",readerPathVBF_high);

    //  VBF low
    m_readerVBF_low = new TMVA::Reader( "!Color:!Silent" );
    m_readerVBF_low->AddVariable("pTt_yy",             &t_pTt_yy);
    m_readerVBF_low->AddVariable("m_jj",               &t_m_jj);
    m_readerVBF_low->AddVariable("jj_DeltaEta",        &t_dEta_jj);
    m_readerVBF_low->AddVariable("Dphi_yy_jj",     &t_dPhi_yy_jj);
    m_readerVBF_low->AddVariable("abs(Zepp)",&t_Zepp);
    m_readerVBF_low->AddVariable("DRmin_y_j",     &t_Drmin_y_j2);
    TString readerPathVBF_low = PathResolverFindCalibFile("HGamAnalysisFramework/TMVA_VBF_low.weights.xml");
    m_readerVBF_low->BookMVA( "BDTG",readerPathVBF_low);


    // Setup VH_had MVA reader
    m_readerVH_had_Moriond2017 = new TMVA::Reader( "!Color:!Silent" );
    m_readerVH_had_Moriond2017->AddVariable("m_jj",               &t_m_jjGeV);
    m_readerVH_had_Moriond2017->AddVariable("pTt_yy",             &t_pTt_yyGeV);
    m_readerVH_had_Moriond2017->AddVariable("Dy_yy_jj",           &t_Dy_yy_jj);
    m_readerVH_had_Moriond2017->AddVariable("cosTS_yyjj",        &t_cosTS_yy_jj);
    TString readerPathVH  = PathResolverFindCalibFile("HGamAnalysisFramework/MVA_config_VH_had.xml");
    m_readerVH_had_Moriond2017->BookMVA( "BDTG", readerPathVH);

    return EL::StatusCode::SUCCESS;
  }
  
  
  //______________________________________________________________________________


  int HGamTruthCategoryTool::ttH_split(xAOD::TruthParticleContainer &gams, xAOD::TruthParticleContainer &electrons, xAOD::TruthParticleContainer &muons, xAOD::JetContainer &jets, xAOD::JetContainer &bjets, xAOD::MissingETContainer &met) 
  {
    int n_electrons = electrons.size();
    int n_muons = muons.size();
    int n_leptons = (n_electrons + n_muons);

    int n_jetsCen25(0), n_jetsCen35(0), n_jetsFwd25(0);
    int n_tags25_70(0), n_tags35_70(0);

    double bweight25_70(1.0), bweight35_70(1.0);

    // count only central jets s.t. |eta| < 2.5
    for ( auto jet:jets ) {
      if (fabs(jet->eta()) > 2.5) { n_jetsFwd25++; continue; }

      // Count 25 GeV Jets and Tags
      n_jetsCen25++;
      // Count 35 GeV Jets and Tags
      if (jet->pt() >= 35.e3) {
        n_jetsCen35++;
      }
    }

    for ( auto jet:bjets) {
      if (fabs(jet->eta()) > 2.5) continue;

      // Count 25 GeV b Jets and Tags
      n_tags25_70++;
      // Count 35 GeV b Jets and Tags
      if (jet->pt() >= 35.e3) {
	n_tags35_70++;
      }
    }
    
    // Leptonic Channel Selection
    if (n_leptons >= 1) {

      double mll = -999.; // require OS dileptons in future?
      if ( n_muons >= 2 ) mll = ( muons[0]->p4() + muons[1]->p4() ).M() * HG::invGeV;
      if ( n_electrons >= 2 ) mll = ( electrons[0]->p4() + electrons[1]->p4() ).M() * HG::invGeV;
      if (fabs(mll-91) < 10) return 0; // fail selection

      bool passPreselTH = ( (n_leptons == 1) && (n_tags25_70 >= 1) );
      if      ( passPreselTH && (n_jetsCen25 <= 3) && (n_jetsFwd25 == 0) )  return 9;
      else if ( passPreselTH && (n_jetsCen25 <= 4) && (n_jetsFwd25 >= 1) )  return 8;
      else if (                 (n_jetsCen25 >= 2) && (n_tags25_70 >= 1) )  return 7;

    // Hadronic Channel Selection
    } else {

      if        ((n_jetsCen35 == 4) && (n_tags35_70 == 1)) {
        return 6;

      } else if ((n_jetsCen35 == 4) && (n_tags35_70 >= 2)) {
        return 5;

      } else if ((n_jetsCen25 == 5) && (n_tags25_70 == 1)) {
        return 4;

      } else if ((n_jetsCen25 == 5) && (n_tags25_70 >= 2)) {
        return 3;

      } else if ((n_jetsCen25 >= 6) && (n_tags25_70 == 1)) {
        return 2;

      } else if ((n_jetsCen25 >= 6) && (n_tags25_70 >= 2)) {
        return 1;
      }
    }

    // should only reach here if failed selection
    return 0;
  }

  int HGamTruthCategoryTool::VH2lep_split(xAOD::TruthParticleContainer &gams, xAOD::TruthParticleContainer &electrons, xAOD::TruthParticleContainer &muons, xAOD::JetContainer &jets, xAOD::JetContainer &bjets, xAOD::MissingETContainer &met) 
   {
      int m_Nelectrons = electrons.size();
      int m_Nmuons = muons.size();
      float m_m_mumu=0, m_m_ee=0, pT_m_ee=0,pT_m_mumu=0;
      if (m_Nelectrons>1) {
          m_m_ee =(electrons[0]->p4() + electrons[1]->p4()).M();
          pT_m_ee = (electrons[0]->p4() + electrons[1]->p4()).Pt();
      }
      if (m_Nmuons>1){
          m_m_mumu =( muons[0]->p4() + muons[0]->p4()).M();
          pT_m_mumu =( muons[0]->p4() + muons[0]->p4()).Pt();
      }

      bool passMuSel = ( m_Nmuons>=2 && m_m_mumu > 70*HG::GeV && m_m_mumu < 110*HG::GeV );
      bool passElSel = ( m_Nelectrons>=2 && m_m_ee > 70*HG::GeV && m_m_ee < 110*HG::GeV );
      bool pTVlt150 = ((passMuSel&&pT_m_mumu<150*HG::GeV)||(passElSel&&pT_m_ee<150*HG::GeV));
      bool pTVgt150 = ((passMuSel&&pT_m_mumu>=150*HG::GeV)||(passElSel&&pT_m_ee>=150*HG::GeV));


      if (pTVlt150)return 1;
      else if (pTVgt150)return 2;
      else return 0;
   }

  int HGamTruthCategoryTool::VHlep_split(xAOD::TruthParticleContainer &gams, xAOD::TruthParticleContainer &electrons, xAOD::TruthParticleContainer &muons, xAOD::JetContainer &jets, xAOD::JetContainer &bjets, xAOD::MissingETContainer &met) 
   {
      int m_Nelectrons = electrons.size();
      int m_Nmuons = muons.size();
      if ((m_Nelectrons+m_Nmuons)!=1) return 0;
      bool passZey(true);
      for ( auto el: electrons ) {
        double mey1 = (el->p4() + gams[0]->p4()).M() * HG::invGeV;
        double mey2 = (el->p4() + gams[1]->p4()).M() * HG::invGeV;
        if ( (fabs(mey1-89) < 5) || (fabs(mey2-89) < 5) ) passZey = false;
      }
    
      if(!passZey) return 0;
      float TST_met = -99, TST_sumet = -99;
      if (met["TST"] != nullptr) {
	//get truth MET:
        TST_met   = met["TST"]->met();
        TST_sumet = met["TST"]->sumet();
      }
      float m_MET_signi=-99999;
      if (TST_met!=0) m_MET_signi=TST_met*HG::invGeV/sqrt(TST_sumet*HG::invGeV);
    
      TLorentzVector metVec;
      metVec.SetPtEtaPhiM(met["TST"]->met(), 0, met["TST"]->phi(), 0);
    
      bool pTVlt150 = true;
      if(((m_Nelectrons==1) && ((electrons[0]->p4() + metVec).Pt() > 150*HG::GeV)) ||
         ((m_Nmuons==1    ) && ((muons[0]->p4() + metVec).Pt() > 150*HG::GeV)))
        pTVlt150 = false;
    
      if(pTVlt150&&m_MET_signi>1)return 1;
      else if(!pTVlt150)return 2;
      else return 0;
   }

  int HGamTruthCategoryTool::VHMET_split(xAOD::TruthParticleContainer &gams, xAOD::TruthParticleContainer &electrons, xAOD::TruthParticleContainer &muons, xAOD::JetContainer &jets, xAOD::JetContainer &bjets, xAOD::MissingETContainer &met) 
   {
    float TST_met = -99, TST_sumet = -99;
    if (met["TST"] != nullptr) {
      TST_met   = met["TST"]->met();
      TST_sumet = met["TST"]->sumet();
    }
    float m_MET_signi=-99999;
    if (TST_met!=0) m_MET_signi=TST_met*0.001/sqrt(TST_sumet*0.001);

    if(TST_met<150*HG::GeV &&TST_met>=80*HG::GeV && m_MET_signi>8) return 1;
    else if(TST_met<250*HG::GeV &&TST_met>=150*HG::GeV && m_MET_signi>9) return 2;
    else if(TST_met>=250*HG::GeV) return 3;
    else return 0;
   }

  bool HGamTruthCategoryTool::Passes_qqH_BSM(xAOD::JetContainer &jets)
   {
    if(jets[0]->pt()>200*HG::GeV) return true;
    return false;
   }

  bool HGamTruthCategoryTool::pass_VHhad( xAOD::JetContainer &jets ) {
    if (jets.size() < 2) return 0;
    float m_jj = (jets[0]->p4() + jets[1]->p4()).M();
    if (m_jj <  60*HG::GeV) return 0;
    if (m_jj > 120*HG::GeV) return 0;
    return true;
  }

  int HGamTruthCategoryTool::Passes_VBF( xAOD::TruthParticleContainer &gams, xAOD::JetContainer &jets)
   {
    if(t_Njets30<2) return 0;
    if(t_dEta_jj<2) return 0;
    if(t_Zepp>5) return 0;

    float pt_yyjj= (gams[0]->p4()+gams[1]->p4()+jets[0]->p4()+jets[1]->p4()).Pt();
    if ( pt_yyjj > 25*HG::GeV) return 2;
    else return 1;
   }

  //______________________________________________________________________________
  int HGamTruthCategoryTool::ggH_split(xAOD::TruthParticleContainer &selphotons, xAOD::JetContainer &jets) {

    int cat=-1;
    if(t_Njets30==0){
      if(fabs(selphotons[0]->eta()) <= 0.95 && fabs(selphotons[1]->eta()) <= 0.95) cat=1;
      else cat = 2;
    }	
    else if(t_Njets30==1){
      double yypt = (selphotons[0]->p4()+selphotons[1]->p4()).Pt();
      if      ( yypt <= 60*HG::GeV)                          cat=3;
      else if ( yypt >  60*HG::GeV && yypt<=120*HG::GeV)     cat=4;
      else if ( yypt > 120*HG::GeV && yypt<=200*HG::GeV)     cat=5;
      else if ( yypt > 200*HG::GeV)                          cat=6;
    }
    else if(t_Njets30>1){
      double yypt = (selphotons[0]->p4()+selphotons[1]->p4()).Pt();
      if      ( yypt <= 60*HG::GeV )                         cat=7;
      else if ( yypt >  60*HG::GeV && yypt <= 120*HG::GeV )  cat=8;
      else if ( yypt > 120*HG::GeV && yypt <= 200*HG::GeV )  cat=9;
      else if ( yypt > 200*HG::GeV )                         cat=10;
    }
    std::cerr << "uh-oh! The categorization tool is being naughty..." << std::endl;
    return cat; // should never reach this point
  }
  

  //______________________________________________________________________________
  void HGamTruthCategoryTool::prepVars(xAOD::TruthParticleContainer &gams, xAOD::TruthParticleContainer &electrons, xAOD::TruthParticleContainer &muons, xAOD::JetContainer &jets, xAOD::JetContainer &bjets, xAOD::MissingETContainer &met) 
   {
    t_Njets30=-999;
    t_pTt_yy      = -999;
    t_m_jj        = -999;
    t_pTt_yyGeV      = -999;
    t_m_jjGeV        = -999;
    t_cosTS_yy_jj = -999;
    t_dEta_jj     = -999;
    t_dPhi_yy_jj  = -999;
    t_Zepp        = -999;
    t_Drmin_y_j2  = -999;
    t_Dy_yy_jj    = -999;

    if (gams.size() < 2) return; 

    for(auto jet:jets)
     {
      if(jet->pt()>=30.e3) t_Njets30++;
     }

    //VH2jet:
    TLorentzVector g1 = gams[0]->p4(), g2 = gams[1]->p4();
    t_pTt_yyGeV = (fabs(g1.Px()*g2.Py() - g2.Px()*g1.Py())/(g1-g2).Pt()*2.0)/1000.;
    if(t_Njets30>1)
     {
      TLorentzVector j1 = jets[0]->p4(),j2 = jets[1]->p4();
      t_m_jjGeV= (j1+j2).M()/1000.;
      t_Dy_yy_jj = fabs((gams[0]->p4() + gams[1]->p4()).Rapidity() - (jets[0]->p4() + jets[1]->p4()).Rapidity());
      //define cosTS_yyjj below:
       {
	TLorentzVector vH, vZ, vg, vl1, vl2;
	vg = jets[0]->p4()+jets[1]->p4();
	vl1 = gams[0]->p4();
	vl2 = gams[1]->p4(); 
	vZ = gams[0]->p4() + gams[1]->p4();
	vH = vZ+vg;

	TVector3 boost = -vH.BoostVector();
	vH.Boost(boost); 
	vZ.Boost(boost);
	vg.Boost(boost);
	vl1.Boost(boost);
	vl2.Boost(boost);

	TLorentzVector q, qbar;
	q.SetPxPyPzE(0, 0, vH.M()/2, vH.M()/2);
	qbar.SetPxPyPzE(0, 0, -vH.M()/2, vH.M()/2);

	t_cosTS_yy_jj = (q-qbar).Dot(vZ)/(vH.M() * vZ.P());
       }
     }


    //VBF cate:
      t_pTt_yy = fabs(g1.Px()*g2.Py() - g2.Px()*g1.Py())/(g1-g2).Pt()*2.0;
    if(t_Njets30>1)
     {
      t_m_jj= (jets[0]->p4()+jets[1]->p4()).M();
      t_dEta_jj = fabs(jets[0]->eta() - jets[1]->eta());
      t_dPhi_yy_jj = fabs((gams[0]->p4() + gams[1]->p4()).DeltaPhi(jets[0]->p4() + jets[1]->p4()));
      t_Zepp = (gams[0]->p4() + gams[1]->p4()).Eta() - (jets[0]->eta() + jets[1]->eta())/2.0;
      //define drmin_yj2:
       {
	double dR2min = 99.0, dR2 = 0.0, eta = 0.0, phi = 0.0;
	for (int i = 0; i < 2; i++) {
	  auto gam = gams[i];
	  eta = gam->eta(); phi = gam->phi();
	  for (int j=0; j<2;j++) {
	    auto jet = jets[j] ;
	    dR2 = xAOD::P4Helpers::deltaR2(*jet, eta, phi, false);
	    if (dR2 < dR2min) dR2min = dR2;
	  }
	}
	if (dR2min == 99)  t_Drmin_y_j2=-99;
	else t_Drmin_y_j2 = sqrt(dR2min);
       }
     }
   }

  //______________________________________________________________________________
  int HGamTruthCategoryTool::getEventCategory(xAOD::TruthParticleContainer &gams, xAOD::TruthParticleContainer &electrons, xAOD::TruthParticleContainer &muons, xAOD::JetContainer &jets, xAOD::JetContainer &bjets, xAOD::MissingETContainer &met) 
   { 
    int catIndex=-999;
    //    ttH Categories
    // --------------------------------
    int CtthCat = ttH_split(gams,electrons,muons,jets,bjets,met);
    if(CtthCat>0) catIndex=M17_ttH+CtthCat-1;

    //   qq->VH Categories  3+2+2
    // --------------------------------
    int CVH2l  = VH2lep_split(gams,electrons,muons,jets,bjets,met);
    if (CVH2l>=2)                catIndex=M17_VHdilep_HIGH;
    if (CVH2l==1)                catIndex=M17_VHdilep_LOW;

    int CVHlep = VHlep_split(gams,electrons,muons,jets,bjets,met);
    if (CVHlep>=2)               catIndex = M17_VHlep_HIGH;
    if (CVHlep==1)               catIndex= M17_VHlep_LOW;

    int CVHMET = VHMET_split(gams,electrons,muons,jets,bjets,met);
    if (CVHMET>=3)               catIndex= M17_VHMET_BSM;    
    if (CVHMET==2)               catIndex= M17_VHMET_HIGH;   
    if (CVHMET==1)               catIndex= M17_VHMET_LOW;    

    //   qq->H 1 category
    // --------------------------------
    if (Passes_qqH_BSM(jets)) catIndex=M17_qqH_BSM;

    //   VH hadronic 2 categories
    // --------------------------------
    if (pass_VHhad(jets))
     {
      float wMVA_VH = getMVAWeight(m_readerVH_had_Moriond2017);
      if ( wMVA_VH >  0.73 )    return catIndex = M17_VHhad_tight;
      if ( wMVA_VH >  0.18 )    return catIndex = M17_VHhad_loose;
     }

    //   VBF categories 2+2
    // --------------------------------
    int CVBF = Passes_VBF( gams, jets);
    if(CVBF==2)
     {
      float wMVA_VBF = getMVAWeight(m_readerVBF_high);
      if ( wMVA_VBF >  0.48 )    catIndex= M17_VBF_HjjHIGH_tight;
      if ( wMVA_VBF > -0.50 )    catIndex= M17_VBF_HjjHIGH_loose;
     } 
    if(CVBF==1)
     {
      float wMVA_VBF = getMVAWeight(m_readerVBF_low);
      if ( wMVA_VBF >  0.87 )    catIndex= M17_VBF_HjjLOW_tight;
      if ( wMVA_VBF > -0.39 )    catIndex= M17_VBF_HjjLOW_loose;
     }

    //  ggH categories 10
    //  -------------------------------
    catIndex = ggH_split( gams, jets );

    return catIndex;
  }

}
