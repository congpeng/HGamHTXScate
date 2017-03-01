#include "HgammaSandbox/MCTRUTH0plots.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamAnalysisFramework/TruthUtils.h"
#include <algorithm>

// this is needed to distribute the algorithm to the workers
ClassImp(MCTRUTH0plots)

bool ComparePt (xAOD::TruthParticle_v1 *i,xAOD::TruthParticle_v1 *j) { return (i->pt()>j->pt()); }


MCTRUTH0plots::MCTRUTH0plots(const char *name)
: HgammaAnalysis(name)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



MCTRUTH0plots::~MCTRUTH0plots()
{
  // Here you delete any memory you allocated during your analysis.
}



EL::StatusCode MCTRUTH0plots::createOutput()
{
  // Here you setup the histograms needed for you analysis. This method
  // gets called after the Handlers are initialized, so that the systematic
  // registry is already filled.
  
  // Initialize HGamCoupling Categorization Tool
  m_catTool = new HG::HGamTruthCategoryTool( event(), store() );
  m_catTool->initialize( *config() );

  //m_weightIndex_nom    = weightTool.spawnIndexRetriever(" nnlops-nominal-pdflhc ")->getIndex();

  //TRUTH HISTOGRAMS
  histoStore()->createTH2F("newstable_ptH",40,0,200,200,0.5,1199.5,";truth #it{p}_{THiggs} [GeV];#it{N}_{stable particles}");
  histoStore()->createTH2F("newpthads_ptH",40,0,200,60,0,400,";truth #it{p}_{THiggs} [GeV];truth #it{H}_{Thadrons} [GeV];");
  histoStore()->createTH1F("newpt_hads", 60, 0, 400,";truth #it{H}_{Thadrons} [GeV]");
  histoStore()->createTH1F("pt_yy", 40, 0, 200,";truth #it{p}_{T#gamma#gamma} [GeV]");
  histoStore()->createTH1F("newpt_yy", 40, 0, 200,";truth #it{p}_{T#gamma#gamma} [GeV]");
  histoStore()->createTH1F("newpt_H", 40, 0, 200,";truth #it{p}_{THiggs} [GeV]");
  histoStore()->createTH1F("newpt_jets", 40, 0, 200,";truth #it{p}_{Tjet} [GeV]");
  histoStore()->createTH1F("m_yy", 50, 110, 140,";truth #it{m}_{#gamma#gamma} [GeV]");
  histoStore()->createTH1F("Rapidity_yy", 32, -8, 8,";truth #it{y}_{#gamma#gamma} ");
  histoStore()->createTH1F("prapidity_yy", 32, -8, 8,";truth #it{#eta}_{#gamma#gamma} ");
  histoStore()->createTH1F("eta_truth_photon",12,-3,3,";truth #eta_{#gamma}");
  histoStore()->createTH1F("pt_truth_photon",40,0,200,";truth #it{p}_{T#gamma} [GeV]");
  histoStore()->createTH1F("quark_n",10,-0.5,9.5,";#it{N}_{truth-quarks}");
  histoStore()->createTH1F("quark_n25",10,-0.5,9.5,";#it{N}_{truth-quarks}");
  histoStore()->createTH1F("bquark_n",10,-0.5,9.5,";#it{N}_{truth-bquarks}");
  histoStore()->createTH1F("bquark_n25",10,-0.5,9.5,";#it{N}_{truth-bquarks}");
  histoStore()->createTH1F("Bhad_n",10,-0.5,9.5,";#it{N}_{truth-Bhadrons}");
  histoStore()->createTH1F("Bhad_n5",10,-0.5,9.5,";#it{N}_{truth-Bhadrons}");
  histoStore()->createTH1F("Dhad_n",10,-0.5,9.5,";#it{N}_{truth-Dhadrons}");
  histoStore()->createTH1F("Dhad_n25",10,-0.5,9.5,";#it{N}_{truth-Dhadrons}");
  histoStore()->createTH1F("ph_n",10,-0.5,9.5,";#it{N}_{truth-photons}");
  histoStore()->createTH1F("newph_n",10,-0.5,9.5,";#it{N}_{truth-photons}");
  histoStore()->createTH1F("newH_n",10,-0.5,9.5,";#it{N}_{Higgs}");
  histoStore()->createTH1F("had_n",40,-0.5,39.5,";#it{N}_{Hadrons}");
  histoStore()->createTH1F("newstable_n",200,0.5,1199.5,";#it{N}_{stable particles}");
  histoStore()->createTH1F("newjets_n",10,-0.5,9.5,";#it{N}_{jets}");
  histoStore()->createTH1F("ph_n25",10,-0.5,9.5,";#it{N}_{truth-photons}");
  histoStore()->createTH1F("elec_n",10,-0.5,9.5,";#it{N}_{truth-electrons}");
  histoStore()->createTH1F("elec_n25",10,-0.5,9.5,";#it{N}_{truth-electrons}");
  histoStore()->createTH1F("muon_n",10,-0.5,9.5,";#it{N}_{truth-muons}");
  histoStore()->createTH1F("muon_n25",10,-0.5,9.5,";#it{N}_{truth-muons}");
  histoStore()->createTH1F("phHiggs_n",10,-0.5,9.5,";#it{N}_{#gamma-Higgs}");
  histoStore()->createTH1F("phHiggs_n25",10,-0.5,9.5,";#it{N}_{#gamma-Higgs}");
  histoStore()->createTH1F("muonsBs_n",10,-0.5,9.5,";#it{N}_{muonsFromBs}");
  histoStore()->createTH1F("muonsBs_n25",10,-0.5,9.5,";#it{N}_{muonsFromBs}");
  histoStore()->createTH2F("m_yy_vs_m_qq", 100, 0, 1000,50,110,140,";truth #it{m}_{qq} [GeV];truth #it{m}_{#gamma#gamma} [GeV]");

  histoStore()->createTH1F("prapidity_qq", 16, -8, 8,";truth #it{#eta}_{qq} ");
  histoStore()->createTH1F("rapidity_qq", 16, -8, 8,";truth #it{y}_{qq} ");
  histoStore()->createTH1F("rapidity_truth_q", 16, -8, 8,";truth #it{y}_{quarks} ");
  histoStore()->createTH1F("rapidity_truth_q_lead", 16, -8, 8,";truth #it{y}_{quarks} ");
  histoStore()->createTH1F("deltaphi_truth_q", 10, -0.5, 4.5,";truth #it{delta#phi}_{quarks} ");
  histoStore()->createTH1F("deltar_truth_q", 10, -0.5, 9.5,";truth #it{deltaR}_{quarks} ");
  histoStore()->createTH1F("drapidity_truth_q", 16, 0, 8,";truth #it{deltay}_{quarks} ");
  histoStore()->createTH1F("pt_qq", 40, 0, 200,";truth #it{p}_{Tqq} [GeV]");
  histoStore()->createTH1F("m_qq", 40, 0, 400,";truth #it{m}_{qq} [GeV]");
  histoStore()->createTH1F("pt_truth_q_lead", 40, 0, 200,";truth #it{p}_{Tquarks} [GeV]");
  histoStore()->createTH1F("prapidity_truth_q_lead", 16, -8, 8,";truth #it{#eta}_{quarks} ");
  histoStore()->createTH1F("pt_truth_q", 40, 0, 200,";truth #it{p}_{Tquarks} [GeV]");
  histoStore()->createTH1F("prapidity_truth_q", 16, -8, 8,";truth #it{#eta}_{quarks} ");
  histoStore()->createTH2F("m_yy_vs_m_bb", 100, 0, 1000,50,110,140,";truth #it{m}_{bb} [GeV];truth #it{m}_{#gamma#gamma} [GeV]");
  histoStore()->createTH1F("rapidity_truth_b", 16, -8, 8,";truth #it{y}_{bquarks} ");
  histoStore()->createTH1F("rapidity_truth_b_lead", 16, -8, 8,";truth #it{y}_{bquarks} ");
  histoStore()->createTH1F("deltaphi_truth_b", 10, -0.5, 4.5,";truth #it{delta#phi}_{bquarks} ");
  histoStore()->createTH1F("deltar_truth_b", 10, -0.5, 9.5,";truth #it{deltaR}_{bquarks} ");
  histoStore()->createTH1F("drapidity_truth_b", 16, 0, 8,";truth #it{deltay}_{bquarks} ");
  histoStore()->createTH1F("pt_bb", 40, 0, 200,";truth #it{p}_{Tbb} [GeV]");
  histoStore()->createTH1F("m_bb", 40, 0, 400,";truth #it{m}_{bb} [GeV]");
  histoStore()->createTH1F("prapidity_bb", 16, -8, 8,";truth #it{#eta}_{bb} ");
  histoStore()->createTH1F("rapidity_bb", 16, -8, 8,";truth #it{y}_{bb} ");
  histoStore()->createTH1F("pt_truth_b_lead", 40, 0, 200,";truth #it{p}_{Tbquarks} [GeV]");
  histoStore()->createTH1F("prapidity_truth_b_lead", 16, -8, 8,";truth #it{#eta}_{bquarks} ");
  histoStore()->createTH1F("pt_truth_b", 40, 0, 200,";truth #it{p}_{Tbquarks} [GeV]");
  histoStore()->createTH1F("prapidity_truth_b", 16, -8, 8,";truth #it{#eta}_{bquarks} ");
 
  histoStore()->createTH1F("rapidity_truth_y_lead", 16, -8, 8,";truth #it{y}_{photon} ");
  histoStore()->createTH1F("prapidity_truth_y_lead", 16, -8, 8,";truth #it{#eta}_{photon} ");
  histoStore()->createTH1F("newpt_truth_y_lead", 40, 0, 200,";truth #it{p}_{Tlead-photon} [GeV]");
  histoStore()->createTH1F("newpt_truth_y_sublead", 40, 0, 200,";truth #it{p}_{Tsub-photon} [GeV]");

  histoStore()->createTH1F("pt_bbyy", 40, 0, 200,";truth #it{p}_{Tbb#gamma#gamma} [GeV]");
  histoStore()->createTH1F("m_bbyy", 40, 100, 500,";truth #it{m}_{bb#gamma#gamma} [GeV]");
  histoStore()->createTH1F("prapidity_bbyy", 16, -8, 8,";truth #it{#eta}_{bb#gamma#gamma} ");
  histoStore()->createTH1F("rapidity_bbyy", 16, -8, 8,";truth #it{y}_{bb#gamma#gamma} ");

  histoStore()->createTH1F("pt_qqyy", 40, 0, 200,";truth #it{p}_{Tqq#gamma#gamma} [GeV]");
  histoStore()->createTH1F("m_qqyy", 40, 100, 500,";truth #it{m}_{qq#gamma#gamma} [GeV]");
  histoStore()->createTH1F("prapidity_qqyy", 16, -8, 8,";truth #it{#eta}_{qq#gamma#gamma} ");
  histoStore()->createTH1F("rapidity_qqyy", 16, -8, 8,";truth #it{y}_{qq#gamma#gamma} ");

  histoStore()->createTH1F("ggF_categories",11,-1,10,";ggF categories");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MCTRUTH0plots::execute()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  // Important to keep this, so that internal tools / event variables
  // are filled properly.
  HgammaAnalysis::execute();

  //float Lum=10.;
  //float xsec=43920;
  //float BR=0.00228;
  //float numberofevents=100000;
  //float w= (Lum*xsec*BR)/numberofevents;
  // float w=1.;

  /* 

  // 1. Grab the truth particles
  const xAOD::TruthParticleContainer *truthPtcls = 0;
  if ( event()->retrieve(truthPtcls, "TruthParticles" ).isFailure() ) {
    Error("execute()", "Failed to retrieve TruthParticle container." );
    return EL::StatusCode::FAILURE;
  }
 
  HG::TruthPtcls bquarkall(SG::VIEW_ELEMENTS);
  for (auto ptcl : *truthPtcls)if ( ptcl->status()==23 && (ptcl->pdgId()==5 ||ptcl->pdgId()==-5)  ) bquarkall.push_back(ptcl);
  
  
  HG::TruthPtcls quarkall(SG::VIEW_ELEMENTS);
  for (auto ptcl : *truthPtcls) if ( ptcl->status()==23 && MCUtils::PID::isQuark(ptcl->pdgId()) ) quarkall.push_back(ptcl);



  HG::TruthPtcls ysall(SG::VIEW_ELEMENTS); 
  HG::TruthPtcls esall(SG::VIEW_ELEMENTS); 
  HG::TruthPtcls musall(SG::VIEW_ELEMENTS);
  HG::TruthPtcls hads(SG::VIEW_ELEMENTS);
  HG::TruthPtcls photonsFromHiggsall(SG::VIEW_ELEMENTS); 
  HG::TruthPtcls HiggsDecay(SG::VIEW_ELEMENTS); 
  HG::TruthPtcls Bhadrons(SG::VIEW_ELEMENTS);
  HG::TruthPtcls Dhadrons(SG::VIEW_ELEMENTS);
  HG::TruthPtcls muonsFromBs(SG::VIEW_ELEMENTS);
  */
  //xAOD::TruthParticleContainer allphotons = truthHandler()->getPhotons();
  //xAOD::TruthParticleContainer selphotons = truthHandler()->applyPhotonSelection(allphotons);
  //xAOD::TruthParticleContainer higgs = truthHandler()->getHiggsBosons();  
  // xAOD::JetContainer jetsok = truthHandler()->getJets();
  //xAOD::JetContainer jets (SG::VIEW_ELEMENTS);
 HG::TruthParticles* ptcls = truthHandler()->getTruthParticles();


 // Truth particles
 xAOD::TruthParticleContainer allphotons   = truthHandler()->getPhotons();
 xAOD::TruthParticleContainer all_electrons = truthHandler()->getElectrons();
 xAOD::TruthParticleContainer all_muons     = truthHandler()->getMuons();
 xAOD::JetContainer           all_jets      = truthHandler()->getJets();
 xAOD::MissingETContainer     all_met       = truthHandler()->getMissingET();
 xAOD::TruthParticleContainer all_higgs     = truthHandler()->getHiggsBosons();
 
 // Apply fiducial selections to all containers
 xAOD::TruthParticleContainer selphotons   = truthHandler()->applyPhotonSelection   (allphotons);
 xAOD::TruthParticleContainer electrons = truthHandler()->applyElectronSelection (all_electrons);
 xAOD::TruthParticleContainer muons     = truthHandler()->applyMuonSelection     (all_muons);
 xAOD::JetContainer           jets      = truthHandler()->applyJetSelection      (all_jets);
 xAOD::JetContainer           bjets     = truthHandler()->applyBJetSelection     (jets);
 xAOD::MissingETContainer     met       = truthHandler()->applyMissingETSelection(all_met);
 
 // remove truth jets that are from electrons or photons
 truthHandler()->removeOverlap(selphotons, jets, electrons, muons);
 truthHandler()->passFiducial(&selphotons);

 // const xAOD::EventInfo *eventInfo = 0;
 //if (m_event->retrieve(eventInfo, "EventInfo").isFailure()){
 // HG::fatal("Cannot access EventInfo");
 // }
 
 //double w = eventInfo->mcEventWeights().at(151); 
 double w =eventHandler()->mcWeight();
 //std::cout<<w<<std::endl;
 /*
ysall = HG::getGoodTruthPhotons(truthPtcls);
  esall = HG::getGoodTruthElectrons(truthPtcls);
  musall = HG::getGoodTruthMuons(truthPtcls);
  hads  =  HG::getHadronsAndTheirDecay(truthPtcls);
  photonsFromHiggsall = HG::getPhotonsFromHiggs(truthPtcls);
  HiggsDecay = HG::getHiggsDecayProducts(truthPtcls);	   
  Bhadrons    = HG::getBHadrons(truthPtcls);
  Dhadrons    = HG::getDHadrons(truthPtcls);
  muonsFromBs = HG::getMuonsFromBs(truthPtcls);
  
  //pt cuts
  HG::TruthPtcls Bs(SG::VIEW_ELEMENTS);
  HG::TruthPtcls ys(SG::VIEW_ELEMENTS);
  HG::TruthPtcls es(SG::VIEW_ELEMENTS); 
  HG::TruthPtcls mus(SG::VIEW_ELEMENTS);
  HG::TruthPtcls photonsFromHiggs(SG::VIEW_ELEMENTS);
  HG::TruthPtcls quark(SG::VIEW_ELEMENTS);
  HG::TruthPtcls bquark(SG::VIEW_ELEMENTS);   

  Bs = HG::getBHadrons(truthPtcls,5.0*HG::GeV); 
  for (auto y : ysall )
    if (y->pt()>25*HG::GeV)ys.push_back(y);
  //for (auto l : jetsok )
  //if (l->pt()>30*HG::GeV)jets.push_back(l);
  for (auto e : esall )
    if (e->pt()>25*HG::GeV)es.push_back(e);
  for (auto m : musall )
    if (m->pt()>25*HG::GeV)mus.push_back(m);
  for (auto h : photonsFromHiggsall )
    if (h->pt()>25*HG::GeV)photonsFromHiggs.push_back(h);

  for (auto q : quarkall ){
    if (HG::minDRrap(q,ys)<0.4) continue;
    if (HG::minDRrap(q,es)<0.4) continue;
    if (q->pt()>25*HG::GeV)quark.push_back(q);
  }
  
  for (auto b : bquarkall ){
    if (HG::minDRrap(b,ys)<0.4) continue;
    if (HG::minDRrap(b,es)<0.4) continue;
    if (b->pt()>25*HG::GeV)bquark.push_back(b);
  }
 */
  //4.Fill histograms (truth)
 
  /*
  histoStore()->fillTH1F("ph_n",(ysall.size()),w);
  histoStore()->fillTH1F("ph_n25",(ys.size()),w); 
  histoStore()->fillTH1F("elec_n",(esall.size()),w);
  histoStore()->fillTH1F("elec_n25",(es.size()),w); 
  histoStore()->fillTH1F("muon_n",(musall.size()),w);
  histoStore()->fillTH1F("muon_n25",(mus.size()),w);
  histoStore()->fillTH1F("phHiggs_n",(photonsFromHiggsall.size()),w); 
  histoStore()->fillTH1F("phHiggs_n25",(photonsFromHiggs.size()),w); 
  histoStore()->fillTH1F("had_n",(hads.size()),w);
  histoStore()->fillTH1F("Bhad_n",(Bhadrons.size()),w);
  histoStore()->fillTH1F("Bhad_n5",(Bs.size()),w);
  histoStore()->fillTH1F("Dhad_n",(Dhadrons.size()),w);
  histoStore()->fillTH1F("muonsBs_n",(muonsFromBs.size()),w); 
  histoStore()->fillTH1F("quark_n",(quarkall.size()),w);
  histoStore()->fillTH1F("bquark_n",(bquarkall.size()),w);
  histoStore()->fillTH1F("quark_n25",(quark.size()),w);
  histoStore()->fillTH1F("bquark_n25",(bquark.size()),w); 

  std::sort(quark.begin(),quark.end(),ComparePt);
  std::sort(bquark.begin(),bquark.end(),ComparePt);


  for ( auto t_gam : ys ) {
    histoStore()->fillTH1F("pt_truth_photon",t_gam->pt()/HG::GeV,w);
    histoStore()->fillTH1F("eta_truth_photon",t_gam->eta(),w);
  }
  
 for ( auto q : quark ) {
          histoStore()->fillTH1F("pt_truth_q",q->pt()/HG::GeV,w);
          histoStore()->fillTH1F("prapidity_truth_q",q->eta(),w);
          histoStore()->fillTH1F("rapidity_truth_q",q->rapidity(),w);
       }

 for ( auto b : bquark ) {
	histoStore()->fillTH1F("pt_truth_b",b->pt()/HG::GeV,w);
	histoStore()->fillTH1F("prapidity_truth_b",b->eta(),w);
	histoStore()->fillTH1F("rapidity_truth_b",b->rapidity(),w);
      }

  //5. the Higgs pT and Rapidity (truth)
  if(ys.size()>1){
    //histoStore()->fillTH1F("pt_truth_y_lead",ys[0]->pt()/HG::GeV,w);
    histoStore()->fillTH1F("prapidity_truth_y_lead",ys[0]->eta(),w);
    histoStore()->fillTH1F("rapidity_truth_y_lead",ys[0]->rapidity(),w);
    TLorentzVector yy = ys[0]->p4()+ys[1]->p4();
    histoStore()->fillTH1F("pt_yy",yy.Pt()/HG::GeV,w);
    histoStore()->fillTH1F("m_yy",yy.M()/HG::GeV,w);
    histoStore()->fillTH1F("Rapidity_yy",yy.Rapidity(),w);
 histoStore()->fillTH1F("prapidity_yy",yy.Eta(),w);
  }  
  */

 int counterstable = 0;
 for(auto j : *ptcls){
   if(HG::isStable(j)) {
     counterstable++;
     // std::cout<<"it is stable"<<std::endl;
   }
 }
 histoStore()->fillTH1F("newstable_n",counterstable,w);
 histoStore()->fillTH1F("newph_n",(selphotons.size()),w);
 histoStore()->fillTH1F("newH_n",(all_higgs.size()),w);
 histoStore()->fillTH1F("newjets_n",(jets.size()),w);
 
 
 
 if(selphotons.size()>1){
    
   TLorentzVector newyy = selphotons[0]->p4()+selphotons[1]->p4();
   double myy = newyy.M();
   if ((selphotons[0]->pt()/myy >= 0.35) && (selphotons[1]->pt()/myy >= 0.25) ){
     if( myy > 105*HG::GeV && myy < 160*HG::GeV ){ 
       // ------------------------------------
       // ------------------------------------
       //    Plot coupling categories here
       // ------------------------------------
       // ------------------------------------
       
       histoStore()->fillTH1F("newpt_truth_y_lead",selphotons[0]->pt()/HG::GeV,w);
       histoStore()->fillTH1F("newpt_truth_y_sublead",selphotons[1]->pt()/HG::GeV,w);
       
       histoStore()->fillTH1F("newpt_yy",newyy.Pt()/HG::GeV,w);
       

       if(all_higgs.size()>=1){
	 histoStore()->fillTH2F("newstable_ptH",all_higgs[0]->pt()/HG::GeV,counterstable,w);
	 
	 histoStore()->fillTH1F("newpt_H",all_higgs[0]->pt()/HG::GeV,w);
       }

       if(jets.size()>=1){
	 histoStore()->fillTH1F("newpt_jets",jets[0]->pt()/HG::GeV,w);
       }
       
     }
   }
    

 }




 /*
 if(hads.size()>=1){
   float sum =0;
   for ( auto q : hads ) {
     sum += q->pt()/HG::GeV;
   }
   histoStore()->fillTH1F("newpt_hads",sum,w);
   if(all_higgs.size()>=1) histoStore()->fillTH2F("newpthads_ptH",all_higgs[0]->pt()/HG::GeV,sum,w);
 }
 

    //6. Fill truth jets histograms
    if (quark.size()>1){
      
       histoStore()->fillTH1F("pt_truth_q_lead",quark[0]->pt()/HG::GeV,w);
       histoStore()->fillTH1F("prapidity_truth_q_lead",quark[0]->eta(),w);
       histoStore()->fillTH1F("rapidity_truth_q_lead",quark[0]->rapidity(),w);
       TLorentzVector q1=quark[0]->p4(), q2 = quark[1]->p4();
       histoStore()->fillTH1F("drapidity_truth_q",fabs(q1.Rapidity()-q2.Rapidity()),w);
       histoStore()->fillTH1F("deltaphi_truth_q",fabs(q1.DeltaPhi(q2)),w);
       histoStore()->fillTH1F("deltar_truth_q",q1.DeltaR(q2),w);
       TLorentzVector qq = q1+q2;
       histoStore()->fillTH1F("pt_qq",qq.Pt()/HG::GeV,w);
       histoStore()->fillTH1F("m_qq",qq.M()/HG::GeV,w);
	histoStore()->fillTH1F("prapidity_qq",qq.Eta(),w);
       histoStore()->fillTH1F("rapidity_qq",qq.Rapidity(),w);
      
    }
    
   

     if (bquark.size()>1){
     
      histoStore()->fillTH1F("pt_truth_b_lead",bquark[0]->pt()/HG::GeV,w);
      histoStore()->fillTH1F("prapidity_truth_b_lead",bquark[0]->eta(),w);
      histoStore()->fillTH1F("rapidity_truth_b_lead",bquark[0]->rapidity(),w);
      TLorentzVector bq1=bquark[0]->p4(), bq2 = bquark[1]->p4();
      histoStore()->fillTH1F("drapidity_truth_b",fabs(bq1.Rapidity()-bq2.Rapidity()),w);
      histoStore()->fillTH1F("deltaphi_truth_b",fabs(bq1.DeltaPhi(bq2)),w);
      histoStore()->fillTH1F("deltar_truth_b",bq1.DeltaR(bq2),w);
      TLorentzVector bb = bq1+bq2;
      histoStore()->fillTH1F("pt_bb",bb.Pt()/HG::GeV,w);
      histoStore()->fillTH1F("m_bb",bb.M()/HG::GeV,w);
   	histoStore()->fillTH1F("prapidity_bb",bb.Eta(),w);
       histoStore()->fillTH1F("rapidity_bb",bb.Rapidity(),w);
    }
    
   if(ys.size()>1 && quark.size()>1){
	 TLorentzVector yy = ys[0]->p4()+ys[1]->p4();
	 TLorentzVector q1=quark[0]->p4(), q2 = quark[1]->p4();
	 TLorentzVector qq = q1+q2;
 	histoStore()->fillTH2F("m_yy_vs_m_qq",qq.M()/HG::GeV,yy.M()/HG::GeV,w);
	TLorentzVector yyqq = qq+yy;
 	histoStore()->fillTH1F("pt_qqyy",yyqq.Pt()/HG::GeV,w);
       histoStore()->fillTH1F("m_qqyy",yyqq.M()/HG::GeV,w);
	histoStore()->fillTH1F("prapidity_qqyy",yyqq.Eta(),w);
       histoStore()->fillTH1F("rapidity_qqyy",yyqq.Rapidity(),w);
	
   }

 if(ys.size()>1 && bquark.size()>1){
	 TLorentzVector yy = ys[0]->p4()+ys[1]->p4();
	 TLorentzVector bq1=bquark[0]->p4(), bq2 = bquark[1]->p4();
	 TLorentzVector bb = bq1+bq2;
 	histoStore()->fillTH2F("m_yy_vs_m_bb",bb.M()/HG::GeV,yy.M()/HG::GeV,w);
	TLorentzVector yybb = bb+yy;
 	histoStore()->fillTH1F("pt_bbyy",yybb.Pt()/HG::GeV,w);
       histoStore()->fillTH1F("m_bbyy",yybb.M()/HG::GeV,w);
	histoStore()->fillTH1F("prapidity_bbyy",yybb.Eta(),w);
       histoStore()->fillTH1F("rapidity_bbyy",yybb.Rapidity(),w);
	
   }
 
 
 */
  
  return EL::StatusCode::SUCCESS;
}

