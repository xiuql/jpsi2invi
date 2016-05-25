// -*- C++ -*-
//
//
// Description: Psi' -> pi+pi-Jpsi,J/psi -> mu+mu- 
//
// Original Author:  Chen Qian <chenq@pku.edu.cn>
//         Created:  [2016-05-23 Mon 14:40] 
//         Inspired by Shi xin's code 
//


//
// system include files
//


#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/LoadFactoryEntries.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"

#include "EventModel/EventHeader.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/Helix.h"
#include "VertexFit/WTrackParameter.h"
#include "VertexFit/VertexFit.h"

#include "ParticleID/ParticleID.h"

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

//
// class declaration
//

class Mumu : public Algorithm {
  
public:
  Mumu(const std::string&, ISvcLocator*);
  ~Mumu(); 
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private:
  // declare r0, z0 cut for charged tracks
  double m_ecms; 
  double m_vr0cut, m_vz0cut;
  double m_distin_pionlep;
  double m_cha_costheta_cut; 
  double m_total_number_of_charged_max; 
  double m_min_emctime;
  double m_max_emctime;
  double m_gammaCosCut; 
  double m_costheta_barrel_max;
  double m_costheta_endcap_min;
  double m_costheta_endcap_max;
  double m_energy_barrel_min;
  double m_energy_endcap_min;
  double m_photon_iso_angle_min;
  double m_pion_polar_angle_max;
  double m_pion_momentum_max;
  double m_prob_pion_min;
  double m_dipion_mass_min; 
  double m_dipion_mass_max;
  double m_pipi_costheta_max;
  double m_pipisys_costheta_max; 
  double m_mumu_mass_min; 
  double m_mumu_mass_max;
  double m_ep_ratio;

  // output file
  std::string m_output_filename;
  TFile* m_fout; 
  
  // define Histograms
  TH1F* h_evtflw; 
  
  // define Trees
  
  // signal 
  TTree* m_tree;

  // common info 
  int m_run;
  int m_event;

  // charged tracks
  int m_ncharged;
  int m_nptrk;
  int m_nmtrk;
  
  // neutral tracks
  int m_nshow;
  int m_ngam;

  // vertex 
  double m_vr0;
  double m_vz0;
  double m_vtx_mrecpipi; // pipi invariant mass
  double m_vtx_mmumu;

  // check MDC and EMC match
  long m_pion_matched;
  long m_lep_matched;

  // jpsi2invi
  int m_ntrk; 
  int m_npho; 
  
  //ncut
  int Ncut0;
  int Ncut1;
  int Ncut2;
  int Ncut3;
  int Ncut4;
  int Ncut5;
  int Ncut6;
  int Ncut7;
  int Ncut8;
  int Ncut9;
  int Ncut10;
  int Ncut11;
  int Ncut12;


  
  // functions
  void book_histogram();
  void book_tree(); 
  void buildJpsiToMumu();
  int selectChargedTracks(SmartDataPtr<EvtRecEvent>,
			  SmartDataPtr<EvtRecTrackCol>,
			  std::vector<int> &,
			  std::vector<int> &); 
  int selectPionPlusPionMinus(SmartDataPtr<EvtRecTrackCol>,
			      std::vector<int>,
			      std::vector<int>);
  int selectMuPlusMuMinus(SmartDataPtr<EvtRecTrackCol>,
			      std::vector<int>,
			      std::vector<int>);
  void calcTrackPID(EvtRecTrackIterator,
		    double& ,
		    double&);
  bool hasGoodPiPiVertex(RecMdcKalTrack *,
			 RecMdcKalTrack *,
			 bool);
  bool hasGoodMuMuVertex(RecMdcKalTrack *,
			 RecMdcKalTrack *,
			 bool);
  int selectNeutralTracks(SmartDataPtr<EvtRecEvent>,
			  SmartDataPtr<EvtRecTrackCol>);
  bool passVertexSelection(CLHEP::Hep3Vector,
			   RecMdcKalTrack* ); 
  CLHEP::Hep3Vector getOrigin();
}; 


//
// module declare
//

DECLARE_ALGORITHM_FACTORY( Mumu )
DECLARE_FACTORY_ENTRIES( Mumu ) {
  DECLARE_ALGORITHM(Mumu);
}

LOAD_FACTORY_ENTRIES( Mumu )

//
// constants
//

const double PION_MASS = 0.139570;
const double MU_MASS = 0.105658;


//
// member functions
//
  
Mumu::Mumu(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator) {
  std::cout<<"test!!!!!!!!!"<<endl;
  declareProperty("OutputFileName", m_output_filename);
  declareProperty("Ecms", m_ecms = 3.686);
  declareProperty("Vr0cut", m_vr0cut=1.0);
  declareProperty("Vz0cut", m_vz0cut=10.0);
  declareProperty("DiffPionLep", m_distin_pionlep=0.8);
  declareProperty("ChaCosthetaCut", m_cha_costheta_cut=0.93);
  declareProperty("TotalNumberOfChargedMax", m_total_number_of_charged_max=50);
  declareProperty("MinEstCut", m_min_emctime=0.0);
  declareProperty("MaxEstCut", m_max_emctime=14.0);
  declareProperty("GammaCosCut",  m_gammaCosCut= 0.93); 
  declareProperty("CosthetaBarrelMax", m_costheta_barrel_max=0.8);
  declareProperty("CosthetaEndcapMin", m_costheta_endcap_min=0.86);
  declareProperty("CosthetaEndcapMax", m_costheta_endcap_max=0.92);
  declareProperty("EnergyBarrelMin", m_energy_barrel_min=0.025); 
  declareProperty("EnergyEndcapMin", m_energy_endcap_min=0.050); 
  declareProperty("PhotonIsoAngleMin", m_photon_iso_angle_min=10);
  declareProperty("PionPolarAngleMax", m_pion_polar_angle_max=0.8);
  declareProperty("PionMomentumMax", m_pion_momentum_max=0.45); 
  declareProperty("ProbPionMin", m_prob_pion_min=0.001);
  declareProperty("DipionMassMin", m_dipion_mass_min=3.0); 
  declareProperty("DipionMassMax", m_dipion_mass_max=3.2); 
  declareProperty("PiPiCosthetaMax", m_pipi_costheta_max=0.95);
  declareProperty("PiPiSysCosthetaMax", m_pipisys_costheta_max=0.90);
  declareProperty("DipionMassMin", m_mumu_mass_min=3.0); 
  declareProperty("DipionMassMax", m_mumu_mass_max=3.2); 
  declareProperty("epratio", m_ep_ratio=0.26);
}


StatusCode Mumu::initialize(){
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << ">>>>>>> in initialize()" << endmsg;
  
    //ncut
  Ncut0=0;
  Ncut1=0;
  Ncut2=0;
  Ncut3=0;
  Ncut4=0;
  Ncut5=0;
  Ncut6=0;
  Ncut7=0;
  Ncut8=0;
  Ncut9=0;
  Ncut10=0;
  Ncut11=0;
  Ncut12=0;
  


  m_fout = new TFile(m_output_filename.c_str(), "RECREATE");
  m_fout->cd(); 

  book_histogram(); 
  book_tree(); 

  log << MSG::INFO << "successfully return from initialize()" <<endmsg;
  return StatusCode::SUCCESS;
}


StatusCode Mumu::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;
 
  h_evtflw->Fill(0); // raw 
  Ncut0++;
  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
  if(!eventHeader) return StatusCode::FAILURE;

  m_run = eventHeader->runNumber();
  m_event = eventHeader->eventNumber();
  
  buildJpsiToMumu();
  m_tree->Fill();

  return StatusCode::SUCCESS; 
}

StatusCode Mumu::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in finalize()" << endmsg;

  m_fout->cd();
  m_tree->Write();
  h_evtflw->Write();
  m_fout->Close();
  
  std::cout<<"raw:                              "<<Ncut0<<endl;
  std::cout<<"N_{Good} = 4:                     "<<Ncut1<<endl;
  std::cout<<"good photon:                      "<<Ncut2<<endl;
  std::cout<<"|cos#theta|<0.8:                  "<<Ncut3<<endl;
  std::cout<<"|p|<0.45:                         "<<Ncut4<<endl;
  std::cout<<"pass PID:                         "<<Ncut5<<endl;
  std::cout<<"cos#theta_{#pi^{+}#pi^{-}}<0.95:  "<<Ncut6<<endl;
  std::cout<<"cos#theta_{#pi#pi sys}<0.9:       "<<Ncut7<<endl;
  std::cout<<"3<M_{#pi#pi}^{rec}<3.2:           "<<Ncut8<<endl;
  std::cout<<"pi^+pi^-:                         "<<Ncut9<<endl;
  std::cout<<"[E/p]^+<0.26&&[E/p]^-<0.26:       "<<Ncut10<<endl;
  std::cout<<"3<M_{#mu#mu}<3.2:                 "<<Ncut11<<endl;
  std::cout<<"mu^+mu^-:                         "<<Ncut12<<endl;

  
  return StatusCode::SUCCESS;
}


Mumu::~Mumu() {
}


void Mumu::book_histogram() {
  h_evtflw = new TH1F("hevtflw", "eventflow", 10, 0, 10);
  if (!h_evtflw) return;
  h_evtflw->GetXaxis()->SetBinLabel(1, "raw");
  h_evtflw->GetXaxis()->SetBinLabel(2, "N_{Good}=2");
  h_evtflw->GetXaxis()->SetBinLabel(3, "N_{#gamma}=0");
  h_evtflw->GetXaxis()->SetBinLabel(4, "|cos#theta|<0.8");
  h_evtflw->GetXaxis()->SetBinLabel(5, "|p|<0.45");
  h_evtflw->GetXaxis()->SetBinLabel(6, "PID"); 
  h_evtflw->GetXaxis()->SetBinLabel(7, "cos#theta_{#pi^{+}#pi^{-}}<0.95");
  h_evtflw->GetXaxis()->SetBinLabel(8, "cos#theta_{#pi#pi sys}<0.9");
  h_evtflw->GetXaxis()->SetBinLabel(9, "3<M_{#pi#pi}^{rec}<3.2");
  h_evtflw->GetXaxis()->SetBinLabel(10, "pi^+pi^-"); 
  h_evtflw->GetXaxis()->SetBinLabel(11, "[E/p]^+<0.26&&[E/p]^-<0.26");
  h_evtflw->GetXaxis()->SetBinLabel(12, "3<M_{#mu#mu}<3.2");
  h_evtflw->GetXaxis()->SetBinLabel(13, "mu^+mu^-");  
}


void Mumu::book_tree() {
  m_tree=new TTree("signal", "signal");
  if (!m_tree) return; 

  //commom info
  m_tree->Branch("run",&m_run,"run/I");
  m_tree->Branch("event",&m_event,"event/I");
	  
  //charged tracks
  m_tree->Branch("ncharged",&m_ncharged,"ncharged/I");
  m_tree->Branch("nptrk",&m_nptrk,"nptrk/I");
  m_tree->Branch("nmtrk",&m_nmtrk,"nmtrk/I");
	  
  //vertex
  m_tree->Branch("vr0",&m_vr0,"vr0/D");
  m_tree->Branch("vz0",&m_vz0,"vz0/D");
  m_tree->Branch("vtx_mrecpipi",&m_vtx_mrecpipi,"vtx_mrecpipi/D");
  m_tree->Branch("vtx_mmumu",&m_vtx_mmumu,"vtx_mmumu/D");
  
	  
  //netual tracks
  m_tree->Branch("nshow",&m_nshow,"nshow/I");
  m_tree->Branch("ngam",&m_ngam,"ngam/I");

  
}


void Mumu::buildJpsiToMumu() {

  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(),"/Event/EvtRec/EvtRecEvent");
  if(!evtRecEvent) return;

  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  if(!evtRecTrkCol) return;

  std::vector<int> iPGood, iMGood; 
  selectChargedTracks(evtRecEvent, evtRecTrkCol, iPGood, iMGood);

  if (m_ncharged != 4) return;
  h_evtflw->Fill(1); // N_{Good} = 4 
  Ncut1++;
  
  selectNeutralTracks(evtRecEvent, evtRecTrkCol);
  //if (m_ngam != 0) return;
  h_evtflw->Fill(2); // good photon
  Ncut2++;
  
  if(selectPionPlusPionMinus(evtRecTrkCol, iPGood, iMGood)!=2) return;
  h_evtflw->Fill(9);               //pipi
  Ncut9++;

  if(selectMuPlusMuMinus(evtRecTrkCol, iPGood, iMGood)!=2) return;
  h_evtflw->Fill(12);              //mumu
  Ncut12++;

}


CLHEP::Hep3Vector Mumu::getOrigin() {
  CLHEP::Hep3Vector xorigin(0,0,0);
  IVertexDbSvc*  vtxsvc;
  Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
  if(vtxsvc->isVertexValid()){
    double *dbv = vtxsvc->PrimaryVertex(); 
    xorigin.setX(dbv[0]);
    xorigin.setY(dbv[1]);
    xorigin.setZ(dbv[2]);
  }
  return xorigin; 
}


bool Mumu::passVertexSelection(CLHEP::Hep3Vector xorigin,
				    RecMdcKalTrack* mdcTrk ) {
  HepVector a = mdcTrk->helix();
  HepSymMatrix Ea = mdcTrk->err();
  HepPoint3D point0(0.,0.,0.);
  VFHelix helixip(point0,a,Ea);
  HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
  helixip.pivot(IP);
  HepVector vecipa = helixip.a();
  
  m_vz0 = vecipa[3];
  m_vr0 = vecipa[0];
  
  if(fabs(m_vz0) >= m_vz0cut) return false;
  if(fabs(m_vr0) >= m_vr0cut) return false;
  
  return true;
}


int Mumu::selectChargedTracks(SmartDataPtr<EvtRecEvent> evtRecEvent,
				   SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
				   std::vector<int> & iPGood,
				   std::vector<int> & iMGood) {

  CLHEP::Hep3Vector xorigin = getOrigin(); 

  std::vector<int> iGood;
  iGood.clear();
  iPGood.clear();
  iMGood.clear();
  
  // loop through charged tracks 
  for(int i = 0; i < evtRecEvent->totalCharged(); i++){
    // get mdcTrk 
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;

    // Good Kalman Track 
    if(!(*itTrk)->isMdcKalTrackValid()) continue;

    if(!(*itTrk)->isMdcTrackValid()) continue;
    RecMdcKalTrack* mdcTrk = (*itTrk)->mdcKalTrack();

    // Good Vertex 
    if (!passVertexSelection(xorigin, mdcTrk) ) continue; 

    // Polar angle cut
    if(fabs(cos(mdcTrk->theta())) > m_cha_costheta_cut) continue;

    iGood.push_back((*itTrk)->trackId());
    if(mdcTrk->charge()>0) iPGood.push_back((*itTrk)->trackId());
    if(mdcTrk->charge()<0) iMGood.push_back((*itTrk)->trackId());

  } // end charged tracks

  m_ncharged = iGood.size();
  m_nptrk = iPGood.size();
  m_nmtrk = iMGood.size(); 
  
  return iGood.size(); 
}

int Mumu::selectPionPlusPionMinus(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
				       std::vector<int> iPGood,
				       std::vector<int> iMGood) {
  int npipi = 0;
  bool evtflw_filled = false;
  
  for(unsigned int i1 = 0; i1 < iPGood.size(); i1++) {
    EvtRecTrackIterator itTrk_p = evtRecTrkCol->begin() + iPGood[i1];
    RecMdcTrack* mdcTrk_p = (*itTrk_p)->mdcTrack();
    if (mdcTrk_p->charge() < 0) continue; // only positive charged tracks

    for(unsigned int i2 = 0; i2 < iMGood.size(); i2++) {
      EvtRecTrackIterator itTrk_m = evtRecTrkCol->begin() + iMGood[i2];
      RecMdcTrack* mdcTrk_m = (*itTrk_m)->mdcTrack();
      if (mdcTrk_m->charge() > 0) continue; // only negative charged tracks

      // polar angle for both pions
      if ( ! ( fabs(cos(mdcTrk_p->theta())) < m_pion_polar_angle_max &&
      	       fabs(cos(mdcTrk_m->theta())) < m_pion_polar_angle_max )) continue;
      if ( !evtflw_filled ) {
		  h_evtflw->Fill(3); // |cos#theta|<0.8
		  Ncut3++;
	  }

      // pion momentum
      if ( ! ( fabs(mdcTrk_p->p()) < m_pion_momentum_max  &&
      	       fabs(mdcTrk_m->p()) < m_pion_momentum_max )) continue;

      if ( !evtflw_filled ) {
		  h_evtflw->Fill(4); //|p|<0.45 
		  Ncut4++;
	  }
      
      // track PID
      double prob_pip, prob_kp, prob_pim, prob_km; 
      calcTrackPID(itTrk_p, prob_pip, prob_kp);  
      calcTrackPID(itTrk_m, prob_pim, prob_km);
      // printf(">>> %f, %f, %f, %f \n", prob_pip, prob_kp, prob_pim, prob_km);

      if(! (prob_pip > prob_kp &&
	    prob_pip > m_prob_pion_min &&
	    prob_pim > prob_km &&
	    prob_pim > m_prob_pion_min) ) continue;

      if ( !evtflw_filled ) {
		  h_evtflw->Fill(5); //PID
		  Ncut5++;
	  }
 
      // apply vertex fit
      RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin()+iPGood[i1]))->mdcKalTrack();
      RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin()+iMGood[i2]))->mdcKalTrack();
	    
      if (! hasGoodPiPiVertex(pipTrk, pimTrk, evtflw_filled) ) continue; 
      
      npipi++;
      evtflw_filled = true;
    }
  } 

  return npipi; 
}


int Mumu::selectMuPlusMuMinus(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
				       std::vector<int> iPGood,
				       std::vector<int> iMGood) {
  int nmumu = 0;
  bool evtflw_filled = false;

  for(unsigned int i1 = 0; i1 < iPGood.size(); i1++) {
    EvtRecTrackIterator itTrk_p = evtRecTrkCol->begin() + iPGood[i1];
    RecMdcTrack* mdcTrk_p = (*itTrk_p)->mdcTrack();

	RecEmcShower *emcTrk_p = (*itTrk_p)->emcShower();
   
   if (mdcTrk_p->charge() < 0) continue; // only positive charged tracks

    for(unsigned int i2 = 0; i2 < iMGood.size(); i2++) {
      EvtRecTrackIterator itTrk_m = evtRecTrkCol->begin() + iMGood[i2];

      RecMdcTrack* mdcTrk_m = (*itTrk_m)->mdcTrack();
	if(!((*itTrk_m)->isEmcShowerValid())||!((*itTrk_p)->isEmcShowerValid())) continue; 
   	  RecEmcShower *emcTrk_m = (*itTrk_m)->emcShower();
      if (mdcTrk_m->charge() > 0) continue; // only negative charged tracks

	  // ratio cut
	  double emcenergy_p = emcTrk_p->energy(); 	  
	  double emcenergy_m = emcTrk_m->energy(); 
	  double pemc_p=fabs(mdcTrk_p->p());
	  double pemc_m=fabs(mdcTrk_m->p());
	  if(emcenergy_p/pemc_p >= m_ep_ratio && emcenergy_m/pemc_m >= m_ep_ratio) continue;

	  h_evtflw->Fill(10);     //epratio cut
	  Ncut10++;

      // apply vertex fit
      RecMdcKalTrack *mupTrk = (*(evtRecTrkCol->begin()+iPGood[i1]))->mdcKalTrack();
      RecMdcKalTrack *mumTrk = (*(evtRecTrkCol->begin()+iMGood[i2]))->mdcKalTrack();
	    
      if (! hasGoodMuMuVertex(mupTrk, mumTrk, evtflw_filled) ) continue; 
      
      nmumu++;
      evtflw_filled = true;
    }
  } 


  return nmumu; 
}

void Mumu::calcTrackPID(EvtRecTrackIterator itTrk_p,
			     double& prob_pip,
			     double& prob_kp) {
  prob_pip = 999.; 
  prob_kp = 999.; 
  ParticleID * pidp = ParticleID::instance();
  pidp->init();
  pidp->setMethod(pidp->methodProbability());
  pidp->setChiMinCut(4);
  pidp->setRecTrack(*itTrk_p);
  // use PID sub-system
  pidp->usePidSys(pidp->useDedx() | pidp->useTof1() | pidp->useTof2());
  pidp->identify(pidp->onlyPionKaonProton());
  pidp->calculate();
  if(pidp->IsPidInfoValid()) {
    prob_pip = pidp->probPion();
    prob_kp  = pidp->probKaon();
  }
}

bool Mumu::hasGoodPiPiVertex(RecMdcKalTrack *pipTrk,
				  RecMdcKalTrack *pimTrk,
				  bool evtflw_filled) {

  HepLorentzVector pcms(0.011*m_ecms, 0., 0., m_ecms);

  HepLorentzVector p4_vtx_pip, p4_vtx_pim,p4_vtx_recpipi;
  WTrackParameter wvpipTrk, wvpimTrk;
  pipTrk->setPidType(RecMdcKalTrack::pion);
  wvpipTrk = WTrackParameter(PION_MASS, pipTrk->getZHelix(), pipTrk->getZError());

  pimTrk->setPidType(RecMdcKalTrack::pion);
  wvpimTrk = WTrackParameter(PION_MASS, pimTrk->getZHelix(), pimTrk->getZError());
  
  HepPoint3D vx(0., 0., 0.);
  HepSymMatrix Evx(3, 0);

  double bx = 1E+6;
  double by = 1E+6;
  double bz = 1E+6;
  Evx[0][0] = bx*bx;
  Evx[1][1] = by*by;
  Evx[2][2] = bz*bz;
  
  VertexParameter vxpar;
  vxpar.setVx(vx);
  vxpar.setEvx(Evx);
  
  VertexFit* vtxfit = VertexFit::instance();
  vtxfit->init();
  vtxfit->AddTrack(0,  wvpipTrk);
  vtxfit->AddTrack(1,  wvpimTrk);
  vtxfit->AddVertex(0, vxpar,0,1);

  if(!vtxfit->Fit(0)) return false;

  vtxfit->Swim(0);
      
  WTrackParameter wpip = vtxfit->wtrk(0);
  WTrackParameter wpim = vtxfit->wtrk(1);
  p4_vtx_pip = vtxfit->pfit(0) ;
  p4_vtx_pim = vtxfit->pfit(1) ;
  p4_vtx_recpipi = pcms - p4_vtx_pip - p4_vtx_pim;

  double cospipi = cos(p4_vtx_pip.vect().angle(p4_vtx_pim.vect()));
  double cos2pisys = (p4_vtx_pip + p4_vtx_pim).cosTheta();

  if( ! (cospipi < m_pipi_costheta_max) ) return false;
  if( !evtflw_filled ) {
	  h_evtflw->Fill(6); // "cos#theta_{#pi^{+}#pi^{-}}<0.95"
	  Ncut6++;
  }

  if( ! (fabs(cos2pisys) < m_pipisys_costheta_max ) ) return false;
  if( !evtflw_filled ){
	  h_evtflw->Fill(7); // cos#theta_{#pi#pi sys}<0.9
	  Ncut7++;
  }	  

  if( ! ( p4_vtx_recpipi.m() >= m_dipion_mass_min &&
	  p4_vtx_recpipi.m() <= m_dipion_mass_max) ) return false;
	  std::cout<<p4_vtx_recpipi.m()<<endl;
  if( !evtflw_filled ){
	  h_evtflw->Fill(8); // 3<M_{#pi#pi}^{rec}<3.2
	  Ncut8++;
  }
  
  
  m_vtx_mrecpipi = p4_vtx_recpipi.m();
  std::cout<<m_vtx_mrecpipi<<endl;
  return true;
}

bool Mumu::hasGoodMuMuVertex(RecMdcKalTrack *mupTrk,
				  RecMdcKalTrack *mumTrk,
				  bool evtflw_filled) {

  HepLorentzVector p4_vtx_mup, p4_vtx_mum,p4_vtx_mumu;
  WTrackParameter wvmupTrk, wvmumTrk;
  mupTrk->setPidType(RecMdcKalTrack::muon);

  wvmupTrk = WTrackParameter(MU_MASS, mupTrk->getZHelix(), mupTrk->getZError());

  mumTrk->setPidType(RecMdcKalTrack::muon);
  wvmumTrk = WTrackParameter(MU_MASS, mumTrk->getZHelix(), mumTrk->getZError());
  
  HepPoint3D vx(0., 0., 0.);
  HepSymMatrix Evx(3, 0);

  double bx = 1E+6;
  double by = 1E+6;
  double bz = 1E+6;
  Evx[0][0] = bx*bx;
  Evx[1][1] = by*by;
  Evx[2][2] = bz*bz;

  VertexParameter vxpar;
  vxpar.setVx(vx);
  vxpar.setEvx(Evx);
  
  VertexFit* vtxfit = VertexFit::instance();
  vtxfit->init();
  vtxfit->AddTrack(0,  wvmupTrk);
  vtxfit->AddTrack(1,  wvmumTrk);
  vtxfit->AddVertex(0, vxpar,0,1);

  //if(!vtxfit->Fit(0)) return false;

  //vtxfit->Swim(0);

  WTrackParameter wmup = vtxfit->wtrk(0);
  WTrackParameter wmum = vtxfit->wtrk(1);

  p4_vtx_mup = vtxfit->pfit(0) ;
  p4_vtx_mum = vtxfit->pfit(1) ;
  p4_vtx_mumu = p4_vtx_mup + p4_vtx_mum;
  
  if(  p4_vtx_mumu.m() < m_mumu_mass_min ||
	  p4_vtx_mumu.m() > m_mumu_mass_max) return false;
	  std::cout<<p4_vtx_mumu.m()<<endl;
  if( !evtflw_filled ){
	  h_evtflw->Fill(11); // 3<M_{#mu#mu}<3.2
	  Ncut11++;
  }

  m_vtx_mmumu = p4_vtx_mumu.m();
  
  return true;
}


int Mumu::selectNeutralTracks(SmartDataPtr<EvtRecEvent> evtRecEvent,
				   SmartDataPtr<EvtRecTrackCol> evtRecTrkCol) {

  std::vector<int> iGam;
  iGam.clear();
  std::vector<int> iShow;
  iShow.clear();

  // loop through neutral tracks
  for(int i=evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++) {
    if (i > m_total_number_of_charged_max) break;

    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i ;
    if(!(*itTrk)->isEmcShowerValid()) continue;
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    iShow.push_back((*itTrk)->trackId());
    
    // TDC window
    if ( !(emcTrk->time() >= m_min_emctime && emcTrk->time() <= m_max_emctime) )
      continue; 

    // Energy threshold
    double abs_costheta(fabs(cos(emcTrk->theta())));
    bool barrel = (abs_costheta < m_costheta_barrel_max); 
    bool endcap = (abs_costheta > m_costheta_endcap_min
		   && abs_costheta < m_costheta_barrel_max);
    double eraw = emcTrk->energy();
    
    if ( !( (barrel && eraw > m_energy_barrel_min)
	    || (endcap && eraw > m_energy_endcap_min)))  continue; 

    // photon isolation: the opening angle between a candidate shower
    // and the closest charged track should be larger than 10 degree 
    CLHEP::Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

    // EMC costheta cut 
    double costhe = cos(emcpos.theta());
    if ( fabs(costhe) >= m_gammaCosCut) continue;
    
    // find the nearest charged track
    double dthe = 200.;
    double dphi = 200.;
    double dang = 200.; 
    for(int j = 0; j < evtRecEvent->totalCharged(); j++) {
      EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
      if(!(*jtTrk)->isExtTrackValid()) continue;
      RecExtTrack *extTrk = (*jtTrk)->extTrack();
      if(extTrk->emcVolumeNumber() == -1) continue;
      CLHEP::Hep3Vector extpos = extTrk->emcPosition();
      double angd = extpos.angle(emcpos);
      double thed = extpos.theta() - emcpos.theta();
      double phid = extpos.deltaPhi(emcpos);
      thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
      phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;

      if(fabs(thed) < fabs(dthe)) dthe = thed;
      if(fabs(phid) < fabs(dphi)) dphi = phid;
      if(angd < dang) dang = angd;	    
    }

    if(dang>=200) continue;
    dthe = dthe * 180 / (CLHEP::pi);
    dphi = dphi * 180 / (CLHEP::pi);
    dang = dang * 180 / (CLHEP::pi);
    if (dang < m_photon_iso_angle_min ) continue; 

    iGam.push_back((*itTrk)->trackId());
  } // end loop neutral tracks     

  m_ngam = iGam.size();
  m_nshow = iShow.size();

  return iGam.size(); 
}
