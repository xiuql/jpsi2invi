// -*- C++ -*-
//
//
// Description: J/psi -> Invisible 
//
// Original Author:  SHI Xin <shixin@ihep.ac.cn>
//         Created:  [2016-03-23 Wed 09:12] 
//         Inspired by Zhu Kai and Zhang Chi's code 
//


//
// system include files
//


// #include "AIDA/IHistogram1D.h"

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/LoadFactoryEntries.h"
#include "GaudiKernel/NTuple.h"
// #include "GaudiKernel/IHistogramSvc.h"
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

//
// class declaration
//

class Jpsi2invi : public Algorithm {
  
public:
  Jpsi2invi(const std::string&, ISvcLocator*);
  ~Jpsi2invi(); 
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

  // output file
  std::string m_output_filename;
  TFile* m_fout; 
  
  // define Histograms
  // IHistogram1D *h_evtflw; 
  TH1F* h_evtflw; 
  
  // define Ntuples
  
  // signal 
  NTuple::Tuple* m_tuple1;

  // common info 
  NTuple::Item<int> m_run;
  NTuple::Item<int> m_event;

  // charged tracks
  NTuple::Item<int> m_ncharged;
  NTuple::Item<int> m_nptrk;
  NTuple::Item<int> m_nmtrk;
  
  // neutral tracks
  NTuple::Item<int> m_nshow;
  NTuple::Item<int> m_ngam;

  // vertex 
  NTuple::Item<double> m_vr0;
  NTuple::Item<double> m_vz0;
  NTuple::Item<double> m_vtx_mrecpipi; // pipi invariant mass

  // check MDC and EMC match
  NTuple::Item<long> m_pion_matched;
  NTuple::Item<long> m_lep_matched;

  // jpsi2invi
  NTuple::Item<int> m_ntrk; 
  NTuple::Item<int> m_npho; 
  
  // functions
  void book_histogram();
  bool book_ntuple_signal(MsgStream); 
  void buildJpsiToInvisible();
  int selectChargedTracks(SmartDataPtr<EvtRecEvent>,
			  SmartDataPtr<EvtRecTrackCol>,
			  std::vector<int> &,
			  std::vector<int> &); 
  int selectPionPlusPionMinus(SmartDataPtr<EvtRecTrackCol>,
			      std::vector<int>,
			      std::vector<int>);
  void calcTrackPID(EvtRecTrackIterator,
		    double& ,
		    double&);
  bool hasGoodPiPiVertex(RecMdcKalTrack *,
			 RecMdcKalTrack *);
  int selectNeutralTracks(SmartDataPtr<EvtRecEvent>,
			  SmartDataPtr<EvtRecTrackCol>);
  bool passVertexSelection(CLHEP::Hep3Vector,
			   RecMdcKalTrack* ); 
  CLHEP::Hep3Vector getOrigin();
}; 


//
// module declare
//

DECLARE_ALGORITHM_FACTORY( Jpsi2invi )
DECLARE_FACTORY_ENTRIES( Jpsi2invi ) {
  DECLARE_ALGORITHM(Jpsi2invi);
}

LOAD_FACTORY_ENTRIES( Jpsi2invi )

//
// constants
//

const double PION_MASS = 0.139570;


//
// member functions
//
  
Jpsi2invi::Jpsi2invi(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator) {
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
}


StatusCode Jpsi2invi::initialize(){
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << ">>>>>>> in initialize()" << endmsg;

  m_fout = new TFile(m_output_filename.c_str(), "RECREATE");
  m_fout->cd(); 


  // SmartDataPtr<IHistogram1D> h1(histoSvc(), "hevtflw");
  // if( h1 ) h_evtflw = h1;
  // else h_evtflw = histoSvc()->book( "hevtflw", "eventflow", 10, 0, 10);

  book_histogram(); 

  if (! book_ntuple_signal(log) ) return StatusCode::FAILURE;

  //--------end booking --------

  log << MSG::INFO << "successfully return from initialize()" <<endmsg;
  return StatusCode::SUCCESS;
}


StatusCode Jpsi2invi::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;

  h_evtflw->Fill(0); // raw 
  SmartDataPtr<Event::EventHeader>eventHeader(eventSvc(),"/Event/EventHeader");
  if(!eventHeader) return StatusCode::FAILURE;

  m_run = eventHeader->runNumber();
  m_event = eventHeader->eventNumber();

  buildJpsiToInvisible();
  // m_tuple1->write();

  return StatusCode::SUCCESS; 

}


StatusCode Jpsi2invi::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in finalize()" << endmsg;

  m_fout->cd();
  m_tuple1->write();
  h_evtflw->Write();
  m_fout->Close();
  
  return StatusCode::SUCCESS;
}


Jpsi2invi::~Jpsi2invi() {
}


void Jpsi2invi::book_histogram() {

  h_evtflw = new TH1F("hevtflw", "eventflow", 10, 0, 10);
  h_evtflw->GetXaxis()->SetBinLabel(1,"raw");
  h_evtflw->GetXaxis()->SetBinLabel(2,"N_{Good}=2");
  h_evtflw->GetXaxis()->SetBinLabel(3,"N_{#gamma}=0");
  h_evtflw->GetXaxis()->SetBinLabel(4,"|cos#theta|<0.8");
  h_evtflw->GetXaxis()->SetBinLabel(5,"|p|<0.45");
  h_evtflw->GetXaxis()->SetBinLabel(6,"b-veto"); 
  h_evtflw->GetXaxis()->SetBinLabel(7,"#Delta #phi(jet,E_{T}^{miss})>0.5");
  h_evtflw->GetXaxis()->SetBinLabel(8,"E_{T}^{miss}>80");
  h_evtflw->GetXaxis()->SetBinLabel(9,"E_{T}^{miss}>125");
}

bool Jpsi2invi::book_ntuple_signal(MsgStream log) {
  StatusCode status;

  NTuplePtr nt1(ntupleSvc(), "FILE1/signal");
  if ( nt1 ) m_tuple1 = nt1;
  else {
    m_tuple1 = ntupleSvc()->book ("FILE1/signal", CLID_ColumnWiseTuple, "signal");
    if ( m_tuple1 ) {

      // common info
      status = m_tuple1->addItem ("run",   m_run);
      status = m_tuple1->addItem ("event", m_event);

      // charged tracks  
      status = m_tuple1->addItem ("ncharged", m_ncharged);
      status = m_tuple1->addItem ("nptrk", m_nptrk);
      status = m_tuple1->addItem ("nmtrk", m_nmtrk);

      // vertex
      status = m_tuple1->addItem ("vr0", m_vr0 );
      status = m_tuple1->addItem ("vz0", m_vz0 );
      status = m_tuple1->addItem ("vtx_mrecpipi", m_vtx_mrecpipi);

      // neutral tracks
      status = m_tuple1->addItem ("nshow",  m_nshow);
      status = m_tuple1->addItem ("ngam",  m_ngam);
      
    } else { 
      log << MSG::ERROR << "    Cannot book N-tuple:" 
	  << long(m_tuple1) << endmsg;
      return false;
    }
  }
  return true; 
}


void Jpsi2invi::buildJpsiToInvisible() {

  SmartDataPtr<EvtRecEvent>evtRecEvent(eventSvc(),"/Event/EvtRec/EvtRecEvent");
  if(!evtRecEvent) return;

  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  if(!evtRecTrkCol) return;

  std::vector<int> iPGood, iMGood; 
  selectChargedTracks(evtRecEvent, evtRecTrkCol, iPGood, iMGood);

  if (m_ncharged != 2) return;
  h_evtflw->Fill(1); // N_{Good} = 2 
  
  selectNeutralTracks(evtRecEvent, evtRecTrkCol);
  if (m_ngam != 0) return;
  h_evtflw->Fill(2); // N_{#gamma} = 0 
    
  selectPionPlusPionMinus(evtRecTrkCol, iPGood, iMGood); 
}


CLHEP::Hep3Vector Jpsi2invi::getOrigin() {
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


bool Jpsi2invi::passVertexSelection(CLHEP::Hep3Vector xorigin,
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


int Jpsi2invi::selectChargedTracks(SmartDataPtr<EvtRecEvent> evtRecEvent,
				   SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
				   std::vector<int> & iPGood,
				   std::vector<int> & iMGood) {

  CLHEP::Hep3Vector xorigin = getOrigin(); 

  // Good Kalman Track
  std::vector<int> iGood;
  iGood.clear();
  // std::vector<int> iPGood;
  iPGood.clear();
  //std::vector<int> iMGood;
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

int Jpsi2invi::selectPionPlusPionMinus(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
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
	//evtflw_filled = true; 
      }

      // pion momentum
      
      // cout << "<<<< mdcTrk_p: " << fabs(mdcTrk_p->p())
      // 	   << "< m pion_momtum max " << m_pion_momentum_max
      // 	   << "  , " << fabs(mdcTrk_m->p()) << endl; 
      
      if ( ! ( fabs(mdcTrk_p->p()) < m_pion_momentum_max  &&
      	       fabs(mdcTrk_m->p()) < m_pion_momentum_max )) continue;

      
      if ( !evtflw_filled ) {
	h_evtflw->Fill(4); //|p|<0.45 
	//evtflw_filled = true; 
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

      // apply vertex fit
      RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin()+iPGood[i1]))->mdcKalTrack();
      RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin()+iMGood[i2]))->mdcKalTrack();
	    
      if (! hasGoodPiPiVertex(pipTrk, pimTrk) ) continue; 

      npipi++;
      evtflw_filled = true;
    }
  } 

  return npipi; 
}

void Jpsi2invi::calcTrackPID(EvtRecTrackIterator itTrk_p,
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
    //prob_p  = pidp->probProton();
  }
}

bool Jpsi2invi::hasGoodPiPiVertex(RecMdcKalTrack *pipTrk,
				  RecMdcKalTrack *pimTrk) {

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

  if( ! ( p4_vtx_recpipi.m() >= m_dipion_mass_min &&
	  p4_vtx_recpipi.m() <= m_dipion_mass_max) ) return false;

  if( ! (cospipi < m_pipi_costheta_max) ) return false;

  if( ! (fabs(cos2pisys) < m_pipisys_costheta_max ) ) return false;

  m_vtx_mrecpipi = p4_vtx_recpipi.m();
  
  return true;
}





int Jpsi2invi::selectNeutralTracks(SmartDataPtr<EvtRecEvent> evtRecEvent,
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
