// -*- C++ -*-
//
//
// Description: J/psi -> Invisible 
//
// Original Author:  SHI Xin <shixin@ihep.ac.cn>
//         Created:  [2016-03-23 Wed 09:12] 
//         Inspired by Zhu's psi' --> J/psi pi pi example code 
//


//
// system include files
//


//
// user include files
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

#include "AIDA/IHistogram1D.h"

#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/Helix.h"

//
// class declaration
//

class Jpsi2invi : public Algorithm {
  
public:
  Jpsi2invi(const std::string& name, ISvcLocator* pSvcLocator);
  ~Jpsi2invi(); 
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private:
  // Declare r0, z0 cut for charged tracks
  double m_vr0cut, m_vz0cut;
  double m_distin_pionlep;

  // Define Histograms
  IHistogram1D *h_vr0;
  IHistogram1D *h_vz0;
  
  // Define Ntuples
  
  // general info 
  NTuple::Tuple* m_tuple8;
  NTuple::Item<long> m_run;
  NTuple::Item<long> m_event;
  NTuple::Item<double> m_vr0;
  // NTuple::Item<double> m_vz0;
  
  // functions
  bool getGeneralInfo();
  bool passPreSelection();
  bool passVertexSelection(); 
  
  
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


  
//
// member functions
//
  
Jpsi2invi::Jpsi2invi(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator) {
  //Declare the properties  
  declareProperty("Vr0cut", m_vr0cut=1.0);
  declareProperty("Vz0cut", m_vz0cut=10.0);
  declareProperty("DiffPionLep", m_distin_pionlep=0.8);

}


StatusCode Jpsi2invi::initialize(){
  MsgStream log(msgSvc(), name());
  std::cout << ">>>>>>>> here in initialize! " << std::endl;
  log << MSG::INFO << ">>>>>>> in initialize()" << endmsg;

  StatusCode status;

  SmartDataPtr<IHistogram1D> h1(histoSvc(),"hvr");
  if( h1 ) h_vr0 = h1;
  else h_vr0 = histoSvc()->book( "hvr" ,  "vr0", 200, -40., 40.);

  SmartDataPtr<IHistogram1D> h2(histoSvc(),"hvz");
  if( h2 ) h_vz0 = h2;
  else h_vz0 = histoSvc()->book( "hvz" ,  "vz0", 200, -40., 40.);

  NTuplePtr nt8(ntupleSvc(), "FILE1/infmom");
  if ( nt8 ) m_tuple8 = nt8;
  else {
    m_tuple8 = ntupleSvc()->book ("FILE1/infmom", CLID_ColumnWiseTuple, 
				  "information with momentum method");
    if ( m_tuple8 )    {
    status = m_tuple8->addItem ("run", m_run );
    status = m_tuple8->addItem ("event", m_event );
    status = m_tuple8->addItem ("vr0", m_vr0 );
    }
    
    else    { 
      log << MSG::ERROR << "    Cannot book N-tuple:" 
	  << long(m_tuple8) << endmsg;
      return StatusCode::FAILURE;
    }
  }

  //--------end booking --------

  log << MSG::INFO << "successfully return from initialize()" <<endmsg;
  return StatusCode::SUCCESS;

}


StatusCode Jpsi2invi::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;

  if (!getGeneralInfo()) return StatusCode::FAILURE; 
  // if (!passPreSelection()) return StatusCode::FAILURE; 
  if (! passVertexSelection())  return StatusCode::FAILURE; 
  m_tuple8->write();
  

}


StatusCode Jpsi2invi::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in finalize()" << endmsg;
  return StatusCode::SUCCESS;
}


Jpsi2invi::~Jpsi2invi() {
}


bool Jpsi2invi::getGeneralInfo() {
  MsgStream log(msgSvc(), name());
  SmartDataPtr<Event::EventHeader>eventHeader(eventSvc(),"/Event/EventHeader");
  if(!eventHeader){
    log << MSG::ERROR << "EventHeader not found" << endreq;
    return false;
  }
  m_run = eventHeader->runNumber();
  m_event = eventHeader->eventNumber();
  return true;
}

bool Jpsi2invi::passPreSelection() {

  if (! passVertexSelection()) return false; 
  return true; 
}

bool Jpsi2invi::passVertexSelection() {
  // check x0, y0, z0, r0
  // suggest cut: |z0|<10 && r0<1 (cm)
  CLHEP::Hep3Vector xorigin(0,0,0);
  IVertexDbSvc*  vtxsvc;
  double *dbv, *vv;
  Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
  if(vtxsvc->isVertexValid()){
    dbv = vtxsvc->PrimaryVertex(); 
    vv = vtxsvc->SigmaPrimaryVertex();  
    xorigin.setX(dbv[0]);
    xorigin.setY(dbv[1]);
    xorigin.setZ(dbv[2]);
  }

  SmartDataPtr<EvtRecEvent>evtRecEvent(eventSvc(),"/Event/EvtRec/EvtRecEvent");
  if(!evtRecEvent) return false; 

  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  if(!evtRecTrkCol) return false;

  // double m_vr0(100), m_vz0(100);
  double m_vz0(100);
  
  EvtRecTrackIterator m_itTrk_begin = evtRecTrkCol->begin();
  for(int i = 0; i < evtRecEvent->totalCharged(); i++){
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    if(! ((*itTrk)->isMdcKalTrackValid()) ) continue;
    RecMdcKalTrack* mdcTrk = (*itTrk)->mdcKalTrack();
    if(mdcTrk->p()<m_distin_pionlep) mdcTrk->setPidType(RecMdcKalTrack::pion);
    else mdcTrk->setPidType(RecMdcKalTrack::muon);

    HepVector a = mdcTrk->helix();
    HepSymMatrix Ea = mdcTrk->err();
    HepPoint3D point0(0.,0.,0.);
    VFHelix helixip(point0,a,Ea);
    HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
    helixip.pivot(IP);
    HepVector vecipa = helixip.a();

    m_vz0 = vecipa[3];
    m_vr0 = vecipa[0];
    
    h_vz0->fill(m_vz0);
    h_vr0->fill(m_vr0);

    if(fabs(m_vz0) >= m_vz0cut) continue;
    if(fabs(m_vr0) >= m_vr0cut) continue;

  }
  
  return true; 
}
