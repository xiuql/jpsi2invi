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

  // Define Ntuples
  
  // general info 
  NTuple::Tuple* m_tuple8;
  NTuple::Item<long> m_run;
  NTuple::Item<long> m_event;
  NTuple::Item<double> m_vr0;
  NTuple::Item<double> m_vz0;
  
  // functions
  bool getGeneralInfo();
  bool passPreSelection();
  bool passVertexSelection(CLHEP::Hep3Vector, RecMdcKalTrack* ); 
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
  log << MSG::INFO << ">>>>>>> in initialize()" << endmsg;

  StatusCode status;

  NTuplePtr nt8(ntupleSvc(), "FILE1/infmom");
  if ( nt8 ) m_tuple8 = nt8;
  else {
    m_tuple8 = ntupleSvc()->book ("FILE1/infmom", CLID_ColumnWiseTuple, 
				  "information with momentum method");
    if ( m_tuple8 )    {
    status = m_tuple8->addItem ("run", m_run );
    status = m_tuple8->addItem ("event", m_event );
    status = m_tuple8->addItem ("vr0", m_vr0 );
    status = m_tuple8->addItem ("vz0", m_vz0 );
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
  if (!passPreSelection()) return StatusCode::FAILURE; 
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
  SmartDataPtr<EvtRecEvent>evtRecEvent(eventSvc(),"/Event/EvtRec/EvtRecEvent");
  if(!evtRecEvent) return false; 

  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  if(!evtRecTrkCol) return false;

  CLHEP::Hep3Vector xorigin = getOrigin(); 

  // loop through charged tracks 
  for(int i = 0; i < evtRecEvent->totalCharged(); i++){
    // get mdcTrk 
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    if(! ((*itTrk)->isMdcKalTrackValid()) ) continue;
    RecMdcKalTrack* mdcTrk = (*itTrk)->mdcKalTrack();
    if(mdcTrk->p()<m_distin_pionlep)
      mdcTrk->setPidType(RecMdcKalTrack::pion);
    else
      mdcTrk->setPidType(RecMdcKalTrack::muon);

    if (!passVertexSelection(xorigin, mdcTrk) ) continue; 

  }
  
  return true; 
}

CLHEP::Hep3Vector Jpsi2invi::getOrigin() {
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

