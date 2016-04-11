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

#include "CLHEP/Vector/ThreeVector.h"
#include "VertexFit/IVertexDbSvc.h"

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
  double m_vr0cut; 

  // Define Histograms
  IHistogram1D *h_vr0;
  IHistogram1D *h_vz0;
  
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
}


StatusCode Jpsi2invi::initialize(){
  MsgStream log(msgSvc(), name());
  std::cout << ">>>>>>>> here in initialize! " << std::endl;
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
  m_tuple8->write();
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

  return true; 
}
