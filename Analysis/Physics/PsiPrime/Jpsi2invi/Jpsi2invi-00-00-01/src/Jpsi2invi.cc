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
#include "EventModel/EventHeader.h"


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

  NTuple::Tuple* m_tuple8; 
  NTuple::Item<long> m_run;
  NTuple::Item<long> m_event;

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
  
  //std::cout << "execute()" << std::endl;


  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;

  StatusCode sc=StatusCode::SUCCESS;
  //save the events passed selection to a new file
  setFilterPassed(false);

  SmartDataPtr<Event::EventHeader>eventHeader(eventSvc(),"/Event/EventHeader");
  if(!eventHeader){
    log << MSG::ERROR << "EventHeader not found" << endreq;
    return sc;
  }
  m_run = eventHeader->runNumber();
  m_event = eventHeader->eventNumber();

  m_tuple8->write();
  
}


StatusCode Jpsi2invi::finalize() {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in finalize()" << endmsg;
  return StatusCode::SUCCESS;
}


Jpsi2invi::~Jpsi2invi() {
}
