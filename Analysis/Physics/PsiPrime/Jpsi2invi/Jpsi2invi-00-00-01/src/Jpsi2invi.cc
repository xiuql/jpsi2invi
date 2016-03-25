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
  
}; 


//
// module declare
//

DECLARE_ALGORITHM_FACTORY( Jpsi2invi )
DECLARE_FACTORY_ENTRIES( Jpsi2invi ) {
  DECLARE_ALGORITHM(Jpsi2invi);
}

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
}


StatusCode Jpsi2invi::execute() {
  
  //std::cout << "execute()" << std::endl;


  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;

}


StatusCode Jpsi2invi::finalize() {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in finalize()" << endmsg;
  return StatusCode::SUCCESS;
}


Jpsi2invi::~Jpsi2invi() {
}
