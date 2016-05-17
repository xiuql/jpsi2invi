//psi'--> J/psi pion pion, J/psi --> di-leptons
//Kai Zhu (zhuk@ihep.ac.cn)
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include <iostream>
#include <stdio.h>

#include "EventModel/EventModel.h"
#include "EventModel/EventHeader.h"
#include "EventModel/Event.h"
#include "TrigEvent/TrigEvent.h"
#include "TrigEvent/TrigData.h"
#include "McTruth/McParticle.h"
#include "EvTimeEvent/RecEsTime.h"


#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"


#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/Helix.h"   
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/WTrackParameter.h"
#include "ParticleID/ParticleID.h" 
#include "MdcRecEvent/RecMdcKalTrack.h"
#include "MucRecEvent/RecMucTrack.h"


#include "PipiJpsiAlg/PipiJpsi.h"
#include "PipiJpsiAlg/PionZeroList.h"


#include <vector>
//const double twopi = 6.2831853;
using namespace std;

const double me  = 0.000511;
const double mpi = 0.13957;
const double mproton = 0.938272;
const double mmu = 0.105658;
const double mpsip = 3.686;
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
const double velc = 29.9792458;  // tof_path unit in cm.
const double PI = 3.1415926;
// const double velc = 299.792458;   // tof path unit in mm
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;


// counter for efficiency
static long m_cout_all(0), m_cout_col(0), m_cout_charge(0); 
static long m_cout_nGood(0), m_cout_mom(0), m_cout_recoil(0), m_cout_everat(0);
/////////////////////////////////////////////////////////////////////////////


PipiJpsi::PipiJpsi(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator) {
  //Declare the properties  
  
  declareProperty("OutputFileName",  m_OutputFileName  );
  declareProperty("Vr0cut", m_vr0cut=1.0);
  declareProperty("Vz0cut", m_vz0cut=10.0);
  declareProperty("PtCut", m_ptcut=0.0);
  declareProperty("MomCut", m_momcut=2.2);
  declareProperty("ChaCosthetaCut", m_cha_costheta_cut=0.93);
  declareProperty("BarEneCut", m_bar_energy_cut=0.025);
  declareProperty("EndEneCut", m_end_energy_cut=0.050);
  declareProperty("DangCut", m_dang_cut=10);
  declareProperty("MinEstCut", m_min_emctime=-0.01);
  declareProperty("MaxEstCut", m_max_emctime=14.01);
  declareProperty("BarCosthetaCut", m_bar_costheta_cut=0.8);
  declareProperty("MinEndCosthetaCut", m_min_end_costheta_cut=0.84);
  declareProperty("MaxEndCosthetaCut", m_max_end_costheta_cut=0.92);
  declareProperty("PipiDangCut",m_pipi_dang_cut=0.98);


  declareProperty("CmsEne",m_cms_ene=3.686);
  declareProperty("DiffPionLep", m_distin_pionlep=0.8);
  declareProperty("CheckDedx", m_checkDedx = false);
  declareProperty("CheckTof",  m_checkTof = false);


  declareProperty("Subsample", m_subsample_flag=false); 
  declareProperty("Trigger", m_trigger_flag=false); 
  declareProperty("Estime", m_est_flag=false);
  declareProperty("DistinEMuon", m_distin_emuon=0.5); // using e/p


  declareProperty("EventRate", m_eventrate=false);
  declareProperty("ChanDet", m_chan_det=1);


  declareProperty("Track", m_track_flag=false); // see head file for more
  declareProperty("Valid", m_valid_flag=true); // whether check kalmanfilter vaild or not
  declareProperty("Fit", m_fit_flag=true); // whether do kinematic fit
  declareProperty("Vtdb", m_vtdb_flag=true); // whether using vertex data base
  declareProperty("Debug", m_debug_flag=false); // whether do debug
  declareProperty("NonFsr", m_non_fsr=false); // whether do debug
  declareProperty("OnlyFour", m_only_four=false); // only save 4 tracks
  declareProperty("Premc", m_premc=false); // whether save information before selection
  declareProperty("Kaon", m_kaon_flag=false); // whether analyze K+ K- J/psi
  
  declareProperty("saveTopo", m_saveTopo = 1);//need to be re-evaluated when running different samples(1 for MC)
  declareProperty("saveMCTruth", m_saveMCTruth = 0);//need to be re-evaluated when running different samples(only 1 for exclusiveMC)

  declareProperty("saveTopoTree", m_saveTopoTree = 0);//need to be re-evaluated when running different samples(1 for MC)
  declareProperty("saveNbInfo", m_saveNbInfo = 1);//need to be re-evaluated when running different samples(1 for MC)
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode PipiJpsi::initialize(){
  MsgStream log(msgSvc(), name());
	
cout<<"This initialize for chenq!"<<endl;
  log << MSG::INFO << "in initialize()" << endmsg;
  
  StatusCode status;
  
      TString s_OutputFileName(m_OutputFileName);

  //Tree and File Modified
	saveFile = new TFile(s_OutputFileName, "recreate");
//	TreeAna = new TTree("save", "save");
	TopoTree = new TTree("topoall","topoall");
	NbInfo = new TTree("nbinfo","nbinfo");
	TreeAna = new TTree("save","save");
	//Overall info
	      TreeAna->Branch("runid", &runid, "runid/I");
	      TreeAna->Branch("evtid", &evtid, "evtid/I");

	if(m_saveNbInfo== 1)
	{
		NbInfo->Branch("m_cout_all", &m_cout_all, "m_cout_all/L");
		NbInfo->Branch("m_cout_col", &m_cout_col, "m_cout_col/L");
		NbInfo->Branch("m_cout_charge", &m_cout_charge, "m_cout_charge/L");
		NbInfo->Branch("m_cout_nGood", &m_cout_nGood, "m_cout_nGood/L");
		NbInfo->Branch("m_cout_mom", &m_cout_mom, "m_cout_mom/L");
		NbInfo->Branch("m_cout_recoil", &m_cout_recoil, "m_cout_recoil/L");
		NbInfo->Branch("m_cout_everat", &m_cout_everat, "m_cout_everat/L");
	}
	
	if(m_saveTopoTree == 1)
	{
		TopoTree->Branch("run", &m_run, "run/L");
		//TopoTree->Branch("rec", &m_rec, "rec/L");
		TopoTree->Branch("indexmc", &m_idxmc, "indexmc/L");

		//TopoTree->Branch("trkidx", m_trkidx, "trkidx[indexmc]/I");
		TopoTree->Branch("pdgid", m_pdgid, "pdgid[indexmc]/L");
		//TopoTree->Branch("drank", m_drank, "drank[indexmc]/I");
		TopoTree->Branch("motheridx", m_motheridx, "motheridx[indexmc]/I");
		//TopoTree->Branch("motherpid", m_motherpid, "motherpid[indexmc]/I");

	}
	
      TreeAna->Branch("f_cos_theta", &f_cos_theta, "f_cos_theta/D");
      TreeAna->Branch("f_how_near", &f_how_near, "f_how_near/D");
      TreeAna->Branch("f_miss_mass", &f_miss_mass, "f_miss_mass/D");
      TreeAna->Branch("f_phi", &f_phi, "f_phi/D");
      TreeAna->Branch("m_cms_lepm", &m_cms_lepm, "m_cms_lepm/D");
      TreeAna->Branch("m_cms_lepp", &m_cms_lepp, "m_cms_lepp/D");
      TreeAna->Branch("m_cos_theta", &m_cos_theta, "m_cos_theta/D");
      TreeAna->Branch("m_dedx_e", &m_dedx_e, "m_dedx_e/D");
      TreeAna->Branch("m_dedx_ghits", &m_dedx_ghits, "m_dedx_ghits/D");
      TreeAna->Branch("m_dedx_mu", &m_dedx_mu, "m_dedx_mu/D");
      TreeAna->Branch("m_dedx_pi", &m_dedx_pi, "m_dedx_pi/D");
      TreeAna->Branch("m_ep_ratio", &m_ep_ratio, "m_ep_ratio/D");
      TreeAna->Branch("m_how_near", &m_how_near, "m_how_near/D");
      TreeAna->Branch("m_miss_mass", &m_miss_mass, "m_miss_mass/D");
      TreeAna->Branch("m_motheridx", &m_motheridx, "m_motheridx/D");
      TreeAna->Branch("m_muc_dep", &m_muc_dep, "m_muc_dep/D");
      TreeAna->Branch("m_pdgid", &m_pdgid, "m_pdgid/D");
      TreeAna->Branch("m_phi", &m_phi, "m_phi/D");
      TreeAna->Branch("m_prob_e", &m_prob_e, "m_prob_e/D");
      TreeAna->Branch("m_prob_mu", &m_prob_mu, "m_prob_mu/D");
      TreeAna->Branch("m_prob_pi", &m_prob_pi, "m_prob_pi/D");
      TreeAna->Branch("m_tof_e", &m_tof_e, "m_tof_e/D");
      TreeAna->Branch("m_tof_mu", &m_tof_mu, "m_tof_mu/D");
      TreeAna->Branch("m_tof_pi", &m_tof_pi, "m_tof_pi/D");
      TreeAna->Branch("m_trig_chan", &m_trig_chan, "m_trig_chan/D");
      TreeAna->Branch("m_trig_cond", &m_trig_cond, "m_trig_cond/D");
      TreeAna->Branch("m_true_cms_lepm", &m_true_cms_lepm, "m_true_cms_lepm/D");
      TreeAna->Branch("m_true_cms_lepp", &m_true_cms_lepp, "m_true_cms_lepp/D");
      TreeAna->Branch("m_dthe", &m_dthe, "m_dthe/D");
      TreeAna->Branch("m_dphi", &m_dphi, "m_dphi/D");
      TreeAna->Branch("m_dang", &m_dang, "m_dang/D");
      TreeAna->Branch("m_eraw", &m_eraw, "m_eraw/D");
      TreeAna->Branch("m_nGam", &m_nGam, "m_nGam/L");
      TreeAna->Branch("m_ptrk", &m_ptrk, "m_ptrk/D");
      TreeAna->Branch("m_chie", &m_chie, "m_chie/D");
      TreeAna->Branch("m_chimu", &m_chimu, "m_chimu/D");
      TreeAna->Branch("m_chipi", &m_chipi, "m_chipi/D");
      TreeAna->Branch("m_chik", &m_chik, "m_chik/D");
      TreeAna->Branch("m_chip", &m_chip, "m_chip/D");
      TreeAna->Branch("m_probPH", &m_probPH, "m_probPH/D");
      TreeAna->Branch("m_normPH", &m_normPH, "m_normPH/D");
      TreeAna->Branch("m_ghit", &m_ghit, "m_ghit/D");
      TreeAna->Branch("m_thit", &m_thit, "m_thit/D");
      TreeAna->Branch("m_ptot_etof", &m_ptot_etof, "m_ptot_etof/D");
      TreeAna->Branch("m_cntr_etof", &m_cntr_etof, "m_cntr_etof/D");
      TreeAna->Branch("m_path_etof", &m_path_etof, "m_path_etof/D");
      TreeAna->Branch("m_tof_etof", &m_tof_etof, "m_tof_etof/D");
      TreeAna->Branch("m_te_etof", &m_te_etof, "m_te_etof/D");
      TreeAna->Branch("m_tmu_etof", &m_tmu_etof, "m_tmu_etof/D");
      TreeAna->Branch("m_tpi_etof", &m_tpi_etof, "m_tpi_etof/D");
      TreeAna->Branch("m_tk_etof", &m_tk_etof, "m_tk_etof/D");
      TreeAna->Branch("m_tp_etof", &m_tp_etof, "m_tp_etof/D");
      TreeAna->Branch("m_ph_etof", &m_ph_etof, "m_ph_etof/D");
      TreeAna->Branch("m_rhit_etof", &m_rhit_etof, "m_rhit_etof/D");
      TreeAna->Branch("m_qual_etof", &m_qual_etof, "m_qual_etof/D");
      TreeAna->Branch("m_ptot_btof1", &m_ptot_btof1, "m_ptot_btof1/D");
      TreeAna->Branch("m_cntr_btof1", &m_cntr_btof1, "m_cntr_btof1/D");
      TreeAna->Branch("m_path_btof1", &m_path_btof1, "m_path_btof1/D");
      TreeAna->Branch("m_tof_btof1", &m_tof_btof1, "m_tof_btof1/D");
      TreeAna->Branch("m_te_btof1", &m_te_btof1, "m_te_btof1/D");
      TreeAna->Branch("m_tmu_btof1", &m_tmu_btof1, "m_tmu_btof1/D");
      TreeAna->Branch("m_tpi_btof1", &m_tpi_btof1, "m_tpi_btof1/D");
      TreeAna->Branch("m_tk_btof1", &m_tk_btof1, "m_tk_btof1/D");
      TreeAna->Branch("m_tp_btof1", &m_tp_btof1, "m_tp_btof1/D");
      TreeAna->Branch("m_ph_btof1", &m_ph_btof1, "m_ph_btof1/D");
      TreeAna->Branch("m_zhit_btof1", &m_zhit_btof1, "m_zhit_btof1/D");
      TreeAna->Branch("m_qual_btof1", &m_qual_btof1, "m_qual_btof1/D");
      TreeAna->Branch("m_ptot_btof2", &m_ptot_btof2, "m_ptot_btof2/D");
      TreeAna->Branch("m_cntr_btof2", &m_cntr_btof2, "m_cntr_btof2/D");
      TreeAna->Branch("m_path_btof2", &m_path_btof2, "m_path_btof2/D");
      TreeAna->Branch("m_tof_btof2", &m_tof_btof2, "m_tof_btof2/D");
      TreeAna->Branch("m_te_btof2", &m_te_btof2, "m_te_btof2/D");
      TreeAna->Branch("m_tmu_btof2", &m_tmu_btof2, "m_tmu_btof2/D");
      TreeAna->Branch("m_tpi_btof2", &m_tpi_btof2, "m_tpi_btof2/D");
      TreeAna->Branch("m_tk_btof2", &m_tk_btof2, "m_tk_btof2/D");
      TreeAna->Branch("m_tp_btof2", &m_tp_btof2, "m_tp_btof2/D");
      TreeAna->Branch("m_ph_btof2", &m_ph_btof2, "m_ph_btof2/D");
      TreeAna->Branch("m_zhit_btof2", &m_zhit_btof2, "m_zhit_btof2/D");
      TreeAna->Branch("m_qual_btof2", &m_qual_btof2, "m_qual_btof2/D");
      TreeAna->Branch("m_pionp_bef", &m_pionp_bef, "m_pionp_bef/D");
      TreeAna->Branch("m_pionm_bef", &m_pionm_bef, "m_pionm_bef/D");
      TreeAna->Branch("m_pionp_pt_bef", &m_pionp_pt_bef, "m_pionp_pt_bef/D");
      TreeAna->Branch("m_pionm_pt_bef", &m_pionm_pt_bef, "m_pionm_pt_bef/D");
      TreeAna->Branch("m_mom_lepm", &m_mom_lepm, "m_mom_lepm/D");
      TreeAna->Branch("m_mom_lepp", &m_mom_lepp, "m_mom_lepp/D");
      TreeAna->Branch("m_mom_pionp", &m_mom_pionp, "m_mom_pionp/D");
      TreeAna->Branch("m_mom_pionm", &m_mom_pionm, "m_mom_pionm/D");
      TreeAna->Branch("m_pipi_dang", &m_pipi_dang, "m_pipi_dang/D");
      TreeAna->Branch("m_cms_index", &m_cms_index, "m_cms_index/D");
      //TreeAna->Branch("m_cms_lepp", &m_cms_lepp, "m_cms_lepp[indexmc]/D");
      //TreeAna->Branch("m_cms_lepm", &m_cms_lepm, "m_cms_lepm[indexmc]/D");
      TreeAna->Branch("m_mass_twopi", &m_mass_twopi, "m_mass_twopi/D");
      TreeAna->Branch("m_mass_jpsi", &m_mass_jpsi, "m_mass_jpsi/D");
      TreeAna->Branch("m_mass_recoil", &m_mass_recoil, "m_mass_recoil/D");
      TreeAna->Branch("m_inv_mass", &m_inv_mass, "m_inv_mass/D");
      TreeAna->Branch("m_tot_e", &m_tot_e, "m_tot_e/D");
      TreeAna->Branch("m_tot_px", &m_tot_px, "m_tot_px/D");
      TreeAna->Branch("m_tot_py", &m_tot_py, "m_tot_py/D");
      TreeAna->Branch("m_tot_pz", &m_tot_pz, "m_tot_pz/D");
      TreeAna->Branch("m_emc_ene", &m_emc_ene, "m_emc_ene/D");
      TreeAna->Branch("m_event_flag", &m_event_flag, "m_event_flag/D");
      TreeAna->Branch("m_pt_lepm", &m_pt_lepm, "m_pt_lepm/D");
      TreeAna->Branch("m_pt_lepp", &m_pt_lepp, "m_pt_lepp/D");
      TreeAna->Branch("m_pt_pionp", &m_pt_pionp, "m_pt_pionp/D");
      TreeAna->Branch("m_pt_pionm", &m_pt_pionm, "m_pt_pionm/D");
      TreeAna->Branch("m_run", &m_run, "m_run/D");
      TreeAna->Branch("m_event", &m_event, "m_event/D");
      TreeAna->Branch("m_index", &m_index, "m_index/D");
      //TreeAna->Branch("m_cos_theta", &m_cos_theta, "m_cos_theta[indexmc]/D");
      //TreeAna->Branch("m_ep_ratio", &m_ep_ratio, "m_ep_ratio[indexmc]/D");
      //TreeAna->Branch("m_phi", &m_phi, "m_phi[indexmc]/D");
      //TreeAna->Branch("m_prob_e", &m_prob_e, "m_prob_e[indexmc]/D");
      //TreeAna->Branch("m_prob_mu", &m_prob_mu, "m_prob_mu[indexmc]/D");
      //TreeAna->Branch("m_prob_pi", &m_prob_pi, "m_prob_pi[indexmc]/D");
      //TreeAna->Branch("m_dedx_e", &m_dedx_e, "m_dedx_e[indexmc]/D");
      //TreeAna->Branch("m_dedx_mu", &m_dedx_mu, "m_dedx_mu[indexmc]/D");
      //TreeAna->Branch("m_dedx_pi", &m_dedx_pi, "m_dedx_pi[indexmc]/D");
      //TreeAna->Branch("m_tof_e", &m_tof_e, "m_tof_e[indexmc]/D");
      //TreeAna->Branch("m_tof_mu", &m_tof_mu, "m_tof_mu[indexmc]/D");
      //TreeAna->Branch("m_tof_pi", &m_tof_pi, "m_tof_pi[indexmc]/D");
      //TreeAna->Branch("m_dedx_ghits", &m_dedx_ghits, "m_dedx_ghits[indexmc]/D");
      //TreeAna->Branch("m_muc_dep", &m_muc_dep, "m_muc_dep[indexmc]/D");
      TreeAna->Branch("m_four_mom", &m_four_mom, "m_four_mom[4][5]/D");
      TreeAna->Branch("m_recoil_four_mom", &m_recoil_four_mom, "m_recoil_four_mom[4][5]/D");
      //TreeAna->Branch("m_miss_mass", &m_miss_mass, "m_miss_mass[indexmc]/D");
      //TreeAna->Branch("m_how_near", &m_how_near, "m_how_near[indexmc]/D");
      TreeAna->Branch("m_vertex", &m_vertex, "m_vertex[4][5]/D");
      TreeAna->Branch("m_num_pi0", &m_num_pi0, "m_num_pi0/D");
      TreeAna->Branch("m_like_pi0_mass", &m_like_pi0_mass, "m_like_pi0_mass/D");
      TreeAna->Branch("m_pipipi0_mass", &m_pipipi0_mass, "m_pipipi0_mass/D");
      TreeAna->Branch("m_omegajpsi_mass", &m_omegajpsi_mass, "m_omegajpsi_mass/D");
      TreeAna->Branch("m_gammajpsi_mass", &m_gammajpsi_mass, "m_gammajpsi_mass/D");
      TreeAna->Branch("f_gammajpsi_mass", &f_gammajpsi_mass, "f_gammajpsi_mass/D");
      TreeAna->Branch("f_chisq_gammajpsi", &f_chisq_gammajpsi, "f_chisq_gammajpsi/D");
      TreeAna->Branch("m_max_gam_ene", &m_max_gam_ene, "m_max_gam_ene/D");
      TreeAna->Branch("m_max_gam_costhe", &m_max_gam_costhe, "m_max_gam_costhe/D");
      TreeAna->Branch("m_max_gam_phi", &m_max_gam_phi, "m_max_gam_phi/D");
      TreeAna->Branch("m_max_gam_cms_ene", &m_max_gam_cms_ene, "m_max_gam_cms_ene/D");
      TreeAna->Branch("m_max_gam_dang", &m_max_gam_dang, "m_max_gam_dang/D");
      TreeAna->Branch("m_like_eta_mass", &m_like_eta_mass, "m_like_eta_mass/D");
      TreeAna->Branch("f_chisq_vt", &f_chisq_vt, "f_chisq_vt/D");
      TreeAna->Branch("f_chisq_jpsi", &f_chisq_jpsi, "f_chisq_jpsi/D");
      TreeAna->Branch("f_mom_lepm", &f_mom_lepm, "f_mom_lepm/D");
      TreeAna->Branch("f_mom_lepp", &f_mom_lepp, "f_mom_lepp/D");
      TreeAna->Branch("f_mom_pionp", &f_mom_pionp, "f_mom_pionp/D");
      TreeAna->Branch("f_mom_pionm", &f_mom_pionm, "f_mom_pionm/D");
      TreeAna->Branch("f_cms_lepp", &f_cms_lepp, "f_cms_lepp/D");
      TreeAna->Branch("f_cms_lepm", &f_cms_lepm, "f_cms_lepm/D");
      TreeAna->Branch("f_pipi_dang", &f_pipi_dang, "f_pipi_dang/D");
      TreeAna->Branch("f_mass_twopi", &f_mass_twopi, "f_mass_twopi/D");
      TreeAna->Branch("f_mass_jpsi", &f_mass_jpsi, "f_mass_jpsi/D");
      TreeAna->Branch("f_mass_recoil", &f_mass_recoil, "f_mass_recoil/D");
      TreeAna->Branch("f_inv_mass", &f_inv_mass, "f_inv_mass/D");
      TreeAna->Branch("f_tot_e", &f_tot_e, "f_tot_e/D");
      TreeAna->Branch("f_tot_px", &f_tot_px, "f_tot_px/D");
      TreeAna->Branch("f_tot_py", &f_tot_py, "f_tot_py/D");
      TreeAna->Branch("f_tot_pz", &f_tot_pz, "f_tot_pz/D");
      TreeAna->Branch("f_pt_lepm", &f_pt_lepm, "f_pt_lepm/D");
      TreeAna->Branch("f_pt_lepp", &f_pt_lepp, "f_pt_lepp/D");
      TreeAna->Branch("f_pt_pionp", &f_pt_pionp, "f_pt_pionp/D");
      TreeAna->Branch("f_pt_pionm", &f_pt_pionm, "f_pt_pionm/D");
      TreeAna->Branch("f_index", &f_index, "f_index/D");
      //TreeAna->Branch("f_cos_theta", &f_cos_theta, "f_cos_theta[indexmc]/D");
      //TreeAna->Branch("f_phi", &f_phi, "f_phi[indexmc]/D");
      TreeAna->Branch("f_four_mom", &f_four_mom, "f_four_mom[4][5]/D");
      TreeAna->Branch("f_recoil_four_mom", &f_recoil_four_mom, "f_recoil_four_mom[4][5]/D");
      //TreeAna->Branch("f_miss_mass", &f_miss_mass, "f_miss_mass[indexmc]/D");
      //TreeAna->Branch("f_how_near", &f_how_near, "f_how_near[indexmc]/D");
      TreeAna->Branch("m_pion_matched", &m_pion_matched, "m_pion_matched/D");
      TreeAna->Branch("m_lep_matched", &m_lep_matched, "m_lep_matched/D");
      TreeAna->Branch("m_Nch", &m_Nch, "m_Nch/D");
      TreeAna->Branch("m_Nty", &m_Nty, "m_Nty/D");
      TreeAna->Branch("m_idxmc", &m_idxmc, "m_idxmc/D");
      //TreeAna->Branch("m_pdgid", &m_pdgid, "m_pdgid[indexmc]/D");
      //TreeAna->Branch("m_motheridx", &m_motheridx, "m_motheridx[indexmc]/D");
      TreeAna->Branch("m_xyzid", &m_xyzid, "m_xyzid/D");
      TreeAna->Branch("m_true_cms_index", &m_true_cms_index, "m_true_cms_index/D");
      //TreeAna->Branch("m_true_cms_lepp", &m_true_cms_lepp, "m_true_cms_lepp[indexmc]/D");
      //TreeAna->Branch("m_true_cms_lepm", &m_true_cms_lepm, "m_true_cms_lepm[indexmc]/D");
      TreeAna->Branch("m_trig_index", &m_trig_index, "trig_index/D");
      //TreeAna->Branch("m_trig_cond", &m_trig_cond, "m_trig_cond[indexmc]/D");
      //TreeAna->Branch("m_trig_chan", &m_trig_chan, "m_trig_chan[indexmc]/D");
      TreeAna->Branch("m_trig_timewindow", &m_trig_timewindow, "m_trig_timewindow/D");
      TreeAna->Branch("m_trig_timetype", &m_trig_timetype, "m_trig_timetype/D");
      TreeAna->Branch("m_est_start", &m_est_start, "m_est_start/D");
      TreeAna->Branch("m_est_status", &m_est_status, "m_est_status/D");
      TreeAna->Branch("m_est_quality", &m_est_quality, "m_est_quality/D");
  
  
 
  //
  //--------end of book--------
  //


  log << MSG::INFO << "successfully return from initialize()" <<endmsg;
  return StatusCode::SUCCESS;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode PipiJpsi::execute() {
  
  //std::cout << "execute()" << std::endl;
cout<<"This execute  for chenq!"<<endl;

	//TopoTree->Fill();

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;
  m_cout_all ++;
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
  if(m_cout_all%10000==0) cout << "events: " << m_cout_all <<  "  run No: "
  			      << m_run << " event No+1: " << m_event+1 << endl;
  if(m_run<0){
    m_Nch = eventHeader->flag1();
    m_Nty = eventHeader->flag2();
  }


  //MC before selection
  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), 
						   "/Event/MC/McParticleCol");
  if(m_run<0 && !m_track_flag){
    if(!mcParticleCol){
      log << MSG::ERROR << "Could not retrieve McParticelCol" << endreq;
      return StatusCode::FAILURE;
    }
    else{
      bool psipDecay(false);
      Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
      for (; iter_mc != mcParticleCol->end(); iter_mc++){
        if ((*iter_mc)->primaryParticle() ) continue;
	if (!(*iter_mc)->decayFromGenerator() ) continue;
        //if ( ((*iter_mc)->mother()).trackIndex()<3 ) continue;
        if ((*iter_mc)->particleProperty()==100443) psipDecay = true;
        if (!psipDecay) continue;
	//if(!(*iter_mc)->leafParticle()) continue;
	if((*iter_mc)->particleProperty() == 211){
	  m_pionp_bef = (*iter_mc)->initialFourMomentum().vect().mag();
	  m_pionp_pt_bef = (*iter_mc)->initialFourMomentum().vect().perp();
	}
	if((*iter_mc)->particleProperty() == -211){ 
	  m_pionm_bef = (*iter_mc)->initialFourMomentum().vect().mag();
	  m_pionm_pt_bef = (*iter_mc)->initialFourMomentum().vect().perp();
	}
      }
    }
  }


  SmartDataPtr<EvtRecEvent>evtRecEvent(eventSvc(),"/Event/EvtRec/EvtRecEvent");
  if(!evtRecEvent){
    log << MSG::ERROR << "EvtRecEvent not found" << endreq;
    return sc;
  }
  log << MSG::DEBUG <<"ncharg, nneu, tottks = " 
      << evtRecEvent->totalCharged() << " , "
      << evtRecEvent->totalNeutral() << " , "
      << evtRecEvent->totalTracks() <<endreq;
    
  SmartDataPtr<EvtRecTrackCol> 
    evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  if(!evtRecTrkCol){
    log << MSG::ERROR << "EvtRecTrackCol" << endreq;
    return sc;
  }


  if(m_trigger_flag && m_run>0){
    SmartDataPtr<TrigData> trigData(eventSvc(),EventModel::Trig::TrigData);
    if (!trigData) {
      log << MSG::FATAL 
          << "Could not find Trigger Data for physics analysis" << endreq;
      return StatusCode::FAILURE;
    }
    // test event rate
    int m_trig_tot(0), m_trig_which(-1);
    if(m_eventrate){
      for(int j=0; j<16; j++){
	if(trigData->getTrigChannel(j)){
	  m_trig_tot ++;
	  m_trig_which = j+1;
	}
      }
      if(m_trig_tot==1 && m_trig_which==m_chan_det) m_cout_everat++;
      return sc;
    }
    /// Print trigger information once:
    log << MSG::DEBUG << "Trigger conditions: " << endreq;
    for(m_trig_index=0; m_trig_index < 48; m_trig_index++){
      log << MSG::DEBUG << "i:" << m_trig_index 
          << "  name:" << trigData->getTrigCondName(m_trig_index) 
          << "  cond:" << trigData->getTrigCondition(m_trig_index) << endreq;
      if(m_trig_index<16) m_trig_chan[m_trig_index] = 
			    trigData->getTrigChannel(m_trig_index);
      m_trig_cond[m_trig_index] = trigData->getTrigCondition(m_trig_index);
    }
    m_trig_timewindow = trigData->getTimeWindow();
    m_trig_timetype = trigData->getTimingType();
  }




  SmartDataPtr<RecEsTimeCol>
    recEstimeCol(eventSvc(),"/Event/Recon/RecEsTimeCol");
  if(m_est_flag & m_run>0){
    if(!recEstimeCol){
      log << MSG::ERROR << "Estime collection is not get!" << endreq;
      return StatusCode::FAILURE;
    }
    log << MSG::DEBUG <<"size of EsTime: " << recEstimeCol->size() << endreq;
    if(recEstimeCol->size() != 1) m_est_status=0;
    else{
      m_est_start = (*(recEstimeCol->begin()))->getTest();
      m_est_status = (*(recEstimeCol->begin()))->getStat();
      m_est_quality = (*(recEstimeCol->begin()))->getQuality();
    }
  }


  m_cout_col ++;
  if(evtRecEvent->totalCharged()<2 || evtRecTrkCol->size()<2 
     || evtRecEvent->totalTracks()>20 || evtRecTrkCol->size()>20)
    return sc;
  m_cout_charge ++;


  // total energy deposited in emc
  m_emc_ene = 0;
  for(int i=0; i< evtRecEvent->totalTracks(); i++){
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    if(!(*itTrk)->isEmcShowerValid()) continue;
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    m_emc_ene += emcTrk->energy();
  }


  // check x0, y0, z0, r0
  // suggest cut: |z0|<10 && r0<1 (cm)
  Hep3Vector xorigin(0,0,0);
  IVertexDbSvc*  vtxsvc;
  double *dbv, *vv;
  if(m_fit_flag & m_vtdb_flag){
    Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
    if(vtxsvc->isVertexValid()){
      dbv = vtxsvc->PrimaryVertex(); 
      vv = vtxsvc->SigmaPrimaryVertex();  
      xorigin.setX(dbv[0]);
      xorigin.setY(dbv[1]);
      xorigin.setZ(dbv[2]);
    }
  }
  // Asign four-momentum with KalmanTrack
  Vint iGood; iGood.clear();
  int igood_index[5] = {-1,-1,-1,-1,-1};
  int m_num[4]={0,0,0,0}; // number of different particles: pi+, pi-, l+, l-
  int nCharge = 0;
  m_pion_matched = 0; m_lep_matched = 0;
  HepLorentzVector m_lv_pionp, m_lv_pionm, m_lv_lepp;
  HepLorentzVector m_lv_lepm, m_lv_ele, m_lv_pos, m_lv_mum, m_lv_mup;
  HepLorentzVector m_lv_add_pi, m_lv_add_lep, m_lv_add_e, m_lv_add_u;//additional track
  double m_vr0(100), m_vz0(100);


  //  cout << "debug event: " <<   m_cout_all << "run: " << m_run << "  event: " << m_event << endl;
  //    cout << "total number of charged from evt: " << evtRecEvent->totalCharged() << endl;
  EvtRecTrackIterator m_itTrk_begin = evtRecTrkCol->begin();
  for(int i = 0; i < evtRecEvent->totalCharged(); i++){
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    if(m_valid_flag && (!(*itTrk)->isMdcKalTrackValid())) continue;
    RecMdcKalTrack* mdcTrk = (*itTrk)->mdcKalTrack(); 
    if(mdcTrk->p()<m_distin_pionlep) mdcTrk->setPidType(RecMdcKalTrack::pion);
    else mdcTrk->setPidType(RecMdcKalTrack::muon);
    if(mdcTrk->p()<m_distin_pionlep && m_kaon_flag) mdcTrk->setPidType(RecMdcKalTrack::kaon);


    HepVector a = mdcTrk->helix();
    HepSymMatrix Ea = mdcTrk->err();
    HepPoint3D point0(0.,0.,0.);
    VFHelix helixip(point0,a,Ea);
    HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
    helixip.pivot(IP);
    HepVector vecipa = helixip.a();


    m_vz0 = vecipa[3];
    m_vr0 = vecipa[0];
    //h_vz0->fill(m_vz0);
    //h_vr0->fill(m_vr0);


    if(fabs(m_vz0) >= m_vz0cut) continue;
    if(fabs(m_vr0) >= m_vr0cut) continue;
    if(mdcTrk->pxy() < m_ptcut) continue;
    if(mdcTrk->p() > m_momcut) continue;
    if(fabs(cos(mdcTrk->theta())) > m_cha_costheta_cut) continue;


    iGood.push_back(i);
    nCharge += mdcTrk->charge();
    // book e/p ratio and dE/dx
    if(mdcTrk->p()<m_distin_pionlep){
      m_pion_matched ++;
      if(mdcTrk->charge()>0){
	if((*itTrk)->isEmcShowerValid()) m_ep_ratio[0] = (*itTrk)->emcShower()->energy()/mdcTrk->p();
	if((*itTrk)->isMucTrackValid()) m_muc_dep[0] = (*itTrk)->mucTrack()->depth();
	f_get_pid( *itTrk, m_prob_e, m_prob_mu, m_prob_pi, 0); 
	f_get_dedx( *itTrk, m_dedx_e, m_dedx_mu, m_dedx_pi, m_dedx_ghits, 0); 
	f_get_tof( *itTrk, m_tof_e, m_tof_mu, m_tof_pi, 0); 
      }
      else{
	  if((*itTrk)->isEmcShowerValid()) m_ep_ratio[1] = (*itTrk)->emcShower()->energy()/mdcTrk->p(); 
	  if((*itTrk)->isMucTrackValid()) m_muc_dep[1] = (*itTrk)->mucTrack()->depth();
	  f_get_pid( *itTrk, m_prob_e, m_prob_mu, m_prob_pi, 1); 
	  f_get_dedx( *itTrk, m_dedx_e, m_dedx_mu, m_dedx_pi, m_dedx_ghits, 1); 
	  f_get_tof( *itTrk, m_tof_e, m_tof_mu, m_tof_pi, 1); 
      }// end of seperate +/-
    }// end of pion
    else{ 
	m_lep_matched ++;
	if(mdcTrk->charge()>0){
	  if((*itTrk)->isEmcShowerValid()) m_ep_ratio[2] = (*itTrk)->emcShower()->energy()/mdcTrk->p(); 
	  if((*itTrk)->isMucTrackValid()) m_muc_dep[2] = (*itTrk)->mucTrack()->depth();
	  f_get_pid( *itTrk, m_prob_e, m_prob_mu, m_prob_pi, 2);
	  f_get_dedx( *itTrk, m_dedx_e, m_dedx_mu, m_dedx_pi, m_dedx_ghits, 2); 
	  f_get_tof( *itTrk, m_tof_e, m_tof_mu, m_tof_pi, 2); 
 	}
	else{
	  if((*itTrk)->isEmcShowerValid()) m_ep_ratio[3] = (*itTrk)->emcShower()->energy()/mdcTrk->p(); 
	  if((*itTrk)->isMucTrackValid()) m_muc_dep[3] = (*itTrk)->mucTrack()->depth();
	  f_get_pid( *itTrk, m_prob_e, m_prob_mu, m_prob_pi, 3);
	  f_get_dedx( *itTrk, m_dedx_e, m_dedx_mu, m_dedx_pi, m_dedx_ghits, 3); 
	  f_get_tof( *itTrk, m_tof_e, m_tof_mu, m_tof_pi, 3); 
 	}// end of seperate +/-
    }// end of lep


    // book four momentum and vertex
    if(mdcTrk->charge()>0){
      if(mdcTrk->p()<m_distin_pionlep){
	mdcTrk->setPidType(RecMdcKalTrack::pion);
	if(m_kaon_flag) mdcTrk->setPidType(RecMdcKalTrack::kaon);
	if(m_num[0] == 0){
	  m_lv_pionp = mdcTrk->p4(xmass[2]);  
	  if(m_kaon_flag) m_lv_pionp = mdcTrk->p4(xmass[3]);  
	  igood_index[0] = i;
	}
	else{
	  m_lv_add_pi = mdcTrk->p4(xmass[2]);  
	  if(m_kaon_flag) m_lv_add_pi = mdcTrk->p4(xmass[3]);  
	  igood_index[4] = i;
	}
	m_num[0] ++;
	m_vertex[0][0] = m_vr0;
	m_vertex[0][1] = m_vz0;
      }
      else{
	if(m_num[2]==0){
	  mdcTrk->setPidType(RecMdcKalTrack::electron);
	  m_lv_pos = mdcTrk->p4(xmass[0]);
	  mdcTrk->setPidType(RecMdcKalTrack::muon);
	  m_lv_mup = mdcTrk->p4(xmass[1]);
	  igood_index[2] = i;
	}
	else{
	  mdcTrk->setPidType(RecMdcKalTrack::electron);
	  m_lv_add_e = mdcTrk->p4(xmass[0]);
	  mdcTrk->setPidType(RecMdcKalTrack::muon);
	  m_lv_add_u = mdcTrk->p4(xmass[1]);
	  igood_index[4] = i;
	}
	m_num[2] ++;
	m_vertex[2][0] = m_vr0;
	m_vertex[2][1] = m_vz0;
      }
    } // end of positive charge
    else{
      if(mdcTrk->p()<m_distin_pionlep){
	mdcTrk->setPidType(RecMdcKalTrack::pion);
	if(m_kaon_flag) mdcTrk->setPidType(RecMdcKalTrack::kaon);
	if(m_num[1]==0){
	  m_lv_pionm = mdcTrk->p4(xmass[2]);
 	  if(m_kaon_flag) m_lv_pionm = mdcTrk->p4(xmass[3]);
	  igood_index[1] = i;
	}
	else{
	  m_lv_add_pi = mdcTrk->p4(xmass[2]);
	  if(m_kaon_flag) m_lv_add_pi = mdcTrk->p4(xmass[3]);
	  igood_index[4] = i;
	}
	m_num[1] ++;
	m_vertex[1][0] = m_vr0;
	m_vertex[1][1] = m_vz0;
      }
      else{
	if(m_num[3]==0){
	  mdcTrk->setPidType(RecMdcKalTrack::electron);
	  m_lv_ele = mdcTrk->p4(xmass[0]); 
	  mdcTrk->setPidType(RecMdcKalTrack::muon);
	  m_lv_mum = mdcTrk->p4(xmass[1]); 
	  igood_index[3] = i;}
	else{
	  mdcTrk->setPidType(RecMdcKalTrack::electron);
	  m_lv_add_e = mdcTrk->p4(xmass[0]); 
	  mdcTrk->setPidType(RecMdcKalTrack::muon);
	  m_lv_add_u = mdcTrk->p4(xmass[1]); 
	  igood_index[4] = i;
	}
	m_num[3] ++;
	m_vertex[3][0] = m_vr0;
	m_vertex[3][1] = m_vz0;
      }
    } // end of negative charge
  }


  int nGood = iGood.size();
  log << MSG::DEBUG << "With KalmanTrack, ngood, totcharge = " 
      << nGood << " , " << nCharge << endreq;
  if(nGood<2 || nGood>5) return sc;
  m_cout_nGood ++;




  if(m_ep_ratio[2] > m_distin_emuon | m_ep_ratio[3] > m_distin_emuon){
    m_lv_lepp = m_lv_pos;
    m_lv_lepm = m_lv_ele;
    m_lv_add_lep = m_lv_add_e;
  }
  else{
    m_lv_lepp = m_lv_mup;
    m_lv_lepm = m_lv_mum;
    m_lv_add_lep = m_lv_add_u;
  }


  HepLorentzVector m_lv_lab(0.04,0,0.004,m_cms_ene);
  if(m_only_four & (nGood!=4)) return sc;
  if(nGood==4){
    if(nCharge) return sc;
    m_event_flag = 4;
  }
  else if(nGood==3){ 
    if(m_num[0]==0 && m_num[1]==1 && m_num[2]==1 && m_num[3]==1){
      if(nCharge != -1) return sc;
      m_lv_pionp = m_lv_lab - m_lv_pionm - m_lv_lepp - m_lv_lepm;
      if(fabs(m_lv_pionp.vect().cosTheta())>m_cha_costheta_cut)return sc;
      m_event_flag = 0;
    }
    if(m_num[0]==1 && m_num[1]==0 && m_num[2]==1 && m_num[3]==1){
      if(nCharge != 1) return sc;
      m_lv_pionm = m_lv_lab - m_lv_pionp - m_lv_lepp - m_lv_lepm;
      if(fabs(m_lv_pionm.vect().cosTheta())>m_cha_costheta_cut)return sc;
      m_event_flag = 1;
    }
    if(m_num[0]==1 && m_num[1]==1 && m_num[2]==0 && m_num[3]==1){
      if(nCharge != -1) return sc;
      m_lv_lepp = m_lv_lab - m_lv_pionp - m_lv_pionm - m_lv_lepm;
      if(fabs(m_lv_lepp.vect().cosTheta())>m_cha_costheta_cut)return sc;
      m_event_flag = 2;
    }
    if(m_num[0]==1 && m_num[1]==1 && m_num[2]==1 && m_num[3]==0){
      if(nCharge != 1) return sc;
      m_lv_lepm = m_lv_lab - m_lv_pionp - m_lv_pionm - m_lv_lepp;
      if(fabs(m_lv_lepm.vect().cosTheta())>m_cha_costheta_cut) return sc;
      m_event_flag = 3;
    }
  } // end of nGood = 3
  else if(nGood==5){
    if(m_num[0]==2 && m_num[1]==1 && m_num[2]==1 && m_num[3]==1){
      if(nCharge != 1) return sc;
      if(fabs((m_lv_lab-m_lv_pionm-m_lv_pionp).m()-3.097) > fabs((m_lv_lab-m_lv_pionm-m_lv_add_pi).m()-3.097) ){
	m_lv_pionp = m_lv_add_pi;
	igood_index[0] = igood_index[4];
      }
      if(fabs(m_lv_pionp.vect().cosTheta())>m_cha_costheta_cut)return sc;
      m_event_flag = 5;
    }
    if(m_num[0]==1 && m_num[1]==2 && m_num[2]==1 && m_num[3]==1){
      if(nCharge != -1) return sc;
      if(fabs((m_lv_lab-m_lv_pionm-m_lv_pionp).m()-3.097) > fabs((m_lv_lab-m_lv_pionp-m_lv_add_pi).m()-3.097) ) {
	m_lv_pionm = m_lv_add_pi;
	igood_index[1] = igood_index[4];
      }
      if(fabs(m_lv_pionm.vect().cosTheta())>m_cha_costheta_cut)return sc;
      m_event_flag = 6;
    }
    if(m_num[0]==1 && m_num[1]==1 && m_num[2]==2 && m_num[3]==1){
      if(nCharge != 1) return sc;
      if(fabs((m_lv_lepp+m_lv_lepm).m()-3.097) > fabs((m_lv_add_lep+m_lv_lepm).m()-3.097) ){
	m_lv_lepp = m_lv_add_lep;
	igood_index[2] = igood_index[4];
      }
      if(fabs(m_lv_lepp.vect().cosTheta())>m_cha_costheta_cut)return sc;
      m_event_flag = 7;
    }
    if(m_num[0]==1 && m_num[1]==1 && m_num[2]==1 && m_num[3]==2){
      if(nCharge != -1) return sc;
      if(fabs((m_lv_lepp+m_lv_lepm).m()-3.097) > fabs((m_lv_lepp+m_lv_add_lep).m()-3.097) ) {
	m_lv_lepm = m_lv_add_lep;
	igood_index[3] = igood_index[4];
      }
      if(fabs(m_lv_lepm.vect().cosTheta())>m_cha_costheta_cut)return sc;
      m_event_flag = 8;
    }
  } // end of nGood = 5
  else{
    if(nCharge) return sc;
    if((m_num[2]!=1) || (m_num[3]!=1)) return sc;
    m_event_flag = 10; // only two tracks
  } // end of nGood = 2
  m_cout_mom ++;
  //  cout << "evt falg: " << m_event_flag << endl;


  // With momentum method calculate the invariant mass of Jpsi
  // actually we use the recoil mass
  HepLorentzVector m_lv_recoil, m_lv_jpsi;
  m_lv_recoil = m_lv_lab - m_lv_pionp - m_lv_pionm;
  m_lv_jpsi = m_lv_lepp + m_lv_lepm;


  m_mass_twopi = (m_lv_pionp + m_lv_pionm).m();
  m_mass_recoil = m_lv_recoil.m();
  m_mass_jpsi = m_lv_jpsi.m();


  // Jpsi mass cut
  if(nGood>2 && ( m_mass_recoil < 2.5 || m_mass_recoil > 4.5) ) return sc;
  if( m_mass_jpsi < 2.5 || m_mass_jpsi > 4.0 ) return sc;
  m_cout_recoil ++;


  HepLorentzVector m_ttm(m_lv_jpsi + m_lv_pionp + m_lv_pionm);
  if(nGood>2 && (m_ttm.m()>5 || m_ttm.m()<3)) return sc;
 
  // dangle between pions, suppress gamma convertion
  m_pipi_dang = m_lv_pionp.vect().cosTheta(m_lv_pionm.vect());


  m_mom_pionp = m_lv_pionp.vect().mag();
  m_mom_pionm = m_lv_pionm.vect().mag();
  m_mom_lepp = m_lv_lepp.vect().mag();
  m_mom_lepm = m_lv_lepm.vect().mag();
  m_pt_lepp = m_lv_lepp.vect().perp();
  m_pt_lepm = m_lv_lepm.vect().perp();
  m_pt_pionp = m_lv_pionp.vect().perp();
  m_pt_pionm = m_lv_pionm.vect().perp();


  Hep3Vector m_boost_jpsi(m_lv_recoil.boostVector());
  if(m_debug_flag) cout << "before boost: " << m_lv_lepp << endl;
  HepLorentzVector m_lv_cms_lepp(boostOf(m_lv_lepp,-m_boost_jpsi));
  HepLorentzVector m_lv_cms_lepm(boostOf(m_lv_lepm,-m_boost_jpsi));
  if(m_debug_flag) cout << "after boost: " << m_lv_lepp << endl;
  m_cms_index = 3;
  m_cms_lepm[0] = m_lv_cms_lepm.vect().mag();
  m_cms_lepp[0] = m_lv_cms_lepp.vect().mag();
  m_cms_lepm[1] = m_lv_cms_lepm.vect().theta();
  m_cms_lepp[1] = m_lv_cms_lepp.vect().theta();
  m_cms_lepm[2] = m_lv_cms_lepm.vect().phi();
  m_cms_lepp[2] = m_lv_cms_lepp.vect().phi();
  if(m_debug_flag) cout << "jpsi four momentum in cms " 
			<< m_lv_cms_lepp + m_lv_cms_lepm << endl;


  m_inv_mass = m_ttm.m();
  m_tot_e = m_ttm.e();
  m_tot_px = m_ttm.px();
  m_tot_py = m_ttm.py();
  m_tot_pz = m_ttm.pz();
  HepLorentzVector m_lv_book(0,0,0,0), m_lv_recoil_book(0,0,0,0);
  for(m_index=0; m_index<4; m_index++){
    switch(m_index){
    case 0: 
      m_lv_book = m_lv_pionp; 
      m_lv_recoil_book = m_lv_lab - m_lv_pionm - m_lv_lepp - m_lv_lepm;
      break;
    case 1: 
      m_lv_book = m_lv_pionm; 
      m_lv_recoil_book = m_lv_lab - m_lv_pionp - m_lv_lepp - m_lv_lepm;
      break;
    case 2: 
      m_lv_book = m_lv_lepp;
      m_lv_recoil_book = m_lv_lab - m_lv_pionp - m_lv_pionm - m_lv_lepm;
      break;
    case 3: 
      m_lv_book = m_lv_lepm;
      m_lv_recoil_book = m_lv_lab - m_lv_pionp - m_lv_pionm - m_lv_lepp;
      break;
    default: m_lv_book.setE(2008); // 2008 is a joke
    }
    m_cos_theta[m_index] = m_lv_book.vect().cosTheta();
    m_phi[m_index] = m_lv_book.vect().phi();
    m_four_mom[m_index][0] = m_lv_book.e();
    m_four_mom[m_index][1] = m_lv_book.px();
    m_four_mom[m_index][2] = m_lv_book.py();
    m_four_mom[m_index][3] = m_lv_book.pz();
    m_recoil_four_mom[m_index][0] = m_lv_recoil_book.e();
    m_recoil_four_mom[m_index][1] = m_lv_recoil_book.px();
    m_recoil_four_mom[m_index][2] = m_lv_recoil_book.py();
    m_recoil_four_mom[m_index][3] = m_lv_recoil_book.pz();
    m_miss_mass[m_index] = m_lv_recoil_book.m2();
    if(m_debug_flag) 
      cout << "recoil mass: " << m_lv_recoil_book.m() << 
	" self: " << m_lv_book.m() << " eve: " << m_event_flag << 
	" in cms: " << m_lv_cms_lepp.m() << endl;
    m_how_near[m_index] = m_lv_book.vect().howNear(m_lv_recoil_book.vect());
  }


  //MC information
  if(m_run<0 && !m_track_flag){
    int m_numParticle(0), m_true_pid(0);
    if(!mcParticleCol){
      log << MSG::ERROR << "Could not retrieve McParticelCol" << endreq;
      return StatusCode::FAILURE;
    }
    else{
      bool psipDecay(false);
      int rootIndex(-1);
      HepLorentzVector m_lv_true_pionp, m_lv_true_pionm, m_lv_true_lepp, m_lv_true_lepm;
      Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
      for (; iter_mc != mcParticleCol->end(); iter_mc++){
	if(m_debug_flag) cout << "pid: " << (*iter_mc)->particleProperty() << endl;
        if ((*iter_mc)->primaryParticle()) continue;
        if (!(*iter_mc)->decayFromGenerator()) continue;
        //if ( ((*iter_mc)->mother()).trackIndex()<3 ) continue;
	long temp_id((*iter_mc)->particleProperty());
	if(temp_id==110441 || temp_id==120443 || temp_id == 100445 || (temp_id>=51 && temp_id<=55)) m_xyzid = temp_id;
        if( (fabs(m_cms_ene-3.686)<0.01 && (*iter_mc)->particleProperty()==100443) || (fabs(m_cms_ene-3.097)<0.01 && (*iter_mc)->particleProperty()==443) || (fabs(m_cms_ene-3.770)<0.01 && (*iter_mc)->particleProperty()==30443) || (fabs(m_cms_ene-4.040)<0.04 && (*iter_mc)->particleProperty()==9000443) ){
          psipDecay = true;
          rootIndex = (*iter_mc)->trackIndex();
        }
        if (!psipDecay) continue;
        int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
        int pdgid = (*iter_mc)->particleProperty();
        m_pdgid[m_numParticle] = pdgid;
        m_motheridx[m_numParticle] = mcidx;
        m_numParticle ++;    


	//if(!(*iter_mc)->leafParticle()) continue;
	if((*iter_mc)->particleProperty() == 211) {
	  m_lv_true_pionp = (*iter_mc)->initialFourMomentum();
	}
	if((*iter_mc)->particleProperty() == -211) {
	  m_lv_true_pionm = (*iter_mc)->initialFourMomentum();
	}
	if((*iter_mc)->particleProperty() == -11 || (*iter_mc)->particleProperty() == -13) {
	  m_lv_true_lepp = (*iter_mc)->initialFourMomentum();
	}
	if((*iter_mc)->particleProperty() ==  11 || (*iter_mc)->particleProperty() ==  13) {
	  m_lv_true_lepm = (*iter_mc)->initialFourMomentum();
	}
      }
      m_true_cms_index = 3;
      Hep3Vector m_true_boost_jpsi((m_lv_true_lepp+m_lv_true_lepm).boostVector());
      HepLorentzVector m_lv_true_cms_lepp(boostOf(m_lv_true_lepp, -m_true_boost_jpsi));
      HepLorentzVector m_lv_true_cms_lepm(boostOf(m_lv_true_lepm, -m_true_boost_jpsi));
      m_true_cms_lepp[0] = m_lv_true_cms_lepp.vect().mag();
      m_true_cms_lepm[0] = m_lv_true_cms_lepm.vect().mag();
      m_true_cms_lepp[1] = m_lv_true_cms_lepp.vect().theta();
      m_true_cms_lepm[1] = m_lv_true_cms_lepm.vect().theta();
      m_true_cms_lepp[2] = m_lv_true_cms_lepp.vect().phi();
      m_true_cms_lepm[2] = m_lv_true_cms_lepm.vect().phi();
      m_idxmc = m_numParticle;
    }
  }


  if(m_subsample_flag &&  m_event_flag==4 && m_mass_recoil>3.08 && m_mass_recoil<3.12 && m_mass_jpsi>3.05 && m_mass_jpsi<3.15 && m_pipi_dang<0.98 && m_cms_lepm[0]<1.7 && m_cms_lepp[0]<1.7 && m_cms_lepm[0]>1.4 && m_cms_lepp[0]>1.4 && m_mom_pionm<0.5 && m_mom_pionp<0.5 && (m_pt_pionm<0.2 || m_pt_pionp<0.2)) setFilterPassed(true);
  // else if(m_event_flag==0 && m_true_pionp>0.35 && m_true_pionp<0.4)
  //  setFilterPassed(true);
  //cout << "passed" << endl;
  //  m_tuple1->write();
  


  // next is good photon selection
  Vint iGam;
  Vp4 m_Gam;
  for(int i=evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++) {
    if(i>50) break;
    EvtRecTrackIterator itTrk = m_itTrk_begin + i;
    if(!(*itTrk)->isEmcShowerValid()) continue;
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    double m_abs_costhe(fabs(cos(emcTrk->theta())));
    if(m_abs_costhe > m_bar_costheta_cut && m_abs_costhe < m_min_end_costheta_cut) continue;
    if(m_abs_costhe > m_max_end_costheta_cut) continue;
    double eraw = emcTrk->energy();
    log << MSG::DEBUG << "i= " << i << "  ene= " << eraw << endreq;
    if(m_abs_costhe < m_bar_costheta_cut){
      if(eraw < m_bar_energy_cut) continue;
    }
    else if(eraw < m_end_energy_cut) continue;
    if(m_run > 0){
      if(emcTrk->time() < m_min_emctime || emcTrk->time() > m_max_emctime) continue;
    }
    Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
    // find the nearest charged track
    double dthe = 200.;
    double dphi = 200.;
    double dang = 200.; 
    for(int j = 0; j < evtRecEvent->totalCharged(); j++) {
      EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
      if(!(*jtTrk)->isExtTrackValid()) continue;
      RecExtTrack *extTrk = (*jtTrk)->extTrack();
      if(extTrk->emcVolumeNumber() == -1) continue;
      Hep3Vector extpos = extTrk->emcPosition();
      //      double ctht = extpos.cosTheta(emcpos);
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
    if(dang < m_dang_cut ) continue;


    iGam.push_back(i);
    if(iGam.size()==1) m_max_gam_dang = dang;
    HepLorentzVector m_lv_temp(emcpos.unit()*=eraw,eraw);
    m_Gam.push_back(m_lv_temp);
  }


  // Finish Good Photon Selection
  m_nGam = iGam.size(); 
  if(m_nGam>0){
    m_max_gam_ene = m_Gam[0].e();
    m_max_gam_costhe = m_Gam[0].cosTheta();
    m_max_gam_phi = m_Gam[0].vect().phi();
  }
  else m_max_gam_ene = 0;
  log << MSG::DEBUG << "num Good Photon " 
      << m_nGam  << " , " <<evtRecEvent->totalNeutral()<<endreq;
  HepLorentzVector m_lv_pi0;
  if(m_nGam>2){
    PionZeroList m_vpi0_list(m_Gam);
    if(m_debug_flag) cout << "run " << m_run << " event " << m_event << endl;
    if(m_debug_flag) cout << "pion0 list before sort " << endl;
    if(m_debug_flag) m_vpi0_list.print();
    m_vpi0_list.sort();
    if(m_debug_flag) cout << "pion0 list after sort, before reduce " << endl;
    if(m_debug_flag) m_vpi0_list.print();
    m_vpi0_list.reduce();
    if(m_debug_flag) cout << "pion0 list after reduce " << endl;
    if(m_debug_flag) m_vpi0_list.print();
    m_num_pi0 = m_vpi0_list.get_num_pi0();
    if(m_num_pi0>0){
      m_like_pi0_mass = (m_vpi0_list.get_pi0_list())[0].m();
      m_lv_pi0 = (m_vpi0_list.get_pi0_list())[0];
    }
    else m_like_pi0_mass = 0;
    if(m_debug_flag) cout << "likest pi0 mass " << m_like_pi0_mass << endl;
  }
  else{
    //  cout << "num of photons less than 2 " << endl;
    m_num_pi0 = 0;
    m_like_pi0_mass = 0;
  }


  //====================================
  // Kinematic 1C fit is next
  // Firstly, we do vertex fit for charged tracks
  //====================================
  if(m_fit_flag && (m_event_flag==0 || m_event_flag==1 || m_event_flag==4 || m_event_flag==5 || m_event_flag==6 || m_event_flag==10)){
    RecMdcKalTrack *m_lepp_track = (*(m_itTrk_begin + igood_index[2]))->mdcKalTrack();
    RecMdcKalTrack *m_lepm_track = (*(m_itTrk_begin + igood_index[3]))->mdcKalTrack();
    WTrackParameter m_lepp_par, m_lepm_par;
    if(m_emc_ene > m_distin_emuon){
      m_lepp_par = WTrackParameter(xmass[0], m_lepp_track->getZHelixE(), m_lepp_track->getZErrorE());
      m_lepm_par = WTrackParameter(xmass[0], m_lepm_track->getZHelixE(), m_lepm_track->getZErrorE());
    }
    else {
      m_lepp_par = WTrackParameter(xmass[1], m_lepp_track->getZHelixMu(), m_lepp_track->getZErrorMu());
      m_lepm_par = WTrackParameter(xmass[1], m_lepm_track->getZHelixMu(), m_lepm_track->getZErrorMu());
    }
    /* Default is pion, for other particles:
       wvppTrk = WTrackParameter(mp, pipTrk->getZHelixP(), pipTrk->getZErrorP());
       wvmupTrk = WTrackParameter(mmu,pipTrk->getZHelixMu(), pipTrk->getZErrorMu());
       wvepTrk = WTrackParameter(me, pipTrk->getZHelixE(), pipTrk->getZErrorE());
       wvkpTrk = WTrackParameter(mk, pipTrk->getZHelixK(), pipTrk->getZErrorK());*/
    HepPoint3D vx(0,0,0);
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
    log << MSG::DEBUG << "before vertex fit init" << endreq;
    VertexFit* vtxfit = VertexFit::instance();
    vtxfit->init();
    log << MSG::DEBUG << "after vertex fit init" << endreq;
    vtxfit->AddTrack(0, m_lepp_par);
    vtxfit->AddTrack(1, m_lepm_par);
    vtxfit->AddVertex(0, vxpar, 0,1);
    f_chisq_vt = 9999.;
    f_chisq_jpsi = 9999.;
    f_chisq_gammajpsi = 9999.;
    if(vtxfit->Fit(0)){
      f_chisq_vt = vtxfit->chisq(0);
      vtxfit->Swim(0);
    }
    // secondly, we do kinematic 1C fit
    log << MSG::DEBUG << "before 1C fit init" << endreq;
    KinematicFit *kmfit = KinematicFit::instance();
    kmfit->init();
    log << MSG::DEBUG << "after 1C fit init" << endreq;
    kmfit->AddTrack(0, vtxfit->wtrk(0));
    kmfit->AddTrack(1, vtxfit->wtrk(1));
    kmfit->AddResonance(0, 3.097, 0, 1);
    if(kmfit->Fit()) f_chisq_jpsi = kmfit->chisq();
    log << MSG::DEBUG << "after 1C fit" << endreq;
    m_lv_lepp = kmfit->pfit(0);
    m_lv_lepm = kmfit->pfit(1);
    if(m_event_flag==0 || m_event_flag==5) m_lv_pionp = m_lv_lab - m_lv_pionm - m_lv_lepp - m_lv_lepm;
    if(m_event_flag==1 || m_event_flag==6) m_lv_pionm = m_lv_lab - m_lv_pionp - m_lv_lepp - m_lv_lepm;
    f_mom_lepp = m_lv_lepp.rho();
    f_mom_lepm = m_lv_lepm.rho();
    f_mom_pionp = m_lv_pionp.rho();
    f_mom_pionm = m_lv_pionm.rho();
    f_pipi_dang = m_lv_pionp.vect().cosTheta(m_lv_pionm.vect());
    f_mass_twopi = (m_lv_pionp + m_lv_pionm).m();
    f_mass_jpsi = (m_lv_lepp + m_lv_lepm).m();
    m_lv_recoil = m_lv_lab - m_lv_pionp - m_lv_pionm;
    f_mass_recoil = m_lv_recoil.m();
    f_inv_mass = (m_lv_lepp + m_lv_lepm + m_lv_pionp + m_lv_pionm).m();
    f_tot_e =  (m_lv_lepp + m_lv_lepm + m_lv_pionp + m_lv_pionm).e();
    f_tot_px =  (m_lv_lepp + m_lv_lepm + m_lv_pionp + m_lv_pionm).px();
    f_tot_py =  (m_lv_lepp + m_lv_lepm + m_lv_pionp + m_lv_pionm).py();
    f_tot_pz =  (m_lv_lepp + m_lv_lepm + m_lv_pionp + m_lv_pionm).pz();
    f_pt_lepp = m_lv_lepp.vect().perp();
    f_pt_lepm = m_lv_lepm.vect().perp();
    f_pt_pionp = m_lv_pionp.vect().perp();
    f_pt_pionm = m_lv_pionm.vect().perp();
    m_boost_jpsi = m_lv_recoil.boostVector();
    m_lv_cms_lepp = boostOf(m_lv_lepp,-m_boost_jpsi);
    m_lv_cms_lepm = boostOf(m_lv_lepm,-m_boost_jpsi);
    f_cms_lepm = m_lv_cms_lepm.rho();
    f_cms_lepp = m_lv_cms_lepp.rho();


    HepLorentzVector f_lv_book(0,0,0,0), f_lv_recoil_book(0,0,0,0);
    for(f_index=0; f_index<4; f_index++){
      switch(f_index){
      case 0: 
	f_lv_book = m_lv_pionp; 
	f_lv_recoil_book = m_lv_lab - m_lv_pionm - m_lv_lepp - m_lv_lepm;
	break;
      case 1: 
	f_lv_book = m_lv_pionm; 
	f_lv_recoil_book = m_lv_lab - m_lv_pionp - m_lv_lepp - m_lv_lepm;
	break;
      case 2: 
	f_lv_book = m_lv_lepp;
	f_lv_recoil_book = m_lv_lab - m_lv_pionp - m_lv_pionm - m_lv_lepm;
	break;
      case 3: 
	f_lv_book = m_lv_lepm;
	f_lv_recoil_book = m_lv_lab - m_lv_pionp - m_lv_pionm - m_lv_lepp;
	break;
      default: f_lv_book.setE(2008); // 2008 is a joke
      }
      f_cos_theta[f_index] = f_lv_book.vect().cosTheta();
      f_phi[f_index] = f_lv_book.vect().phi();
      f_four_mom[f_index][0] = f_lv_book.e();
      f_four_mom[f_index][1] = f_lv_book.px();
      f_four_mom[f_index][2] = f_lv_book.py();
      f_four_mom[f_index][3] = f_lv_book.pz();
      f_recoil_four_mom[f_index][0] = f_lv_recoil_book.e();
      f_recoil_four_mom[f_index][1] = f_lv_recoil_book.px();
      f_recoil_four_mom[f_index][2] = f_lv_recoil_book.py();
      f_recoil_four_mom[f_index][3] = f_lv_recoil_book.pz();
      f_miss_mass[f_index] = f_lv_recoil_book.m2();
      f_how_near[f_index] = f_lv_book.vect().howNear(f_lv_recoil_book.vect());
    } // end of book matrix variables
    if(m_num_pi0>0){
      m_pipipi0_mass = (m_lv_pi0 + m_lv_pionp + m_lv_pionm).m();
      m_omegajpsi_mass = (m_lv_pi0 + m_lv_pionp + m_lv_pionm + m_lv_lepp + m_lv_lepm).m();
    }
    else {
      m_pipipi0_mass = 0;
      m_omegajpsi_mass = 0;
    } 
    if(nGood==2 && m_nGam>0) m_gammajpsi_mass = (m_Gam[0] + m_lv_lepp + m_lv_lepm).m();
    else m_gammajpsi_mass = 0;
    if(m_event_flag==10 & m_nGam>1){ // do 4C fit here
      kmfit->init();
      kmfit->AddTrack(0, vtxfit->wtrk(0));
      kmfit->AddTrack(1, vtxfit->wtrk(1));
      RecEmcShower *gTrk_1 = (*(m_itTrk_begin+iGam[0]))->emcShower();
      RecEmcShower *gTrk_2 = (*(m_itTrk_begin+iGam[1]))->emcShower();
      kmfit->AddTrack(2, 0., gTrk_1);
      kmfit->AddTrack(3, 0., gTrk_2);
      kmfit->AddFourMomentum(0, m_lv_lab);
      if(kmfit->Fit(0)){
	f_chisq_gammajpsi = kmfit->chisq(0);
	f_gammajpsi_mass = (kmfit->pfit(0) + kmfit->pfit(1) + kmfit->pfit(2)).m();
      } 
    }// end of 4C fit of gamma gamma Jpsi
  }// end of 1C fit of J/psi
  log << MSG::DEBUG << "after 1C fit" << endreq;




  //
  // check dedx infomation
  //
  if(m_checkDedx) {
    int m_dedx_cout(0);
    for(int i = 0; i < nGood; i++) {
      EvtRecTrackIterator  itTrk = evtRecTrkCol->begin() + iGood[i];
      if(!(*itTrk)->isMdcDedxValid())continue;
      RecMdcKalTrack *mdcTrk = (*itTrk)->mdcKalTrack();
      RecMdcDedx *dedxTrk = (*itTrk)->mdcDedx();
  
      m_ptrk = mdcTrk->p();
      m_chie = dedxTrk->chiE();
      m_chimu = dedxTrk->chiMu();
      m_chipi = dedxTrk->chiPi();
      m_chik = dedxTrk->chiK();
      m_chip = dedxTrk->chiP();
      m_ghit = dedxTrk->numGoodHits();
      m_thit = dedxTrk->numTotalHits();
      m_probPH = dedxTrk->probPH();
      m_normPH = dedxTrk->normPH();


    }
  }


  //
  // check TOF infomation
  //
  if(m_checkTof) {
    int m_endcap_cout(0), m_layer1_cout(0), m_layer2_cout(0);
    TofHitStatus *status = new TofHitStatus; 
    for(int i = 0; i < nGood; i++) {
      EvtRecTrackIterator  itTrk = evtRecTrkCol->begin() + iGood[i];
      if(!(*itTrk)->isTofTrackValid()) continue;


      RecMdcKalTrack *mdcTrk = (*itTrk)->mdcKalTrack();
      SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();


      double ptrk = mdcTrk->p();


      for( SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
	   iter_tof != tofTrkCol.end(); iter_tof++ ) { 
        status->setStatus((*iter_tof)->status());
        if(!(status->is_barrel())){//endcap
          if( !(status->is_counter()) ) continue; // ? 
          if( status->layer()!=0 ) continue;//layer1
          double path = (*iter_tof)->path(); // ? the unit is cm
          double tof  = (*iter_tof)->tof();  // the unit is ns/100
          double ph   = (*iter_tof)->ph();
          double rhit = (*iter_tof)->zrhit();
          double qual = 0.0 + (*iter_tof)->quality();
          double cntr = 0.0 + (*iter_tof)->tofID();
          double texp[5];
          for(int j = 0; j < 5; j++) {
            double gb = xmass[j]/ptrk;
            double beta = sqrt(1+gb*gb);
            texp[j] = path*beta/velc; // the unit here is ns
          }
          m_cntr_etof  = cntr;
          m_ptot_etof  = ptrk;
	  m_path_etof = path;
          m_ph_etof    = ph;
          m_rhit_etof  = rhit;
          m_qual_etof  = qual;
	  m_tof_etof = tof;
          m_te_etof    = tof - texp[0];
          m_tmu_etof   = tof - texp[1];
          m_tpi_etof   = tof - texp[2];
          m_tk_etof    = tof - texp[3];
          m_tp_etof    = tof - texp[4];


        }
        else {//barrel
          if( !(status->is_counter()) ) continue; // ? 
          if(status->layer()==1){ //layer1
            double path=(*iter_tof)->path(); // ? 
            double tof  = (*iter_tof)->tof();
            double ph   = (*iter_tof)->ph();
            double rhit = (*iter_tof)->zrhit();
            double qual = 0.0 + (*iter_tof)->quality();
            double cntr = 0.0 + (*iter_tof)->tofID();
            double texp[5];
            for(int j = 0; j < 5; j++) {
	      double gb = xmass[j]/ptrk;
	      double beta = sqrt(1+gb*gb);
	      texp[j] = path*beta/velc;
            }
 
            m_cntr_btof1  = cntr;
            m_ptot_btof1 = ptrk;
	    m_path_btof1 = path;
            m_ph_btof1   = ph;
            m_zhit_btof1  = rhit;
            m_qual_btof1  = qual;
	    m_tof_btof1 = tof;
            m_te_btof1    = tof - texp[0];
            m_tmu_btof1   = tof - texp[1];
            m_tpi_btof1   = tof - texp[2];
            m_tk_btof1    = tof - texp[3];
            m_tp_btof1    = tof - texp[4];



			}


          if(status->layer()==2){//layer2
            double path=(*iter_tof)->path(); // ? 
            double tof  = (*iter_tof)->tof();
            double ph   = (*iter_tof)->ph();
            double rhit = (*iter_tof)->zrhit();
            double qual = 0.0 + (*iter_tof)->quality();
            double cntr = 0.0 + (*iter_tof)->tofID();
            double texp[5];
            for(int j = 0; j < 5; j++) {
	      double gb = xmass[j]/ptrk;
	      double beta = sqrt(1+gb*gb);
	      texp[j] = path*beta/velc;
            }
 
            m_cntr_btof2  = cntr;
            m_ptot_btof2 = ptrk;
	    m_path_btof2 = path;
            m_ph_btof2   = ph;
            m_zhit_btof2  = rhit;
            m_qual_btof2  = qual;
	    m_tof_btof2 = tof;
            m_te_btof2    = tof - texp[0];
            m_tmu_btof2   = tof - texp[1];
            m_tpi_btof2   = tof - texp[2];
            m_tk_btof2    = tof - texp[3];
            m_tp_btof2    = tof - texp[4];


          } 
        }


      } 
    } // loop all charged track
    delete status; 
  }  // check tof
 
TreeAna->Fill();


  return sc;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode PipiJpsi::finalize() {

cout<<"This finalize for chenq!"<<endl;


    NbInfo->Fill();
	saveFile->cd();
	TreeAna->Write();
	if(m_saveTopoTree == 1)	TopoTree->Write();
	if(m_saveNbInfo == 1) NbInfo->Write();
	saveFile->Close();
	
	
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in finalize()" << endmsg;
  if(m_eventrate) cout << "all event: " << m_cout_all << endl 
		       << "only channel " << m_chan_det << ": " 
		       << m_cout_everat << endl;
  cout << "all event: " << m_cout_all << endl 
       << "all collection point is OK: " << m_cout_col << endl 
       << "charged tracks >=3: " << m_cout_charge << endl 
       << "good charged tracks [3,4]: " << m_cout_nGood << endl 
       << "after momentum assign: " << m_cout_mom << endl 
       << "after recoild mass cut: " << m_cout_recoil << endl;
	   cout<<"test0512"<<endl;
  return StatusCode::SUCCESS;
}// end of finalize()


 
void PipiJpsi::f_get_pid(EvtRecTrack *recTrk, double m_prob_e[500], double m_prob_mu[500], double m_prob_pi[500], int n){
  ParticleID *pid = ParticleID::instance();
  pid->init();
  pid->setRecTrack(recTrk);
  pid->setMethod(pid->methodProbability());
  pid->setChiMinCut(4.0);
  // use PID sub-system
  pid->usePidSys( pid->useDedx() );
  pid->usePidSys( pid->useTof1() );
  pid->usePidSys( pid->useTof2() );
  pid->usePidSys( pid->useTofE() );
  pid->usePidSys( pid->useEmc() );
  pid->usePidSys( pid->useMuc() );
  pid->identify( pid->onlyElectron() | pid->onlyMuon() | pid->onlyPion() );
  pid->calculate();
  m_prob_e[n] = pid->probElectron();
  m_prob_mu[n]= pid->probMuon(); 
  m_prob_pi[n] = pid->probPion();
}// enf of f_get_pid()




void PipiJpsi::f_get_dedx(EvtRecTrack *itTrk, double m_dedx_e[500], double m_dedx_mu[500], double m_dedx_pi[500],  double m_dedx_ghits[500], int n){
  m_dedx_e[n] = m_dedx_mu[n] = m_dedx_pi[n] = -99.;  
  m_dedx_ghits[n] = 0;
  if((itTrk)->isMdcDedxValid()){
    RecMdcDedx *dedxTrk = (itTrk)->mdcDedx();
    m_dedx_e[n] = dedxTrk->chiE();
    m_dedx_mu[n] = dedxTrk->chiMu();
    m_dedx_pi[n] = dedxTrk->chiPi();
    m_dedx_ghits[n] = dedxTrk->numGoodHits();
  }
}// end of f_get_dedx()




void PipiJpsi::f_get_tof(EvtRecTrack *itTrk, double m_tof_e[500], double m_tof_mu[500], double m_tof_pi[500], int n){
  int m_multi = 0;
  m_tof_e[n] = m_tof_mu[n] = m_tof_pi[n] = 0.;  
  if(!(itTrk)->isTofTrackValid()) return;
  RecMdcKalTrack *mdcTrk = (itTrk)->mdcKalTrack();
  SmartRefVector<RecTofTrack> tofTrkCol = (itTrk)->tofTrack();
  double ptrk = mdcTrk->p();
     
  TofHitStatus *status = new TofHitStatus; 
  for( SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin(); iter_tof != tofTrkCol.end(); iter_tof++ ) { 
    status->setStatus((*iter_tof)->status());
    if(!(status->is_barrel())){//endcap
      if( !(status->is_counter()) ) continue; // ? 
      if( status->layer()!=0 ) continue;//layer1
      double path = (*iter_tof)->path(); // ? the unit is cm
      double tof  = (*iter_tof)->tof();  // the unit is ns/100
      double texp[5];
      for(int j = 0; j < 5; j++) {
	double gb = xmass[j]/ptrk;
	double beta = sqrt(1+gb*gb);
	texp[j] = path*beta/velc; // the unit here is ns
      }
      m_tof_e[n]  += (tof - texp[0]);
      m_tof_mu[n] += (tof - texp[1]);
      m_tof_pi[n] += (tof - texp[2]);
    }// end of endcap
    else {//barrel
      if( !(status->is_counter()) ) continue; // ? 
      if(status->layer()==1){ //layer1
	double path=(*iter_tof)->path(); // ? 
	double tof  = (*iter_tof)->tof();
	double texp[5];
	for(int j = 0; j < 5; j++) {
	  double gb = xmass[j]/ptrk;
	  double beta = sqrt(1+gb*gb);
	  texp[j] = path*beta/velc;
	}
	m_tof_e[n]  += (tof - texp[0]);
	m_tof_mu[n] += (tof - texp[1]);
	m_tof_pi[n] += (tof - texp[2]);
	m_multi ++;
      }// end of layer1
      if(status->layer()==2){//layer2
	double path=(*iter_tof)->path(); // ? 
	double tof  = (*iter_tof)->tof();
	double texp[5];
	for(int j = 0; j < 5; j++) {
	  double gb = xmass[j]/ptrk;
	  double beta = sqrt(1+gb*gb);
	  texp[j] = path*beta/velc;
	}
	m_tof_e[n]  += (tof - texp[0]);
	m_tof_mu[n] += (tof - texp[1]);
	m_tof_pi[n] += (tof - texp[2]);
	m_multi ++;
      }// end of layer2
    }// end of barrel
  }// end of loop in TOF
  delete status;
  
  if(m_multi<1) { m_tof_e[n] = m_tof_mu[n] = m_tof_pi[n] = -99.; }
  else {
    m_tof_e[n]  /= m_multi;
    m_tof_mu[n] /= m_multi;
    m_tof_pi[n] /= m_multi;
  }
}// end of f_get_tof()


