//psi'--> J/psi pion pion, J/psi --> di-leptons
//Kai Zhu (zhuk@ihep.ac.cn)
#ifndef Physics_Analysis_PipiJpsi_H
#define Physics_Analysis_PipiJpsi_H 

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "trkInfo.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogramFactory.h"
#include "EvtRecEvent/EvtRecTrack.h"
// Specify the namespace
using AIDA::IHistogram1D;
//#include "TH1F.h"

#include <string>
#include <TTree.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TClonesArray.h>

class PipiJpsi : public Algorithm {

public:
  PipiJpsi(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();  

private:
  void f_get_pid(EvtRecTrack *recTrk, double m_prob_e[500], double m_prob_mu[500], double m_prob_pi[500], int n);
  void f_get_dedx(EvtRecTrack *recTrk, double m_dedx_e[500], double m_dedx_mu[500], double m_dedx_pi[500], double m_dedx_ghits[500], int n);
  void f_get_tof(EvtRecTrack *recTrk, double m_tof_e[500], double m_tof_mu[500], double m_tof_pi[500], int n);

  // Declare r0, z0 cut for charged tracks
  double m_vr0cut, m_vz0cut, m_ptcut, m_momcut, m_cha_costheta_cut;
  double m_bar_energy_cut, m_end_energy_cut, m_dang_cut, m_bar_costheta_cut,
  m_min_end_costheta_cut, m_max_end_costheta_cut, m_min_emctime, m_max_emctime;

  // 
  bool m_checkDedx, m_checkTof, m_eventrate, m_fit_flag, m_debug_flag, m_non_fsr, m_vtdb_flag, m_only_four;
  int m_chan_det;
  // declare the track angle cut and track's transverse momentum
  double m_cosThetaCut;  

  // declare the cut angle between two pions to drop gamma conversion
  double m_pipi_dang_cut;

  // declare whether pick up sub-sample
  bool m_subsample_flag, m_trigger_flag, m_est_flag, m_premc;

  // declare energy/moentum to distinguish e and muon
  double m_distin_emuon;

  double m_cms_ene, m_distin_pionlep;

  // use PAT/TST instead of KalmanFilter to check the number of tracks
  // another function of this flag is to deal with all the mud on the 
  //track on which Xiao Jiao walks 
  bool m_track_flag, m_valid_flag;

  // whether analyze K+ K- J/psi instead of pi+ pi- J/psi
  bool m_kaon_flag; 

  // define Histogram here
 /* 
 IHistogram1D *h_vr0;
  IHistogram1D *h_vz0;
 */ 
  
//-------------------------file and MC -----------------------------------
  
  std::string m_OutputFileName;
		TFile *saveFile;
		TTree *TreeAna;
		TTree *genTree;
		TTree *TopoTree;
		TTree *NbInfo;
		
  int m_saveTopo;
  int m_saveMCTruth;
  int m_saveTopoTree;
  int m_saveNbInfo;
		
  int runid;
  int evtid;
  //int nevt;


  double m_dthe;
  double m_dphi;
  double m_dang;
  double m_eraw;
  long m_nGam;

  double m_ptrk;
  double m_chie;
  double m_chimu;
  double m_chipi;
  double m_chik;
  double m_chip;
  double m_probPH;
  double m_normPH;
  double m_ghit;
  double m_thit;


  double m_ptot_etof;
  double m_cntr_etof;
  double m_path_etof;
  double m_tof_etof;
  double m_te_etof;
  double m_tmu_etof;
  double m_tpi_etof;
  double m_tk_etof;
  double m_tp_etof;
  double m_ph_etof;
  double m_rhit_etof;
  double m_qual_etof;


  double m_ptot_btof1;
  double m_cntr_btof1;
  double m_path_btof1;
  double m_tof_btof1;
  double m_te_btof1;
  double m_tmu_btof1;
  double m_tpi_btof1;
  double m_tk_btof1;
  double m_tp_btof1;
  double m_ph_btof1;
  double m_zhit_btof1;
  double m_qual_btof1;


  double m_ptot_btof2;
  double m_cntr_btof2;
  double m_path_btof2;
  double m_tof_btof2;
  double m_te_btof2;
  double m_tmu_btof2;
  double m_tpi_btof2;
  double m_tk_btof2;
  double m_tp_btof2;
  double m_ph_btof2;
  double m_zhit_btof2;
  double m_qual_btof2;


  double m_pionp_bef; // pions' momentum from MC truth   
  double m_pionm_bef;   
  double m_pionp_pt_bef; // pions' momentum from MC truth   
  double m_pionm_pt_bef;   
  
  // with the method of momentum selection

  // momentum
  double m_mom_lepm;		    
  double m_mom_lepp;	    
  double m_mom_pionp;		    
  double m_mom_pionm;		    
  double m_pipi_dang;		    
  long m_cms_index; //in p; theta; phi sequence
  double m_cms_lepp[500];		    
  double m_cms_lepm[500];
  double m_mass_twopi;       
  double m_mass_jpsi;   
  double m_mass_recoil;   
  double m_inv_mass;   
  double m_tot_e;   
  double m_tot_px;   
  double m_tot_py;   
  double m_tot_pz;   
  double m_emc_ene;
  long m_event_flag; 
  // 3 or 4 or 5 tracks, 
  // 4=>4 tracks, 
  // 0=> miss pi+, 1=> miss pi-, 2=> miss lepton+, 3=> miss lepton-
  // 5=> more pi+, 6=> more pi-, 7=> more lepton+, 8=> more lepton-
  double m_pt_lepm;		    
  double m_pt_lepp;	    
  double m_pt_pionp;		    
  double m_pt_pionm;	    
  long m_run;
  long m_event;
  long m_index;
  double m_cos_theta[500]; // in pi+ pi- l+ l- sequence
  double m_ep_ratio[500]; 
  double m_phi[500];
  double m_prob_e[500];
  double m_prob_mu[500];
  double m_prob_pi[500];
  double m_dedx_e[500];
  double m_dedx_mu[500];
  double m_dedx_pi[500];
  double m_tof_e[500];
  double m_tof_mu[500];
  double m_tof_pi[500];
  double m_dedx_ghits[500];
  double m_muc_dep[500];
  double m_four_mom[4][5];
  double m_recoil_four_mom[4][5]; // each track from other three
  double m_miss_mass[500]; // the invariant mass of missing track
  double m_how_near[500]; // match and resolution
  double m_vertex[4][5]; // 4x2 matrix, (pi+, pi-, l+, l-)x(vr, vz)  

  // veto pi0 to veto pi0 pi0 J/psi, eta J/psi, all photon cuts
  double m_num_pi0;
  double m_like_pi0_mass;
  double m_pipipi0_mass;
  double m_omegajpsi_mass;
  double m_gammajpsi_mass;//combined with the most energic photon
  double f_gammajpsi_mass;// after 4C fit
  double f_chisq_gammajpsi;// after 4C fit
  double m_max_gam_ene;
  double m_max_gam_costhe;
  double m_max_gam_phi;
  double m_max_gam_cms_ene;
  double m_max_gam_dang;
  double m_like_eta_mass;

  // 1C fit of J/psi
  double f_chisq_vt;
  double f_chisq_jpsi;
  double f_mom_lepm;		    
  double f_mom_lepp;	    
  double f_mom_pionp;		    
  double f_mom_pionm;		    
  double f_cms_lepp;		    
  double f_cms_lepm;
  double f_pipi_dang;		    
  double f_mass_twopi;       
  double f_mass_jpsi;   
  double f_mass_recoil;   
  double f_inv_mass;   
  double f_tot_e;   
  double f_tot_px;   
  double f_tot_py;   
  double f_tot_pz;   
  double f_pt_lepm;		    
  double f_pt_lepp;	    
  double f_pt_pionp;		    
  double f_pt_pionm;	    
  long f_index;
  double f_cos_theta[500]; // in pi+ pi- l+ l- sequence
  double f_phi[500];
  double f_four_mom[4][5];
  double f_recoil_four_mom[4][5]; // each track from other three
  double f_miss_mass[500];
  double f_how_near[500]; // match and resolution

  // we don't need smear, actually, assume one track of pi pi is missing,
  // the recoil mass is actually the invariant of pair-lepton
  // following two lines are applied to check MDC and EMC match
  long m_pion_matched;
  long m_lep_matched;
  // flag from EventHeader
  long m_Nch;
  long m_Nty;
  // book MCtruth
  long   m_idxmc;
  //long   m_run;
  long  m_pdgid[500];
  long  m_motheridx[500];
  long    m_xyzid;
  long    m_true_cms_index;
  double m_true_cms_lepp[500];		    
  double m_true_cms_lepm[500];

  // trigger and t0
   long  m_trig_index;
   long  m_trig_cond[500];
   long m_trig_chan[500];
   long m_trig_timewindow;
   long m_trig_timetype;

   double  m_est_start;
   long  m_est_status;
   double  m_est_quality;

};


#endif 
