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

class PipiJpsi : public Algorithm {

public:
  PipiJpsi(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();  

private:
  void f_get_pid(EvtRecTrack *recTrk, NTuple::Array<double> m_prob_e, NTuple::Array<double> m_prob_mu, NTuple::Array<double> m_prob_pi, int n);
  void f_get_dedx(EvtRecTrack *recTrk, NTuple::Array<double> m_dedx_e, NTuple::Array<double> m_dedx_mu, NTuple::Array<double> m_dedx_pi, NTuple::Array<double> m_dedx_ghits, int n);
  void f_get_tof(EvtRecTrack *recTrk, NTuple::Array<double> m_tof_e, NTuple::Array<double> m_tof_mu, NTuple::Array<double> m_tof_pi, int n);

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
  IHistogram1D *h_vr0;
  IHistogram1D *h_vz0;

  // define Ntuples here
  NTuple::Tuple*  m_tuple2;      // fake photon
  NTuple::Item<double>  m_dthe;
  NTuple::Item<double>  m_dphi;
  NTuple::Item<double>  m_dang;
  NTuple::Item<double>  m_eraw;
  NTuple::Item<long>  m_nGam;

  NTuple::Tuple* m_tuple3;    // dE/dx
  NTuple::Item<double> m_ptrk;
  NTuple::Item<double> m_chie;
  NTuple::Item<double> m_chimu;
  NTuple::Item<double> m_chipi;
  NTuple::Item<double> m_chik;
  NTuple::Item<double> m_chip;
  NTuple::Item<double> m_probPH;
  NTuple::Item<double> m_normPH;
  NTuple::Item<double> m_ghit;
  NTuple::Item<double> m_thit;

  NTuple::Tuple* m_tuple4;   // endcap tof
  NTuple::Item<double> m_ptot_etof;
  NTuple::Item<double> m_cntr_etof;
  NTuple::Item<double> m_path_etof;
  NTuple::Item<double> m_tof_etof;
  NTuple::Item<double> m_te_etof;
  NTuple::Item<double> m_tmu_etof;
  NTuple::Item<double> m_tpi_etof;
  NTuple::Item<double> m_tk_etof;
  NTuple::Item<double> m_tp_etof;
  NTuple::Item<double> m_ph_etof;
  NTuple::Item<double> m_rhit_etof;
  NTuple::Item<double> m_qual_etof;

  NTuple::Tuple* m_tuple5;  // barrel inner tof
  NTuple::Item<double> m_ptot_btof1;
  NTuple::Item<double> m_cntr_btof1;
  NTuple::Item<double> m_path_btof1;
  NTuple::Item<double> m_tof_btof1;
  NTuple::Item<double> m_te_btof1;
  NTuple::Item<double> m_tmu_btof1;
  NTuple::Item<double> m_tpi_btof1;
  NTuple::Item<double> m_tk_btof1;
  NTuple::Item<double> m_tp_btof1;
  NTuple::Item<double> m_ph_btof1;
  NTuple::Item<double> m_zhit_btof1;
  NTuple::Item<double> m_qual_btof1;

  NTuple::Tuple* m_tuple6;  // barrel outer tof
  NTuple::Item<double> m_ptot_btof2;
  NTuple::Item<double> m_cntr_btof2;
  NTuple::Item<double> m_path_btof2;
  NTuple::Item<double> m_tof_btof2;
  NTuple::Item<double> m_te_btof2;
  NTuple::Item<double> m_tmu_btof2;
  NTuple::Item<double> m_tpi_btof2;
  NTuple::Item<double> m_tk_btof2;
  NTuple::Item<double> m_tp_btof2;
  NTuple::Item<double> m_ph_btof2;
  NTuple::Item<double> m_zhit_btof2;
  NTuple::Item<double> m_qual_btof2;

  // book mctruth before any selection to check the generator
  NTuple::Tuple *m_tuple7;
  NTuple::Item<double> m_pionp_bef; // pions' momentum from MC truth   
  NTuple::Item<double> m_pionm_bef;   
  NTuple::Item<double> m_pionp_pt_bef; // pions' momentum from MC truth   
  NTuple::Item<double> m_pionm_pt_bef;   
  
  // with the method of momentum selection
  NTuple::Tuple* m_tuple8; 
  // momentum
  NTuple::Item<double> m_mom_lepm;		    
  NTuple::Item<double> m_mom_lepp;	    
  NTuple::Item<double> m_mom_pionp;		    
  NTuple::Item<double> m_mom_pionm;		    
  NTuple::Item<double> m_pipi_dang;		    
  NTuple::Item<long>   m_cms_index; //in p; theta; phi sequence
  NTuple::Array<double> m_cms_lepp;		    
  NTuple::Array<double> m_cms_lepm;
  NTuple::Item<double> m_mass_twopi;       
  NTuple::Item<double> m_mass_jpsi;   
  NTuple::Item<double> m_mass_recoil;   
  NTuple::Item<double> m_inv_mass;   
  NTuple::Item<double> m_tot_e;   
  NTuple::Item<double> m_tot_px;   
  NTuple::Item<double> m_tot_py;   
  NTuple::Item<double> m_tot_pz;   
  NTuple::Item<double> m_emc_ene;
  NTuple::Item<long>   m_event_flag; 
  // 3 or 4 or 5 tracks, 
  // 4=>4 tracks, 
  // 0=> miss pi+, 1=> miss pi-, 2=> miss lepton+, 3=> miss lepton-
  // 5=> more pi+, 6=> more pi-, 7=> more lepton+, 8=> more lepton-
  NTuple::Item<double> m_pt_lepm;		    
  NTuple::Item<double> m_pt_lepp;	    
  NTuple::Item<double> m_pt_pionp;		    
  NTuple::Item<double> m_pt_pionm;	    
  NTuple::Item<long> m_run;
  NTuple::Item<long> m_event;
  NTuple::Item<long> m_index;
  NTuple::Array<double> m_cos_theta; // in pi+ pi- l+ l- sequence
  NTuple::Array<double> m_ep_ratio; 
  NTuple::Array<double> m_phi;
  NTuple::Array<double> m_prob_e;
  NTuple::Array<double> m_prob_mu;
  NTuple::Array<double> m_prob_pi;
  NTuple::Array<double> m_dedx_e;
  NTuple::Array<double> m_dedx_mu;
  NTuple::Array<double> m_dedx_pi;
  NTuple::Array<double> m_tof_e;
  NTuple::Array<double> m_tof_mu;
  NTuple::Array<double> m_tof_pi;
  NTuple::Array<double> m_dedx_ghits;
  NTuple::Array<double> m_muc_dep;
  NTuple::Matrix<double> m_four_mom;
  NTuple::Matrix<double> m_recoil_four_mom; // each track from other three
  NTuple::Array<double> m_miss_mass; // the invariant mass of missing track
  NTuple::Array<double> m_how_near; // match and resolution
  NTuple::Matrix<double> m_vertex; // 4x2 matrix, (pi+, pi-, l+, l-)x(vr, vz)  

  // veto pi0 to veto pi0 pi0 J/psi, eta J/psi, all photon cuts
  NTuple::Item<double> m_num_pi0;
  NTuple::Item<double> m_like_pi0_mass;
  NTuple::Item<double> m_pipipi0_mass;
  NTuple::Item<double> m_omegajpsi_mass;
  NTuple::Item<double> m_gammajpsi_mass;//combined with the most energic photon
  NTuple::Item<double> f_gammajpsi_mass;// after 4C fit
  NTuple::Item<double> f_chisq_gammajpsi;// after 4C fit
  NTuple::Item<double> m_max_gam_ene;
  NTuple::Item<double> m_max_gam_costhe;
  NTuple::Item<double> m_max_gam_phi;
  NTuple::Item<double> m_max_gam_cms_ene;
  NTuple::Item<double> m_max_gam_dang;
  NTuple::Item<double> m_like_eta_mass;

  // 1C fit of J/psi
  NTuple::Item<double> f_chisq_vt;
  NTuple::Item<double> f_chisq_jpsi;
  NTuple::Item<double> f_mom_lepm;		    
  NTuple::Item<double> f_mom_lepp;	    
  NTuple::Item<double> f_mom_pionp;		    
  NTuple::Item<double> f_mom_pionm;		    
  NTuple::Item<double> f_cms_lepp;		    
  NTuple::Item<double> f_cms_lepm;
  NTuple::Item<double> f_pipi_dang;		    
  NTuple::Item<double> f_mass_twopi;       
  NTuple::Item<double> f_mass_jpsi;   
  NTuple::Item<double> f_mass_recoil;   
  NTuple::Item<double> f_inv_mass;   
  NTuple::Item<double> f_tot_e;   
  NTuple::Item<double> f_tot_px;   
  NTuple::Item<double> f_tot_py;   
  NTuple::Item<double> f_tot_pz;   
  NTuple::Item<double> f_pt_lepm;		    
  NTuple::Item<double> f_pt_lepp;	    
  NTuple::Item<double> f_pt_pionp;		    
  NTuple::Item<double> f_pt_pionm;	    
  NTuple::Item<long> f_index;
  NTuple::Array<double> f_cos_theta; // in pi+ pi- l+ l- sequence
  NTuple::Array<double> f_phi;
  NTuple::Matrix<double> f_four_mom;
  NTuple::Matrix<double> f_recoil_four_mom; // each track from other three
  NTuple::Array<double> f_miss_mass;
  NTuple::Array<double> f_how_near; // match and resolution

  // we don't need smear, actually, assume one track of pi pi is missing,
  // the recoil mass is actually the invariant of pair-lepton
  // following two lines are applied to check MDC and EMC match
  NTuple::Item<long> m_pion_matched;
  NTuple::Item<long> m_lep_matched;
  // flag from EventHeader
  NTuple::Item<long> m_Nch;
  NTuple::Item<long> m_Nty;
  // book MCtruth
  NTuple::Item<long>   m_idxmc;
  NTuple::Array<long>  m_pdgid;
  NTuple::Array<long>  m_motheridx;
  NTuple::Item<long>    m_xyzid;
  NTuple::Item<long>    m_true_cms_index;
  NTuple::Array<double> m_true_cms_lepp;		    
  NTuple::Array<double> m_true_cms_lepm;

  // trigger and t0
   NTuple::Item<long>  m_trig_index;
   NTuple::Array<long>  m_trig_cond;
   NTuple::Array<long> m_trig_chan;
   NTuple::Item<long> m_trig_timewindow;
   NTuple::Item<long> m_trig_timetype;

   NTuple::Item<double>  m_est_start;
   NTuple::Item<long>  m_est_status;
   NTuple::Item<double>  m_est_quality;

};


#endif 
