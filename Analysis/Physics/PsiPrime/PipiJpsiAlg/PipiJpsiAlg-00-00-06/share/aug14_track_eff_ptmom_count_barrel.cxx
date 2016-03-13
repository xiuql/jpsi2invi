#include "TCut.h"
#include "cut_pipijpsi.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TText.h"
#include "TGraphErrors.h"
#include <iostream>


TCut m_cut_three_yes = "invjpsi>3.0 && invjpsi<3.2 && pipidang<0.98  && cmslepm<1.7 && cmslepp<1.7 && cmslepm>1.4 && cmslepp>1.4 && mompionm<0.5 && mompionp<0.5 && Npi0<1 && maxene<0.3 && epratio<0.9 && (((eveflag==1||eveflag==6) && missmass<0.25) || eveflag==4)";
TCut m_cut_four_yes = "invjpsi>3.0 && invjpsi<3.2 && pipidang<0.98  && cmslepm<1.7 && cmslepp<1.7 && cmslepm>1.4 && cmslepp>1.4 && mompionm<0.5 && mompionp<0.5 && Npi0<1 && maxene<0.3 && epratio<0.9 &&  eveflag==4";

TCut m_cut_three_no = "invjpsi>3.0 && invjpsi<3.2 && pipidang<0.98  && cmslepm<1.7 && cmslepp<1.7 && cmslepm>1.4 && cmslepp>1.4 && mompionm<0.5 && mompionp<0.5 && epratio<1.0 && (eveflag==1||eveflag==6 || eveflag==4)";
TCut m_cut_four_no = "invjpsi>3.0 && invjpsi<3.2 && pipidang<0.98  && cmslepm<1.7 && cmslepp<1.7 && cmslepm>1.4 && cmslepp>1.4 && mompionm<0.5 && mompionp<0.5 &&  epratio<1.0 &&  eveflag==4";

void f_eff(TTree *m_tree, TCut &cut_all, TCut &cut_four, Double_t *m_eff, Double_t *m_err, Double_t *m_xcen, Double_t *m_xerr){
  TH1F *hpionm = new TH1F("hpionm","hpionm",7,0.05,0.4);
  TH1F *hfourm = new TH1F("hfourm","hfourm",7,0.05,0.4);
  m_tree->Draw("tppionratiom*mompionm>>hpionm",cut_all+"fabs(costhe[1])<0.7"); 
  m_tree->Draw("tppionratiom*mompionm>>hfourm",cut_four+"fabs(costhe[1])<0.7"); 

  for(Int_t i=1; i<8; i++){
    m_eff[i-1] = hpionm->GetBinContent(i) ? hfourm->GetBinContent(i)/hpionm->GetBinContent(i) : 0;
    m_err[i-1] = m_eff[i-1]*(hpionm->GetBinContent(i)) ? sqrt((m_eff[i-1])*(1-m_eff[i-1])/((hpionm->GetBinContent(i)))) : 0 ;
    m_xcen[i-1] = hfourm->GetBinCenter(i);
    m_xerr[i-1] = hfourm->GetBinWidth(i)/2;
    cout << "pt mom: " << 0.05+(i-1)*(0.4-0.05)/7 << "   eff: " <<m_eff[i-1] << "   err: " << m_err[i-1] << endl; 
  }
}

//void aug14_track_eff_ptmom_count_barrel(TString file_fir_name="../../651_psip_pipijpsi_data_check.root", TString file_sec_name="../../651_psip_pipijpsi_dst_uu.root", TCut m_cut_run="", TCut m_cut_run="", TString tree_name="infmom", Int_t m_bins=7){
void aug14_track_eff_ptmom_count_barrel(TString file_fir_name="../../651_pipijpsi_inclusive.root", TString file_sec_name="../../651_pipijpsi_inclusive.root", TCut m_cut_run="", TCut m_cut_run="", TString tree_name="infmom", Int_t m_bins=7){

  gROOT->SetStyle("BES");
  TFile *m_file = new TFile(file_fir_name);
  TTree *m_tree = (TTree*)m_file->Get(tree_name);

  Double_t m_data_eff[7], m_data_err[7], m_mc_eff[7], m_mc_err[7];
  Double_t m_xcen[7], m_xerr[7];
  f_eff(m_tree, m_cut_three_yes, m_cut_four_yes, m_data_eff, m_data_err, m_xcen, m_xerr);

  TFile *m_file2 = new TFile(file_sec_name);
  TTree *m_tree2 = (TTree*)m_file2->Get(tree_name);
  f_eff(m_tree2, m_cut_three_no, m_cut_four_no, m_mc_eff, m_mc_err, m_xcen, m_xerr);

  // create the TGraphErrors and draw it
  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,643);
  TGraphErrors *gr_data = new TGraphErrors(m_bins, m_xcen, m_data_eff, m_xerr, m_data_err);
  gr_data->SetTitle("Tracking effciency of #pi^{-}");
  gr_data->GetXaxis()->SetTitle("transverse momentum (GeV)");
  gr_data->GetYaxis()->SetTitle("efficiency");
  gr_data->GetYaxis()->SetRangeUser(0.65,1.0);
  gr_data->SetMarkerColor(2);
  gr_data->SetMarkerStyle(20);

  TGraphErrors *gr_mc = new TGraphErrors(m_bins, m_xcen, m_mc_eff, m_xerr, m_mc_err);
  gr_mc->SetMarkerColor(4);
  gr_mc->SetMarkerStyle(21);
  gr_mc->GetXaxis()->SetTitle("transverse momentum (GeV)");
  gr_mc->GetYaxis()->SetTitle("efficiency");
  gr_data->Draw("AP");
  gr_mc->Draw("P");
  c1->Update();

  TText m_text;
  m_text.SetTextSize(0.08);
  m_text.SetTextColor(2);
  m_text.DrawTextNDC(0.5,0.4,"bg 2.4%");

  m_text.SetTextColor(4);
  m_text.DrawTextNDC(0.5,0.3,"bg 12%");

  m_text.SetTextColor(1);
  m_text.DrawTextNDC(0.5,0.2,"651, Barrel");
  m_text.DrawTextNDC(0.5,0.5,"Inclusive MC");

  c1->Print("figs/nov30_track_eff_ptmom_count_barrel_651_bg_study.eps");
  c1->Print("figs/nov30_track_eff_ptmom_count_barrel_651_bg_study.pdf");
}

