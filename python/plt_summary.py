#!/usr/bin/env python
"""
Plot summary histograms 
"""

__author__ = "SHI Xin <shixin@ihep.ac.cn>"
__copyright__ = "Copyright (c) SHI Xin"
__created__ = "[2016-07-25 Mon 09:22]" 

import os
import sys 
import ROOT 
from tools import check_outfile_path



def main():
    set_root_style 
    f1 = ROOT.TFile('run/hist/jpsi2invi_data_merged_1.root')
    
    h1_mrecpipi = f1.Get('h_mrecpipi')
    h1_mrecpipi.Sumw2()
    h1_mrecpipi.SetXTitle('M(recoil(#pi^{+}#pi^{-})) (GeV/c^{2})') 
    h1_mrecpipi.SetYTitle('Events/(0.0012 GeV/c^{2})') 

    h1_mpipi = f1.Get('h_mpipi')
    h1_mpipi.Sumw2()
    h1_mpipi.SetXTitle('M(#pi^{+}#pi^{-}) (GeV/c^{2})') 
    h1_mpipi.SetYTitle('Events/(0.0012 GeV/c^{2})') 

    h1_pip_p = f1.Get('h_pip_p')
    h1_pip_p.Sumw2()
    h1_pip_p.SetXTitle('P(#pi^{+}) (GeV/c)')
    h1_pip_p.SetYTitle('Events/(0.005 GeV/c)') 
    
    h1_pim_p = f1.Get('h_pim_p')
    h1_pim_p.Sumw2()
    h1_pim_p.SetXTitle('P(#pi^{-}) (GeV/c)')
    h1_pim_p.SetYTitle('Events/(0.005 GeV/c)') 

    h1_pip_costhe = f1.Get('h_pip_costhe')
    h1_pip_costhe.SetXTitle('Cos#theta_{#pi^{+}}')
    h1_pip_costhe.SetYTitle('Events/0.02')

    h1_pim_costhe = f1.Get('h_pim_costhe')
    h1_pim_costhe.SetXTitle('Cos#theta_{#pi^{-}}')
    h1_pim_costhe.SetYTitle('Events/0.02')

    h1_cospipi = f1.Get('h_cospipi')
    h1_cospipi.SetXTitle('Cos#theta_{#pi^{+}#pi^{-}}')
    h1_cospipi.SetYTitle('Events/0.01')

    h1_cos2pisys = f1.Get('h_cos2pisys')
    h1_cos2pisys.SetXTitle('Cos#theta_{#pi#pi sys.}')
    h1_cos2pisys.SetYTitle('Events/0.02')

    h1_ngam = f1.Get('h_ngam')
    h1_ngam.SetXTitle('N_{#gamma}')
    h1_ngam.SetYTitle('Events')
    
    # c = ROOT.TCanvas('c', 'c', 800, 800)
    # h1_mrecpipi.Draw()
    # figfile = 'doc/fig/jpsi2invi_data_mrecpipi.pdf'
    # check_outfile_path(figfile)
    # c.SaveAs(figfile)

    outfile = 'run/summary/jpsi2invi_data.root'
    check_outfile_path(outfile)
    
    fout = ROOT.TFile(outfile, 'RECREATE')
    h1_mrecpipi.Write()
    h1_mpipi.Write()
    h1_pip_p.Write()
    h1_pim_p.Write()
    h1_pip_costhe.Write()
    h1_pim_costhe.Write()
    h1_cospipi.Write()
    h1_cos2pisys.Write()
    h1_ngam.Write()
    fout.Close()

    
def set_root_style(stat=0, grid=0): 
  ROOT.gROOT.Reset()

  ROOT.gStyle.SetTitleFillColor(0)   
  ROOT.gStyle.SetTitleBorderSize(0)  
    
  ROOT.gStyle.SetCanvasBorderMode(0) 
  ROOT.gStyle.SetCanvasColor(0) 
  ROOT.gStyle.SetCanvasDefX(0)  
  ROOT.gStyle.SetCanvasDefY(0)  
  ROOT.gStyle.SetFrameBorderMode(0)  
  ROOT.gStyle.SetFrameBorderSize(1)  
  ROOT.gStyle.SetFrameFillColor(0)  
  ROOT.gStyle.SetFrameFillStyle(0)  
  ROOT.gStyle.SetFrameLineColor(1)  
  ROOT.gStyle.SetFrameLineStyle(1)  
  ROOT.gStyle.SetFrameLineWidth(1)  

  #ROOT.gStyle.SetPadTopMargin(PadTopMargin)   
  ROOT.gStyle.SetPadLeftMargin(0.10)   
  ROOT.gStyle.SetPadRightMargin(0.25)   

  ROOT.gStyle.SetLabelSize(0.03, 'XYZ')   
  ROOT.gStyle.SetTitleSize(0.04, 'XYZ')   
  ROOT.gStyle.SetTitleOffset(1.2, 'Y')   

  ROOT.gStyle.SetPadBorderMode(0)   
  ROOT.gStyle.SetPadColor(0)   
  ROOT.gStyle.SetPadTickX(1)  
  ROOT.gStyle.SetPadTickY(1)  
  ROOT.gStyle.SetPadGridX(grid)  
  ROOT.gStyle.SetPadGridY(grid)  

  ROOT.gStyle.SetOptStat(stat)  
  ROOT.gStyle.SetStatColor(0)  
  ROOT.gStyle.SetStatBorderSize(1)

  
    
if __name__ == '__main__':
    main()
