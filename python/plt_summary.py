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
from tools import check_outfile_path, set_root_style



def main():
    set_root_style(stat=0, grid=0) 
    ROOT.gStyle.SetPadLeftMargin(0.15)

    c = ROOT.TCanvas('c', 'c', 800, 800)    

    f1 = ROOT.TFile('run/hist/jpsi2invi_data_merged_1.root')

    h1_mrecpipi = draw_mrecpipi(c, f1)
    h1_mpipi = draw_mpipi(c, f1)
    h1_pip_p = draw_pip_p(c, f1) 
    h1_pim_p = draw_pim_p(c, f1) 
    h1_pip_costhe = draw_pip_costhe(c, f1) 
    h1_pim_costhe = draw_pim_costhe(c, f1)
    h1_cospipi = draw_cospipi(c, f1) 
    h1_cos2pisys = draw_cos2pisys(c, f1)
    h1_ngam = draw_ngam(c, f1) 
    
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


def draw_mrecpipi(c, f1):
    h1 = f1.Get('h_mrecpipi')
    h1.Sumw2()
    h1.SetXTitle('M(recoil(#pi^{+}#pi^{-})) (GeV/c^{2})') 
    h1.SetYTitle('Events/(0.0012 GeV/c^{2})')
    h1.GetXaxis().SetLabelSize(0.03) 
    h1.GetYaxis().SetLabelSize(0.03) 
    h1.GetYaxis().SetTitleOffset(1.8) 
    h1.SetMarkerStyle(ROOT.kFullDotLarge)

    h1.Draw()
    figfile = 'doc/fig/jpsi2invi_data_mrecpipi.pdf'
    check_outfile_path(figfile)
    c.SaveAs(figfile)
    return h1


def draw_mpipi(c, f1):
    h1 = f1.Get('h_mpipi')
    h1.Sumw2()
    h1.SetXTitle('M(#pi^{+}#pi^{-}) (GeV/c^{2})') 
    h1.SetYTitle('Events/(0.0012 GeV/c^{2})')
    h1.GetYaxis().SetTitleOffset(1.8) 
    h1.GetXaxis().SetLabelSize(0.03) 
    h1.GetYaxis().SetLabelSize(0.03) 
    h1.SetMarkerStyle(ROOT.kFullDotLarge)

    h1.Draw()
    figfile = 'doc/fig/jpsi2invi_data_mpipi.pdf'
    c.SaveAs(figfile)
    return h1


def draw_pip_p(c, f1):
    h1 = f1.Get('h_pip_p')
    h1.Sumw2()
    h1.SetXTitle('P(#pi^{+}) (GeV/c)')
    h1.SetYTitle('Events/(0.005 GeV/c)') 
    h1.GetYaxis().SetTitleOffset(1.8) 
    h1.GetXaxis().SetLabelSize(0.03) 
    h1.GetYaxis().SetLabelSize(0.03) 
    h1.SetMarkerStyle(ROOT.kFullDotLarge)

    h1.Draw()
    figfile = 'doc/fig/jpsi2invi_data_pip_p.pdf'
    c.SaveAs(figfile)
    return h1

def draw_pim_p(c, f1):
    h1 = f1.Get('h_pim_p')
    h1.Sumw2()
    h1.SetXTitle('P(#pi^{-}) (GeV/c)')
    h1.SetYTitle('Events/(0.005 GeV/c)') 
    h1.GetYaxis().SetTitleOffset(1.8) 
    h1.GetXaxis().SetLabelSize(0.03) 
    h1.GetYaxis().SetLabelSize(0.03) 
    h1.SetMarkerStyle(ROOT.kFullDotLarge)

    h1.Draw()
    figfile = 'doc/fig/jpsi2invi_data_pim_p.pdf'
    c.SaveAs(figfile)
    return h1


def draw_pip_costhe(c, f1):
    h1 = f1.Get('h_pip_costhe')
    h1.Sumw2()
    h1.SetXTitle('Cos#theta_{#pi^{+}}')
    h1.SetYTitle('Events/0.02')
    h1.GetYaxis().SetTitleOffset(1.8) 
    h1.GetXaxis().SetLabelSize(0.03) 
    h1.GetYaxis().SetLabelSize(0.03) 
    h1.SetMarkerStyle(ROOT.kFullDotLarge)

    h1.Draw()
    figfile = 'doc/fig/jpsi2invi_data_pip_costhe.pdf'
    c.SaveAs(figfile)
    return h1

def draw_pim_costhe(c, f1):
    h1 = f1.Get('h_pim_costhe')
    h1.Sumw2()
    h1.SetXTitle('Cos#theta_{#pi^{-}}')
    h1.SetYTitle('Events/0.02')
    h1.GetYaxis().SetTitleOffset(1.8) 
    h1.GetXaxis().SetLabelSize(0.03) 
    h1.GetYaxis().SetLabelSize(0.03) 
    h1.SetMarkerStyle(ROOT.kFullDotLarge)

    h1.Draw()
    figfile = 'doc/fig/jpsi2invi_data_pim_costhe.pdf'
    c.SaveAs(figfile)
    return h1


def draw_cospipi(c, f1):
    h1 = f1.Get('h_cospipi')
    h1.Sumw2()
    h1.SetXTitle('Cos#theta_{#pi^{+}#pi^{-}}')
    h1.SetYTitle('Events/0.01')
    h1.GetYaxis().SetTitleOffset(1.8) 
    h1.GetXaxis().SetLabelSize(0.03) 
    h1.GetYaxis().SetLabelSize(0.03) 
    h1.SetMarkerStyle(ROOT.kFullDotLarge)

    h1.Draw()
    figfile = 'doc/fig/jpsi2invi_data_cospipi.pdf'
    c.SaveAs(figfile)
    return h1


def draw_cos2pisys(c, f1):
    h1 = f1.Get('h_cos2pisys')
    h1.Sumw2()
    h1.SetXTitle('Cos#theta_{#pi#pi sys.}')
    h1.SetYTitle('Events/0.02')
    h1.GetYaxis().SetTitleOffset(1.8) 
    h1.GetXaxis().SetLabelSize(0.03) 
    h1.GetYaxis().SetLabelSize(0.03) 
    h1.SetMarkerStyle(ROOT.kFullDotLarge)

    h1.Draw()
    figfile = 'doc/fig/jpsi2invi_data_cos2pisys.pdf'
    c.SaveAs(figfile)
    return h1

def draw_ngam(c, f1):
    h1 = f1.Get('h_ngam')
    h1.Sumw2()
    h1.SetXTitle('N_{#gamma}')
    h1.SetYTitle('Events')
    h1.GetYaxis().SetTitleOffset(1.8) 
    h1.GetXaxis().SetLabelSize(0.03) 
    h1.GetYaxis().SetLabelSize(0.03) 
    h1.SetMarkerStyle(ROOT.kFullDotLarge)

    h1.Draw()
    figfile = 'doc/fig/jpsi2invi_data_ngam.pdf'
    c.SaveAs(figfile)
    return h1

    
if __name__ == '__main__':
    main()
