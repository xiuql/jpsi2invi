#!/usr/bin/env python
"""
Event Selection 
"""

__author__ = "SHI Xin <shixin@ihep.ac.cn>"
__copyright__ = "Copyright (c) SHI Xin"
__created__ = "[2016-06-28 Tue 09:09]" 

import sys
import os
import math 
import ROOT 
from progressbar import Bar, Percentage, ProgressBar
from time import time 
from tools import duration, check_outfile_path

#TEST=True 
TEST=False

# Global constants 
JPSI_MASS = 3.096916; 

# Global histograms

h_evtflw = ROOT.TH1F('hevtflw', 'eventflow', 10, 0, 10) 
h_evtflw.GetXaxis().SetBinLabel(1, 'raw')
h_evtflw.GetXaxis().SetBinLabel(2, 'N_{#gamma}=0')
h_evtflw.GetXaxis().SetBinLabel(3, '|cos#theta|<0.8')
h_evtflw.GetXaxis().SetBinLabel(4, '|p|<0.45') 
h_evtflw.GetXaxis().SetBinLabel(5, 'PID') 
h_evtflw.GetXaxis().SetBinLabel(6, 'cos#theta_{#pi^{+}#pi^{-}}<0.95') 
h_evtflw.GetXaxis().SetBinLabel(7, 'cos#theta_{#pi#pi sys}<0.9') 
h_evtflw.GetXaxis().SetBinLabel(8, '3<M_{#pi#pi}^{rec}<3.2') 

h_mrecpipi = ROOT.TH1D('h_mrecpipi', 'mrecpipi', 100, 3.04, 3.16)
h_mpipi = ROOT.TH1D('h_mpipi', 'mpipi', 100, 0.2, 0.7) 
h_pip_p = ROOT.TH1D('h_pip_p', 'pip_p', 100, 0.0, 0.5) 
h_pim_p = ROOT.TH1D('h_pim_p', 'pim_p', 100, 0.0, 0.5) 
h_pip_costhe = ROOT.TH1D('h_pip_costhe', 'pip_costhe', 100, -1.0, 1.0)
h_pim_costhe = ROOT.TH1D('h_pim_costhe', 'pim_costhe', 100, -1.0, 1.0)
h_cospipi = ROOT.TH1D('h_cospipi', 'cospipi', 100, 0.0, 1.0)
h_cos2pisys = ROOT.TH1D('h_cos2pisys', 'cos2pisys', 100, -1.0, 1.0)
h_ngam = ROOT.TH1D('h_ngam', 'ngam', 100, 0, 20)


def usage():
    sys.stdout.write('''
NAME
    sel_events.py 

SYNOPSIS

    ./sel_events.py infile outfile 

AUTHOR 
    SHI Xin <shixin@ihep.ac.cn> 

DATE
    July 2016 
\n''')

    
def main():

    if TEST:
        sys.stdout.write('Run in TEST mode! \n')
    
    args = sys.argv[1:]
    if len(args) < 2:
        return usage()
    
    infile = args[0]
    outfile = args[1]
    check_outfile_path(outfile)

    fin = ROOT.TFile(infile)
    t = fin.Get('tree')
    entries = t.GetEntriesFast()

    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=entries).start()
    time_start = time()

    for jentry in xrange(entries):
        pbar.update(jentry+1)
        # get the next tree in the chain and verify
        ientry = t.LoadTree(jentry)
        if ientry < 0:
            break
        # copy next entry into memory and verify

        if TEST and ientry > 10000:
            break
        
        nb = t.GetEntry(jentry)
        if nb<=0:
            continue

        fill_histograms(t)
        
        if select_jpsi_to_invisible(t): 
            h_mrecpipi.Fill(t.vtx_mrecpipi)
 
    fout = ROOT.TFile(outfile, "RECREATE")
    write_histograms() 
    fout.Close()
    pbar.finish()
    
    dur = duration(time()-time_start)
    sys.stdout.write(' \nDone in %s. \n' % dur) 


def fill_histograms(t):
    cut_ngam = (t.ngam == 0)
    cut_trkp_costhe = (abs(math.cos(t.trkp_theta)) < 0.8)
    cut_trkm_costhe = (abs(math.cos(t.trkm_theta)) < 0.8)
    cut_trkp_p = (abs(t.trkp_p) < 0.45) 
    cut_trkm_p = (abs(t.trkm_p) < 0.45)
    cut_cospipi =  (t.vtx_cospipi < 0.95)
    cut_cos2pisys = (t.vtx_cos2pisys < 0.9)
    cut_pi_PID = (t.prob_pip > t.prob_kp and t.prob_pip > 0.001 and
                  t.prob_pim > t.prob_km and t.prob_pim > 0.001)
    cut_mjpsi_win = (t.vtx_mrecpipi > 3.0 and t.vtx_mrecpipi < 3.2)
    cut_mjpsi_sig = (abs(t.vtx_mrecpipi - JPSI_MASS)<0.015)

    if (cut_ngam and cut_trkp_costhe and cut_trkm_costhe and cut_trkp_p and cut_trkm_p and
        cut_cospipi and cut_cos2pisys and cut_pi_PID and cut_mjpsi_sig):
        h_mpipi.Fill(t.vtx_mpipi)

    if (cut_ngam and cut_trkp_costhe and cut_trkm_costhe                and cut_trkm_p and
        cut_cospipi and cut_cos2pisys and cut_pi_PID and cut_mjpsi_sig):
        h_pip_p.Fill(t.trkp_p)

    if (cut_ngam and cut_trkp_costhe and cut_trkm_costhe and cut_trkp_p                and
        cut_cospipi and cut_cos2pisys and cut_pi_PID and cut_mjpsi_sig):
        h_pim_p.Fill(t.trkm_p)

    if (cut_ngam and                     cut_trkm_costhe and cut_trkp_p and cut_trkm_p and
        cut_cospipi and cut_cos2pisys and cut_pi_PID and cut_mjpsi_sig):
        h_pip_costhe.Fill(math.cos(t.trkp_theta))

    if (cut_ngam and cut_trkp_costhe                     and cut_trkp_p and cut_trkm_p and
        cut_cospipi and cut_cos2pisys and cut_pi_PID and cut_mjpsi_sig):
        h_pim_costhe.Fill(math.cos(t.trkm_theta))

        
    if (cut_ngam and cut_trkp_costhe and cut_trkm_costhe and cut_trkp_p and cut_trkm_p and
                        cut_cos2pisys and cut_pi_PID and cut_mjpsi_sig):
        h_cospipi.Fill(t.vtx_cospipi)

    if (cut_ngam and cut_trkp_costhe and cut_trkm_costhe and cut_trkp_p and cut_trkm_p and
        cut_cospipi                   and cut_pi_PID and cut_mjpsi_sig):
        h_cos2pisys.Fill(t.vtx_cos2pisys)

    if (             cut_trkp_costhe and cut_trkm_costhe and cut_trkp_p and cut_trkm_p and
        cut_cospipi and cut_cos2pisys and cut_pi_PID and cut_mjpsi_sig):
        h_ngam.Fill(t.ngam)

    
def write_histograms():
    h_evtflw.Write()
    h_mrecpipi.Write()
    h_mpipi.Write()
    h_pip_p.Write()
    h_pim_p.Write()
    h_pip_costhe.Write()
    h_pim_costhe.Write()
    h_cospipi.Write()
    h_cos2pisys.Write()
    h_ngam.Write()

    
def select_jpsi_to_invisible(t):
    h_evtflw.Fill(0) 

    if not (t.ngam == 0):
        return False
    h_evtflw.Fill(1) 
    
    if not ( abs(math.cos(t.trkp_theta)) < 0.8 and abs(math.cos(t.trkm_theta)) < 0.8):
        return False
    h_evtflw.Fill(2) 

    if not (abs(t.trkp_p) < 0.45 and abs(t.trkm_p) < 0.45):
        return False 
    h_evtflw.Fill(3) 

    if not (t.prob_pip > t.prob_kp and t.prob_pip > 0.001 and
            t.prob_pim > t.prob_km and t.prob_pim > 0.001):
        return False
    h_evtflw.Fill(4)

    if not (t.vtx_cospipi < 0.95):
        return False
    h_evtflw.Fill(5)

    if not (t.vtx_cos2pisys < 0.9):
        return False
    h_evtflw.Fill(6)

    if not (t.vtx_mrecpipi > 3.0 and t.vtx_mrecpipi < 3.2):
        return False
    h_evtflw.Fill(7)
    
    return True
    
    
if __name__ == '__main__':
    main()
