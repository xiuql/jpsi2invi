#!/usr/bin/env python
"""
Run Plotter 
"""

__author__ = "SHI Xin <shixin@ihep.ac.cn>"
__copyright__ = "Copyright (c) SHI Xin"
__created__ = "[2016-06-20 Mon 08:08]" 

import sys
import os
import ROOT 

def usage():
    sys.stdout.write('''
NAME
    run_plotter.py 

SYNOPSIS

    ./run_plotter.py  

AUTHOR 
    SHI Xin <shixin@ihep.ac.cn> 

DATE
    June 2016 
\n''')

    
def main():

    infile = 'run/data/jpsi2invi_data-1.root'
    outfile = 'plotter.root'



    
    myfile = ROOT.TFile(infile)
    fout = ROOT.TFile(outfile, "RECREATE")
    # retrieve the ntuple of interest
    t = myfile.Get('tree')
    entries = t.GetEntriesFast()
    print entries 

    h_mrecpipi = ROOT.TH1D('h_mrecpipi', 'mrecpipi', 100, 3.04, 3.16) 
    for jentry in range(entries):
        # get the next tree in the chain and verify
        ientry = t.LoadTree(jentry)

        print 'jentry: %s, ientry: %s' %(jentry, ientry)

        if ientry < 0:
            break

        # copy next entry into memory and verify
        nb = t.GetEntry(jentry)
        if nb<=0:
            continue

        # use the values directly from the tree
        print t.event

        h_mrecpipi.Fill(t.vtx_mrecpipi)

        if ientry > 10:
            break

    h_mrecpipi.Write()
    fout.Close()

    
if __name__ == '__main__':
    main()
