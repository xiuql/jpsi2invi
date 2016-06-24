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

    ./run_plotter.py infile outfile 

AUTHOR 
    SHI Xin <shixin@ihep.ac.cn> 

DATE
    June 2016 
\n''')

    
def main():

    args = sys.argv[1:]
    if len(args) < 2:
        return usage()
    
    infile = args[0]
    outfile = args[1]

    fin = ROOT.TFile(infile)
    t = fin.Get('tree')
    entries = t.GetEntriesFast()

    h_mrecpipi = ROOT.TH1D('h_mrecpipi', 'mrecpipi', 100, 3.04, 3.16) 
    for jentry in range(entries):
        # get the next tree in the chain and verify
        ientry = t.LoadTree(jentry)
        if ientry < 0:
            break
        # copy next entry into memory and verify
        nb = t.GetEntry(jentry)
        if nb<=0:
            continue

        h_mrecpipi.Fill(t.vtx_mrecpipi)

        if ientry > 100:
            break

    fout = ROOT.TFile(outfile, "RECREATE")
    h_mrecpipi.Write()
    fout.Close()


if __name__ == '__main__':
    main()
