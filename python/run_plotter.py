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
from progressbar import Bar, Percentage, ProgressBar
from time import time 

TEST=True 
#TEST=False

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

    if TEST:
        sys.stdout.write('Run in TEST mode! \n')
    
    args = sys.argv[1:]
    if len(args) < 2:
        return usage()
    
    infile = args[0]
    outfile = args[1]

    fin = ROOT.TFile(infile)
    t = fin.Get('tree')
    entries = t.GetEntriesFast()

    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=entries).start()
    time_start = time()
    
    h_mrecpipi = ROOT.TH1D('h_mrecpipi', 'mrecpipi', 100, 3.04, 3.16) 
    for jentry in xrange(entries):
        pbar.update(jentry+1)
        # get the next tree in the chain and verify
        ientry = t.LoadTree(jentry)
        if ientry < 0:
            break
        # copy next entry into memory and verify
        nb = t.GetEntry(jentry)
        if nb<=0:
            continue

        h_mrecpipi.Fill(t.vtx_mrecpipi)

        if TEST and ientry > 100000:
            break

    fout = ROOT.TFile(outfile, "RECREATE")
    h_mrecpipi.Write()
    fout.Close()
    pbar.finish()
    
    dur = duration(time()-time_start)
    sys.stdout.write(' \nDone in %s. \n' % dur) 


def duration(seconds):
    seconds = long(round(seconds))
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    years, days = divmod(days, 365.242199)
 
    minutes = long(minutes)
    hours = long(hours)
    days = long(days)
    years = long(years)
 
    duration = []
    if years > 0:
        duration.append('%d year' % years + 's'*(years != 1))
    else:
        if days > 0:
            duration.append('%d day' % days + 's'*(days != 1))
        if hours > 0:
            duration.append('%d hour' % hours + 's'*(hours != 1))
        if minutes > 0:
            duration.append('%d minute' % minutes + 's'*(minutes != 1))
        if seconds > 0:
            duration.append('%d second' % seconds + 's'*(seconds != 1))
    return ' '.join(duration)

    
if __name__ == '__main__':
    main()
