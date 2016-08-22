#!/bin/bash
# Script to submit PBS jobs 
# Author: SHI Xin <shixin@ihep.ac.cn>
# Date: [2016-08-15 Mon 09:25] 

#PBS -N jmcpsip12
#PBS -q "besq@torqsrv"
#PBS -o $HOME/bes/jpsi2invi/v0.1/run/log/mc_psip12/jpsi2invi_mc_psip12.log 
##PBS -t 1-394%200
#PBS -t 41,124,162,148,115,158,178,181,182,79,136,48,151,22,110,127,220,141,206,101,138,146,200,209,39,34,37,123,175,108,159,107,97,102,135,207,105,168,166,36,184,185,13,201,25,125,55,214,21,20,40,47,92,106,31 


date

hostname

cd $HOME/bes/jpsi2invi/v0.1/
source setup.sh 
cd pbs 
boss.exe jobOptions_jpsi2invi_mc_psip12.pbs

date


