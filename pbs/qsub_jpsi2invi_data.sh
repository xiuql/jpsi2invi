#!/bin/bash
# Script to submit PBS jobs 
# Auther: SHI Xin <shixin@ihep.ac.cn>
# Date: [2016-06-02 Thu 13:29]  

#PBS -N jdata
##PBS -j oe
#PBS -q "besq@torqsrv"
#PBS -o $HOME/bes/jpsi2invi/v0.1/run/log/jpsi2invi_data-test.log 
##PBS -t 1-10 


date

hostname

cd $HOME/bes/jpsi2invi/v0.1/
source setup.sh 
cd pbs 
boss.exe jobOptions_jpsi2invi_data.pbs

date


