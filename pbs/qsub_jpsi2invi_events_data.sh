#!/bin/bash
# Script to submit PBS jobs 
# Author: SHI Xin <shixin@ihep.ac.cn>
# Date: [2016-07-14 Thu 09:31] 

#PBS -N sel
#PBS -q "besq@torqsrv"
#PBS -o $HOME/bes/jpsi2invi/v0.1/run/log/jpsi2invi_events_data.log 
#PBS -t 1-633%200

date

hostname

cd $HOME/bes/jpsi2invi/v0.1/
source setup.sh

./python/sel_events.py run/data/jpsi2invi_data-$PBS_ARRAYID.root run/events/jpsi2invi_data-$PBS_ARRAYID.root

date




