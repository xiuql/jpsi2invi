#!/usr/bin/env bash

# Main driver to submit jobs 
# Author Xin Shi <shixin@ihep.ac.cn>
# Created [2016-03-28 Mon 07:58]


if [[ $# -eq 0 ]]; then 
    printf "NAME\n\tsubmit.sh - Main driver to submit jobs\n"
    printf "\nSYNOPSIS\n"
    printf "\n\t%-5s\n" "./submit.sh [OPTION]" 
    printf "\nOPTIONS\n" 
    printf "\n\t%-5s  %-40s\n"  "0.1"  "run data sample" 
fi


option=$1

case $option in 
    0.1) echo "Running on data sample..."
	 # boss.exe jobOptions_jpsi2invi.txt
	 # ./python/get_samples.py  /bes3fs/offline/data/664p03/psip/dst $HOME/bes/jpsi2invi/v0.1/run/samples/data_664p03_psip.txt 40G 
	 qsub pbs/qsub_jpsi2invi_data.sh  
       ;;
esac



