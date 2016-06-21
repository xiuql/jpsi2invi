#!/usr/bin/env bash

# Main driver to submit jobs 
# Author SHI Xin <shixin@ihep.ac.cn>
# Created [2016-03-28 Mon 07:58]


usage() {
    printf "NAME\n\tsubmit.sh - Main driver to submit jobs\n"
    printf "\nSYNOPSIS\n"
    printf "\n\t%-5s\n" "./submit.sh [OPTION]" 
    printf "\nOPTIONS\n" 
    printf "\n\t%-5s  %-40s\n"  "0.1"      "[run data sample]" 
    printf "\n\t%-5s  %-40s\n"  "0.1.1"    "Run with a few samples" 
    printf "\n\t%-5s  %-40s\n"  "0.1.2"    "Split data sample with each group 40G"
    printf "\n\t%-5s  %-40s\n"  "0.1.3"    "Submit PBS jobs on data"
    printf "\n\t%-5s  %-40s\n"  "0.1.4"    "Check PBS jobs on data."
    printf "\n\t%-5s  %-40s\n"  "0.1.4.1"  "arXiv prod files."
    printf "\n\t%-5s  %-40s\n"  "0.1.5"    "Merge root files."
    printf "\n\t%-5s  %-40s\n"  "0.1.6"    "run Plotter on data."
    printf "\nAUTHOR\n"
    printf "\n\t%-5s\n" "SHI Xin <shixin@ihep.ac.cn>" 
}


if [[ $# -eq 0 ]]; then
    usage
fi


option=$1

case $option in 
    0.1) echo "Running on data sample..."
	 #qsub pbs/qsub_jpsi2invi_data.sh  
	 ;;

    0.1.1) echo "Run with a few samples ..."
	   boss.exe jobOptions_jpsi2invi.txt
	   ;;
    
    0.1.2) echo "Split data sample with each group 20G ..."
	   ./python/get_samples.py  /bes3fs/offline/data/664p03/psip/dst $HOME/bes/jpsi2invi/v0.1/run/samples/data_664p03_psip.txt 20G 
	   ;;

    0.1.3) echo "Submit PBS jobs on data..."
	   qsub pbs/qsub_jpsi2invi_data.sh  
	   ;;

    0.1.4) echo "Check PBS jobs on data..."
	   ./python/chk_pbsjobs.py $HOME/bes/jpsi2invi/v0.1/run/data  633
	   ;;
    
    0.1.4.1) echo  "Arxiv production files to v1 ..."
	   mkdir run/data/v1
	   mv run/data/*.root run/data/v1

	   mkdir run/log/v1
	   mv run/log/*.log* run/log/v1

	   mkdir run/samples/v1
	   mv run/samples/*.txt run/samples/v1 
	   ;; 

    0.1.5) echo  "Merge root files..."
	   ./python/mrg_rootfiles.py  $HOME/bes/jpsi2invi/v0.1/run/data 
	   ;; 

    0.1.6) echo  "run Plotter on data..."
	   #./bin/runPlotter -i run/data/jpsi2invi_data-1.root -o plotter_data.root 
	   ./python/run_plotter.py 
	   ;;
    
esac

