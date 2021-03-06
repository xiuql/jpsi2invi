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
    printf "\n\t%-5s  %-40s\n"  "0.1.2"    "Split data sample with each group 20G"
    printf "\n\t%-5s  %-40s\n"  "0.1.3"    "Submit PBS jobs on data"
    printf "\n\t%-5s  %-40s\n"  "0.1.4"    "Check PBS jobs on data."
    printf "\n\t%-5s  %-40s\n"  "0.1.5"    "Select events."
    printf "\n\t%-5s  %-40s\n"  "0.1.6"    "Submit events jobs on data."
    printf "\n\t%-5s  %-40s\n"  "0.1.7"    "Check events jobs on data."
    printf "\n\t%-5s  %-40s\n"  "0.1.8"    "Merge events files." 
    printf "\n\t%-5s  %-40s\n"  "0.1.9"    "Plot summary with data."
    printf "\n\t%-5s  %-40s\n"  "0.2"      "[run on MC sample]"
    printf "\n\t%-5s  %-40s\n"  "0.2.1"    "Run with a few samples"
    printf "\n\t%-5s  %-40s\n"  "0.2.2"    "Split psi(2S) MC sample with each group 20G"
    printf "\n\t%-5s  %-40s\n"  "0.2.3"    "Submit PBS jobs on psi(2S) MC sample"     
    printf "\n\t%-5s  %-40s\n"  "0.2.4"    "Check PBS jobs on psi(2S) MC sample"     
    printf "\nAUTHOR\n"
    printf "\n\t%-5s\n" "SHI Xin <shixin@ihep.ac.cn>"
    printf "\nDATE\n"
    printf "\n\t%-5s\n" "AUGUST 2016"     
}


if [[ $# -eq 0 ]]; then
    usage
fi


option=$1

case $option in 
    0.1) echo "Running on data sample..."
	 ;;

    0.1.1) echo "Run with a few events ..."
	   boss.exe jobOptions_jpsi2invi.txt
	   ;;
    
    0.1.2) echo "Split data sample with each group 20G ..."
	   ./python/get_samples.py  /bes3fs/offline/data/664p03/psip/dst $HOME/bes/jpsi2invi/v0.1/run/samples/data_664p03_psip.txt 20G
	   # made 633 groups 
	   ;;

    0.1.3) echo "Submit PBS jobs on data..."
	   mkdir run/data
	   mkdir run/log/data  
	   qsub pbs/qsub_jpsi2invi_data.sh  
	   ;;

    0.1.4) echo "Check PBS jobs on data..."
	   ./python/chk_pbsjobs.py $HOME/bes/jpsi2invi/v0.1/run/data  633
	   ;;
    
    0.1.5) echo  "Select events on data..."
	   ./python/sel_events.py  run/data/jpsi2invi_data-1.root  run/events/jpsi2invi_data-1.root 
	   ;; 

    0.1.6) echo "Submit selection PBS jobs on data..."
	   mkdir run/events
	   mkdir run/log/events  
	   qsub pbs/qsub_jpsi2invi_events_data.sh  
	   ;;

    0.1.7) echo "Check PBS jobs on events data..."
	   ./python/chk_pbsjobs.py run/events  633
	   ;;

    0.1.8) echo  "Merge root files..."
	   ./python/mrg_rootfiles.py  run/events run/hist 
	   ;; 

    0.1.9) echo  "Plot summary with data..."
	   ./python/plt_summary.py 
	   ;; 


    
    0.2) echo "Running on MC sample..."
	 ;;

    0.2.1) echo "Run with a few events ..."
	   boss.exe jobOptions_jpsi2invi.txt
	   ;;
    0.2.2) echo "Split psi(2S) MC sample with each group 20G ..."
	   ./python/get_samples.py  /bes3fs/offline/data/664p03/psip/12mc/dst $HOME/bes/jpsi2invi/v0.1/run/samples/mc_664p03_psip_12mc.txt 20G
	   # made 394 groups 
	   ;;

    0.2.3) echo "Submit PBS jobs on psi(2S) MC sample..."
	   mkdir run/mc_psip12
	   mkdir run/log/mc_psip12  
	   qsub pbs/qsub_jpsi2invi_mc_psip12.sh  
	   ;;

    0.2.4) echo "Check PBS jobs on psi(2S) MC sample..."
	   ./python/chk_pbsjobs.py $HOME/bes/jpsi2invi/v0.1/run/mc_psip12  394 
	   ;;
    
esac

