#!/usr/bin/env bash

# Main driver to build programs
# Author Xin Shi <shixin@ihep.ac.cn>
# Created [2016-03-28 Mon 08:19]


if [[ $# -eq 0 ]]; then 
    printf "NAME\n\tbuild.sh - Main driver to build programs\n"
    printf "\nSYNOPSIS\n"
    printf "\n\t%-5s\n" "./build.sh [OPTION]" 
    printf "\nOPTIONS\n" 
    printf "\n\t%-5s  %-40s\n"  "0.1"  "build Jpsi2invi " 
fi

option=$1

case $option in 
    0.1) echo "Building Jpsi2invi module..."
	 cd Analysis/Physics/PsiPrime/Jpsi2invi/Jpsi2invi-00-00-01/cmt 
	 gmake  
       ;;
esac

