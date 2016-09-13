#!/usr/bin/env bash


cd besenv
source setupCVS.sh
source setupCMT.sh
cmt br cmt config
source setup.sh

cd ../TestRelease/TestRelease-00-00-81/cmt/
cmt br cmt config
source setup.sh
cd ../../../


