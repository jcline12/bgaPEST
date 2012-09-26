#!/bin/sh
cInd=$1

export WINEPREFIX=$_CONDOR_SCRATCH_DIR

wine unzip data.zip
wine unzip Python27.zip > nul
cd data

wine ../Python27/python condor_single_run.py $cInd

python zip_results.py $cInd
