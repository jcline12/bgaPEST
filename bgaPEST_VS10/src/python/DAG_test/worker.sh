#!/bin/sh
cInd=$1

tar xzf data.tar
cd data

python con_test.py

tar czf results_$cInd.tar *.out
