unzip data.zip
unzip Python27.zip > nul
cd data

..\Python27\python condor_single_run.py %1 %2

..\Python27\python zip_results.py %1
