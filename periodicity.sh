#!/bin/bash -x




FOLDER=$1


cd $FOLDER
ls *.dat | xargs -n 1 realfft

ls *.fft | xargs -n 1 accelsearch -zmax 0 

python ACCEL_sift.py > $FOLDER/periodicity_cands.txt








prepfold -nsub 64 -p $PERIOD -dm $DM -o 