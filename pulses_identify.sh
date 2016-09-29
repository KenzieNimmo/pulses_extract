#!/bin/bash -x

if [ $# -ne 1 ]; then
   echo "Pipeline to process Arecibo data of FRB121102."
   echo "The pipeline identify interesting pulses, store them in a SinglePulse.hdf5 database and produces diagnostic plots."
   echo ""
   echo "Usage: sh pulses_identify.sh fits_filename"
   echo "fits_filename contains the full path of the file."
   echo ""
   echo "The script will store the output in the path defined by FITS_PATH,"
   echo "in a directory named as the fits file base."
   exit
fi

echo "Pipeline pulses_identify.sh starting..."
date

FITS=$1
FITS_PATH=/psr_temp/hessels/AO-FRB/P3054/FRB_pipeline/output/
FITS_BASENAME=${FITS##*/}; FITS_BASENAME=${FITS_BASENAME%_subs_0001.fits}
SCRIPT_DIR="$( cd -P "$( dirname "$0" )" && pwd )"

cd $FITS_PATH 
mkdir $FITS_BASENAME
cd $FITS_BASENAME
mkdir TEMP
cd TEMP

prepsubband -nsub 64 -noscales -nooffsets -nobary -lodm 461.0 -numdms 200 -numout 88473600 -dmstep 1.0 -o ${FITS_BASENAME}_TOPO $FITS  #-zero_dm

sh ${SCRIPT_DIR}/periodicity.sh $FITS_BASENAME $FITS

single_pulse_search.py -t 6.0 -b -m 150 *.dat
cp ${FITS_BASENAME}_TOPO_DM561.00.dat ..
cp ${FITS_BASENAME}_TOPO_DM561.00.inf ..

python ${SCRIPT_DIR}/pulse_extract.py -fits $FITS -store_events -idL ${FITS_BASENAME} -store_dir $FITS_PATH/$FITS_BASENAME \
  -folder $FITS_PATH/$FITS_BASENAME/TEMP -plot_pulses


#cd $FITS_PATH/$FITS_BASENAME
#rm -r TEMP

date
echo "Pipeline pulses_identify.sh finished"


