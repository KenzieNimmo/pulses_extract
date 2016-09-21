#!/bin/bash -x

if [ $# -ne 1 ]; then
   echo "Pipeline to process Arecibo data of FRB121102."
   echo ""
   echo "Usage: sh FRB_pipeline.sh fits_filename"
   echo "fits_filename contains the full path of the file."
   echo ""
   echo "The script will store the output in the same fits file folder,"
   echo "in a directory named as the fits file."
   exit
fi

FITS=$1
FITS_PATH=${FITS%/*}
FITS_BASENAME=${FITS##*/}; FITS_BASENAME=${FITS_BASENAME%.*}
SCRIPT_DIR="$( cd -P "$( dirname "$0" )" && pwd )"

cd $FITS_PATH 
mkdir $FITS_BASENAME
cd $FITS_BASENAME
mkdir TEMP
cd TEMP

prepsubband -nsub 64 -noscales -nooffsets -nobary -lodm 461.0 -numdms 200 -dmstep 1.0 -o ${FITS_BASENAME}_TOPO $FITS

single_pulse_search.py -t 6.0 -b -m 150 *.dat
mv ${FITS_BASENAME}_TOPO_DM561.00.dat ..

python ${SCRIPT_DIR}/pulse_extract.py $FITS -store_events -idL ${FITS_BASENAME} -store_dir $FITS_PATH/$FITS_BASENAME -folder $FITS_PATH/$FITS_BASENAME/TEMP -plot_pulses


#cd $FITS_PATH/$FITS_BASENAME
#rm -r TEMP



