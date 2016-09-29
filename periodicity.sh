#!/bin/bash -x

if [ $# -ne 2 ]; then
   echo "Script for periodicity search on Arecibo FRB121102 data."
   echo ""
   echo "Usage: sh periodicity.sh output_folder fits_filename"
   echo "fits_filename contains the full path of the file."
   echo "output_folder is the directory to store the plots."
   echo "The script assumes that meta files are in a TEMP subdirectory."
   exit
fi

FOLDER=$1
FITS=$2
SCRIPT_DIR="$( cd -P "$( dirname "$0" )" && pwd )"

cd $FOLDER/TEMP

ls *.dat | xargs -n 1 realfft
ls *.fft | xargs -n 1 accelsearch -zmax 0 
python $SCRIPT_DIR/periodic_candidates_plot.py -folder $FOLDER -fits $FITS
