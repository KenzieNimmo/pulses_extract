#!/bin/bash -x

if [ $# -ne 1 ]; then
   echo "Pipeline to process Arecibo data of FRB121102."
   echo "The pipeline use pulses stored in a SinglePulse.hdf5 database to chop raw data around the pulses."
   echo ""
   echo "Usage: sh chop_raw_data.sh database_filename"
   echo "database_filename contains the full path of the SinglePulse.hdf5 file."
   echo ""
   echo "The script will store the output in the path defined by FITS_PATH,"
   echo "in a directory named fits."
   exit
fi

DB=$1
FITS_PATH=/psr_temp/hessels/AO-FRB/P3054/FRB_pipeline/output/
FITS_BASENAME=${DB##*/}; FITS_BASENAME=${FITS_BASENAME%SinglePulse.hdf5}
SCRIPT_DIR="$( cd -P "$( dirname "$0" )" && pwd )"

cd $FITS_PATH/$FITS_BASENAME
mkdit fits

python ${SCRIPT_DIR}/pulse_extract.py -pulses_database $DB -pulses_checked $FITS_PATH/${FITS_BASENAME}/pulses_list.txt \
  -store_dir $FITS_PATH/$FITS_BASENAME/fits -extract_raw -raw_basename $FITS_BASENAME

