#!/bin/bash -x

if [ $# -ne 2 ]; then
   echo "Pipeline to process Arecibo data of FRB121102."
   echo "The pipeline use pulses stored in a SinglePulse.hdf5 database to chop raw data around the pulses."
   echo ""
   echo "Usage: sh chop_raw_data.sh database_filename raw_fits_files_basename"
   echo "database_filename contains the full path of the SinglePulse.hdf5 file."
   echo "raw_fits_files_basename contains the basename of the raw fits files."
   echo ""
   echo "The script will store the output in the path defined by FITS_PATH,"
   echo "in a directory named fits."
   exit
fi

DB=$1
FITS=$2
FITS_BASENAME=${FITS%_0???.fits}
DB_PATH=${DB%SinglePulses.hdf5}
SCRIPT_DIR="$( cd -P "$( dirname "$0" )" && pwd )"

cd $DB_PATH
mkdir fits

python ${SCRIPT_DIR}/pulse_extract.py -pulses_database $DB -pulses_checked $DB_PATH/pulses_list.txt \
  -store_dir $DB_PATH -extract_raw $FITS_BASENAME

