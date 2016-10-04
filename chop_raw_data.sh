#!/bin/bash -x

if [ $# -ne 2 ]; then
   echo "Pipeline to process Arecibo data of FRB121102."
   echo "The pipeline uses pulses stored in a HDF5 database to chop raw data around the pulses."
   echo ""
   echo "Usage: bash chop_raw_data.sh database_filename"
   echo "NB: use the command bash to run the pipeline"
   exit
fi

#Check that bash is used
if [ ! "$BASH_VERSION" ] ; then
    echo "Execute the script using the bash command. Exiting..."
    exit 1
fi

echo "Pipeline chop_raw_data.sh starting..."
date

#Setting variables
SUB_DIR="/psr_archive/hessels/hessels/AO-FRB/subbanded_data"
GENERAL_OUT_DIR="/psr_archive/hessels/hessels/AO-FRB/pipeline_products"
RAW_DIR="/psr_archive/hessels/hessels/AO-FRB/raw_data"





#Check that SinglePulse.hdf5 exists
if [ ! -e $SUB_DIR/$FITS_NAME ]; then
  echo ""
  echo "ATTENTION! Subbanded fits file $FITS_NAME not found. Exiting..."
  exit 1
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

#RFI masks
cd fits
for fits in `ls *.fits`; do 
  rfifind -blocks 10 -noweights -noscales -nooffsets -o ${fits%.fits} $fits
done

date
echo "Pipeline chop_raw.sh finished"
