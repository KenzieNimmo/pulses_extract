#!/bin/bash -x

if [ $# -ne 2 ]; then
   echo "Pipeline to process Arecibo data of FRB121102."
   echo "The pipeline uses pulses stored in a HDF5 database to chop raw data around the pulses."
   echo ""
   echo "Usage: bash chop_raw_data.sh OBS_ID"
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

SCRIPT_DIR="$( cd -P "$( dirname "$0" )" && pwd )"
OBS_ID="$1"
OUT_DIR="$GENERAL_OUT_DIR/$OBS_ID"
DB_FILE="$OBS_ID.hdf5"
CAL_FILE="${OBS_ID:0: -1}$((${OBS_ID: -1} - 1))_cal_0001.fits"
FITS_NAME="${OBS_ID}_subs_0001.fits"

#Check that database exists
if [ ! -e $OUT_DIR/pulses/$DB_FILE ]; then
  echo ""
  echo "ATTENTION! HDF5 database $DB_FILE not found. Exiting..."
  exit 1
fi
#Check that the calibration file exists
if [ ! -e $RAW_DIR/$CAL_FILE ]; then
  echo ""
  echo "ATTENTION! Calibration file $CAL_FILE not found. Exiting..."
  exit 1
fi
#Check that subbanded fits file exists
if [ ! -e $SUB_DIR/$FITS_NAME ]; then
  echo ""
  echo "ATTENTION! Subbanded fits file $FITS_NAME not found. Exiting..."
  exit 1
fi

cp $RAW_DIR/$CAL_FILE $OUT_DIR/obs_data
cp $SUB_DIR/$FITS_NAME $OUT_DIR/obs_data

#Create raw fits files
python ${SCRIPT_DIR}/pulses_extract.py -db_name $DB_FILE -pulses_database -pulses_checked ${OUT_DIR}/pulses/${OBS_ID}_pulses.txt \
  -store_dir $OUT_DIR/pulses -extract_raw $RAW_DIR/$OBS_ID

#RFI masks
cd $OBS_ID/pulses
for puls in `ls`; do 
  if [ -d $puls ]; then
    cd $puls
    for fits in `ls *.fits`; do 
      rfifind -blocks 10 -noweights -noscales -nooffsets -o ${fits%.fits} $fits
    done
    cd $OBS_ID/pulses
  fi    
done

date
echo "Pipeline chop_raw.sh finished"
