#!/bin/bash -x

if [ $# -ne 1 ]; then
   echo "Pipeline to process Arecibo data of FRB121102."
   echo "The pipeline identify interesting pulses, store them in a SinglePulse.hdf5 database and produces diagnostic plots."
   echo ""
   echo "Usage: bash pulses_identify.sh fits_filename"
   echo "NB: use the command bash to run the pipeline"
   exit
fi

#Check that bash is used
if [ ! "$BASH_VERSION" ] ; then
    echo "Execute the script using the bash command. Exiting..."
    exit 1
fi

echo "Pipeline pulses_identify.sh starting..."
date

#Setting variables
SUB_DIR="/psr_archive/hessels/hessels/AO-FRB/subbanded_data"
GENERAL_OUT_DIR="/psr_archive/hessels/hessels/AO-FRB/pipeline_products"
RAW_DIR="/psr_archive/hessels/hessels/AO-FRB/raw_data"

ORIGINAL_FITS_FILE="$1"
FITS_NAME=${ORIGINAL_FITS_FILE##*/}
FITS_ID=${FITS_NAME%_subs_0001.fits}
SCRIPT_DIR="$( cd -P "$( dirname "$0" )" && pwd )"
CAL_FILE="${FITS_ID:0: -1}$((${FITS_ID: -1} - 1))_cal_0001.fits"
OUT_DIR="$GENERAL_OUT_DIR/$FITS_ID"
FITS_FILE=$OUT_DIR/obs_data/$FITS_NAME

#Check that subbanded fits file and calibration file exist
if [ ! -e $SUB_DIR/$FITS_NAME ]; then
  echo ""
  echo "ATTENTION! Subbanded fits file $FITS_NAME not found. Exiting..."
  exit 1
fi
if [ ! -e $RAW_DIR/$CAL_FILE ]; then
  echo ""
  echo "ATTENTION! Calibration file $CAL_FILE not found. Exiting..."
  exit 1
fi

#Set up the output folder
mkdir $OUT_DIR
cd $OUT_DIR

mkdir obs_data
mkdir stat_plots
mkdir pulses
mkdir periodic_cands
mkdir TEMP
cd TEMP

#Run the search
cp $RAW_DIR/$CAL_FILE $OUT_DIR/obs_data
cp $SUB_DIR/$FITS_NAME $OUT_DIR/obs_data

prepsubband -nsub 64 -noscales -nooffsets -nobary -lodm 461.0 -numdms 200 -numout 88473600 -dmstep 1.0 -o ${FITS_ID}_TOPO -zerodm $FITS_FILE

bash ${SCRIPT_DIR}/periodicity.sh $OUT_DIR $FITS_FILE

single_pulse_search.py -t 6.0 -b -m 150 *.dat
cp ${FITS_ID}_TOPO_DM561.00.dat $OUT_DIR/obs_data
cp ${FITS_ID}_TOPO_DM561.00.inf $OUT_DIR/obs_data

python ${SCRIPT_DIR}/pulse_extract.py -db_name ${FITS_ID}.hdf5 -fits $FITS_FILE -store_events -idL ${FITS_ID}_TOPO -store_dir $OUT_DIR/pulses \
  -folder $OUT_DIR/TEMP -plot_pulses -create_pulse_dirs

#Remove meta products
cd $OUT_DIR
rm -rf TEMP

date
echo "Pipeline pulses_identify.sh finished"


