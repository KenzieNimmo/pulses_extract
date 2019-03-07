#!/bin/bash -x

if [ $# -ne 1 ] && [ $# -ne 2 ]; then
   echo "Pipeline to process Effelsberg data of FRB180814 (R2) on DRAGNET."
   echo "The pipeline identify interesting pulses, store them in a SinglePulse.hdf5 database and produces diagnostic plots."
   echo ""
   echo "Usage: bash pulses_identify.sh fil_filename"
   echo "NB: use the command bash to run the pipeline"
   echo "Use the argument -single_core to avoid the pipeline to run in parallel"
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
 SUB_DIR="/data2/nimmo/FRB_R2_2019/R2_data"
 GENERAL_OUT_DIR="/data2/nimmo/FRB_R2_2019/pipeline_products"

 ORIGINAL_FIL_FILE="$1"
 FIL_NAME=${ORIGINAL_FIL_FILE##*/}
 FIL_ID=${FIL_NAME%.out.fil}
 SCRIPT_DIR="~/pulses_extract/src"
 OUT_DIR="$GENERAL_OUT_DIR/$FIL_ID"
 FIL_FILE="$SUB_DIR/$FIL_NAME"

 #Check that the filterbank file exists
 if [ ! -e $FIL_FILE ]; then
   echo ""
   echo "ATTENTION! Filterbank file $FIL_NAME not found. Exiting..."
   exit 1
 fi

 #Set up the output folder
 mkdir $OUT_DIR
 cd $OUT_DIR
 mkdir obs_data
 mkdir pulses
 mkdir periodic_cands
 mkdir TEMP
 cd TEMP

 #Create .dat files
 echo ".dat files creating..."
 SECONDS=0
 prepsubband -nsub 640 -nobary -lodm 120 -numdms 400 -dmstep 0.3 -o ${FIL_ID}_TOPO -zerodm $FIL_FILE
 duration=$SECONDS
 echo ".dat files created. Time taken: $(($duration / 60)) m"

 #Parallelisation parameters
 if [ "$2" == "-single_core" ]; then
   n_cores=1
 else
   n_cores=`lscpu -p | egrep -v '#' | sort -u -t, -k 2,4 | wc -l`
 fi

 #Create .singlepulse files
 echo ".singlepulse files creating..."
 SECONDS=0
 ls *.dat | awk '{printf("single_pulse_search.py -t 6.0 -b %s\n",$1)}' > jobs.txt
 bash $SCRIPT_DIR/parallel.sh jobs.txt $n_cores >/dev/null
 duration=$SECONDS
 echo ".singlepulse files created. Time taken: $(($duration / 60)) m"

 #Create database and plots
 echo "Database and plots creating..."
 SECONDS=0
 #Single pulse candidates
 python ${SCRIPT_DIR}/pulses_extract.py -db_name ${FIL_ID}.hdf5 -fits $FIL_FILE -store_events -idL ${FIL_ID}_TOPO -store_dir $OUT_DIR/pulses \
   -folder $OUT_DIR/TEMP -plot_pulses -plot_statistics -parameters_id FRB180814_PSRIX > /dev/null
 duration=$SECONDS
 echo "Database and plots created. Time taken: $(($duration / 60)) m"

 #Remove meta products
 cd $OUT_DIR
 rm -rf TEMP

 echo "Pipeline pulses_identify.sh finished"
 date
