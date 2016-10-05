#!/bin/bash

# Maximum number of jobs
if [ -z $2 ]; then
    MAXJOBS=4
else
    MAXJOBS=$2
fi

# Get job ID from scrint
function jobidfromstring()
{
    local STRING;
    local RET;

    STRING=$1;
    RET="$(echo $STRING | sed 's/^[^0-9]*//' | sed 's/[^0-9].*$//')"

    echo $RET;
}

# Ready to spawn job
function clearToSpawn
{
    local JOBCOUNT="$(jobs -r | grep -c .)"
    if [ $JOBCOUNT -lt $MAXJOBS ] ; then
        echo 1;
        return 1;
    fi

    echo 0;
    return 0;
}

# Loop over commands
JOBLIST=""
while read line; do
  while [ `clearToSpawn` -ne 1 ] ; do
      sleep 0.5
  done

  # Run command
  #echo "Running" $line
  echo $line | sh &

  LASTJOB=`jobidfromstring $(jobs %%)`
  JOBLIST="$JOBLIST $LASTJOB"
done < $1
# Wait for jobs to finish
for JOB in $JOBLIST ; do
    wait %$JOB
done
