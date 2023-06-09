#!/usr/bin/env bash

# This script is for computing postprocess data and ergoregions on a sequence of GECO solutions.
# Usage: enter
#   >> ./run-geco-postprocess
# from the sequence directory (e.g. tss_L08)

# From docker terminal:
# first change geco-postprocess-data to
# /home/fenics/shared/geco/bin/geco-postprocess-data below. Then, run
# >> /home/fenics/shared/geco/bin/run-geco-postprocess.sh


# For each step_ directory, compute geco-postprocess-data in that directory.

# Directories
CURRENTDIR=$( pwd )
ALLSTEPS=$( ls adaptive_solver/ | grep '^step_') #grep '9[0-9]$')  #'[0-9][0,5]$') #

# Optional: evaluate only subset of directories. 
START=${2:-0}
STOP=${3:-1}

COUNTER=0
for s in `seq $START $STOP`;
do
    STEPS[$COUNTER]='step_'$(printf "%03d" $s)
    let COUNTER=COUNTER+1
done

if [ "$#" -gt 1 ]; then
    ALLSTEPS=${STEPS[@]}
fi

# Run postprocessing code in desired steps
for STEP in $ALLSTEPS
do
    cd $CURRENTDIR/adaptive_solver/$STEP
    # If operation is data, check for existing ppdata file
    if [ "$1" = "data" ] && [ ! -f "ppdata.csv" ]; then
	geco-postprocess-data
    fi
    # Otherwise, execute operation.
    if [ "$1" != "data" ]; then
	geco-postprocess-$1
    fi
done

# Print a nice message
echo
echo "Postprocessing complete!"
