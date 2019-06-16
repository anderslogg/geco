#!/usr/bin/env bash

# This script is for computing postprocess data and ergoregions on a sequence of GECO solutions.
# Usage: enter
#   >> ./run-geco-postprocess
# from the sequence directory (e.g. tss_L08)

# From docker terminal:
# first change geco-postprocess-data to
# /home/fenics/shared/geco/bin/geco-postprocess-data below. Then, run
# >> /home/fenics/shared/geco/bin/run-geco-postprocess.sh

######################################################################
######################################################################

# Handle Optional arguments

force_flag=0

print_usage() {
  printf "Usage: ..."
}

while getopts 'f' flag; do
  case "${flag}" in
      f) force_flag=1 ;;
      --) shift
	  break ;;      
      *) print_usage
         exit 1 ;;
  esac
  shift
done

if [ ${force_flag} -eq 1 ]; then
    echo "Force flag set. Overwriting postprocessing files..."
fi

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
    
    echo "...entering $PWD"
    
    # If operation is data, check for existing ppdata file
    if [[ ("$1" = "data") && ( (! -f "ppdata.csv") || (${force_flag} -eq 1) ) ]]; then
	echo "   computing postprocess data..."
	geco-postprocess-data
    fi
    # Otherwise, execute operation.
    if [ "$1" != "data" ]; then
	echo "   executing postprocess $1..."
	geco-postprocess-$1
    fi
done

# Print a nice message
echo
echo "Postprocessing complete!"
