#!/usr/bin/env bash

# This script is for computing postprocess data on a sequence of GECO solutions.
# Usage: enter
#   >> ./run-geco-postprocess
# from the sequence directory (e.g. tss_L08)

# For each step_ directory, compute geco-postprocess-data in that directory.

# Directories
CURRENTDIR=$( pwd )
ALLSTEPS=$( ls adaptive_solver/ | grep '^step_' )

for STEP in $ALLSTEPS
do
    cd $CURRENTDIR/adaptive_solver/$STEP
    geco-postprocess-data
done

# Print a nice message
echo
echo "Postprocessing complete!"
