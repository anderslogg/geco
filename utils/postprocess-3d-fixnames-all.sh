#!/bin/bash

# Fix names in files to enable time series plotting in Paraview.
# Before this fix, the fields were named f_25 and then changed
# to f_27 or similar. This is not needed after fixing (some of)
# the scripts to rename the FEniCS Function before saving.

A="postprocess-3d-density-torus000000.vtu"
B="postprocess-3d-ergoregion-torus000000.vtu"

for f in step*; do
    echo ""
    echo "Renaming in $f"
    echo "--------------------------"
    cd $f
    #sed -i s/f_25/RHO/g $A
    #sed -i.bak s/f_27/RHO/g $A
    #sed -i.bak s/f_46/GTT/g $B
    #sed -i.bak s/f_48/GTT/g $B
    #sed -i.bak s/RHO/GTT/g $B
    grep Scalars $A
    grep Scalars $B
    cd ..
done
