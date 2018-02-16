#!/bin/bash

for f in step*; do
    echo ""
    echo "Postprocessing $f"
    echo "--------------------------"
    cd $f
    #geco-postprocess-3d-density-torus
    geco-postprocess-3d-ergoregion-torus
    cd ..
done
