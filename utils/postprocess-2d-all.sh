#!/bin/bash

for f in step*; do
    echo ""
    echo "Postprocessing $f"
    echo "--------------------------"
    cd $f
    geco-postprocess-2d-density
    geco-postprocess-2d-ergoregion
    cd ..
done
