#!/bin/bash

#prefix="2d-ergoregion"
#prefix="2d-density"
#prefix="3d-ergoregion-torus"
prefix="3d-ergoregion-box"
#prefix="3d-density-torus"

echo "<?xml version=\"1.0\"?>"
echo "<VTKFile type=\"Collection\" version=\"0.1\">"
echo "  <Collection>"

COUNTER=0

for f in step*; do
    cd $f

    # Alternative 1: use E0 for the time variable for annotation in Paraview
    #E0=`geco-postprocess-print | grep E0 | awk '{ printf("%.4f\n", $3) }'`
    #echo "    <DataSet timestep=\"$E0\" part=\"0\" file=\"../$f/postprocess-$prefix""000000.vtu\" />"

    # Alternative 2: just label by counter = step
    echo "    <DataSet timestep=\"$COUNTER\" part=\"0\" file=\"../$f/postprocess-$prefix""000000.vtu\" />"
    cd ..
    (( COUNTER++ ))
done

echo "  </Collection>"
echo "</VTKFile>"
