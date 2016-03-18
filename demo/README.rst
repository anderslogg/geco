==========
GECo Demos
==========

This directory contains a collection of demo programs that solve the
Vlasov-Poisson, Einstein-Vlasov and Rotating Einstein-Vlasov equations
with axisymmetry using a variety of material model ansatses.

Running the demos
=================

To run the demos, you must -- unless you have already installed the
GECo library and have it in your path -- update the Python path by
running the command

    export PYTHONPATH="$PWD/..:$PYTHONPATH"

Then demos can be run as illustrated by this example:

    python vp_disk.py

Computed solutions can then be found in the `solutions` directory.

Advanced usage
==============

To run a demo using extra high resolution of the computational mesh,
add the --hires flag:

    python vp_disk.py --hires

One may also use the run script provided in this directory to run
demos with extra high resolution. This command will in addition store
the computed solution in a directory named by the demo program and
current timestamp. It will also store a log file from the computation
to that same directory.
