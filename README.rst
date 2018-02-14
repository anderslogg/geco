============================================
GECo - Gothenburg Einstein solver Collection
============================================

------------
Introduction
------------

GECo is a collection of solvers for the Einstein-Vlasov equations
developed at Chalmers University of Technology / University of
Gothenburg.

------------
Installation
------------

GECo is installed like any other Python module:

    [sudo] pip install .

Alternatively, if you want to run a code against GECo without
installing it, enter the GECo directory (this directory) and
run the following command:

    export PYTHONPATH=`pwd`:$PYTHONPATH

------------
Dependencies
------------

GECo depends on the very latest version of FEniCS.

-----
Usage
-----

GECo can be used either as a programmable Python library (see demos in
the demo subdirectory) or via its graphical user interface (command 'geco').

-------
License
-------

To be decided.

-------
Authors
-------

GECo is written by

   Ellery Ames <ellery@kth.se>
   Anders Logg <logg@chalmers.se>

with valuable theoretical / modeling / algorithmic input from

   Håkan Andreasson <hand@chalmers.se>
