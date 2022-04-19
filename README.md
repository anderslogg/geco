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

    [sudo -H] pip install [--upgrade] .

Alternatively, if you want to run a code against GECo without
installing it, enter the GECo directory (this directory) and
run the following command:

    export PYTHONPATH=`pwd`:$PYTHONPATH

Tip for working with different directories inside a Docker container:

    alias make='pushd . && cd ~/shared/geco && sudo -H pip install --upgrade . && popd'

------------
Dependencies
------------

GECo depends on FEniCS version 2017.2.0.

-----
Usage
-----

GECo can be used either as a programmable Python library (see demos in
the demo subdirectory) or via its graphical user interface (command 'geco').


-----------
Citing GECo
-----------

If you use GECo in your work please cite the following paper: 

@article{GECo2016,
    author = {Ames, Ellery and Andr{\'e}asson, H{\aa}kan and Logg, Anders},
    title = {{On axisymmetric and stationary solutions of the self-gravitating  Vlasov system}},
    journal = {Class. Quantum Grav.},
    year = {2016},
    volume = {33},
    number = {15},
    pages = {155008},
    doi = {10.1088/0264-9381/33/15/155008}
}


-------
License
-------

Copyright 2019 Anders Logg, Ellery Ames, Håkan Andréasson

GECo is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

GECo is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with GECo. If not, see <https://www.gnu.org/licenses/>.

-------
Authors
-------

GECo is written by

   Ellery Ames <ellery.ames@humboldt.edu>
   Anders Logg <logg@chalmers.se>

with valuable theoretical / modeling / algorithmic input from

   Håkan Andreasson <hand@chalmers.se>

and contributions to postprocessing capabilities by 

    Eric Malekos. 
