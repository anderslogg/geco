# Usage

GECo can be used either as a programmable Python library (see demos in
the demo subdirectory) or via its graphical user interface (command 'geco').

The simplest way to get started is by using the demos located in `geco/demos/`. 

## Running the demos
Before running the demos you should complete the set-up steps indicated in the [Installation](./installation.md).

To run the demos navigate to `geco/demos/` and run for example
```bash
python3 ev_torus.py
```
The solution files are stored in `geco/demos/solutions/ev/` by default, and this location can be customized by modifying the `output.solution_directory` parameter (for example `solver.parameters["output"]["solution_directory"] = custom_solution_dir`). 
For additional parameters see the `SolverBase` class. 

The generated output files depend on whether the gravitational interaction is Newtonian or general relativistic. 
For the relativistic demo `ev_torus` and similar, GECo generates the following
output files:
- Gravitational fields NU, MU, BB, WW each in `.xdmf` and `.xml.gz` format
- Spatial density RHO in `.xdmf` and `.xml.gz` format
- Energy-momentum components P00, P03, P33, P11 in `.xdmf` and `.xml.gz` format.

For details on the meaning of these quantities see [On Axisymmetric and Stationary Solutions of the Self-Gravitating Vlasov System](https://iopscience.iop.org/article/10.1088/0264-9381/33/15/155008)

## Postprocessing
Several postprocessing scripts are collected in `geco/bin/`. 
To run postprocessing on a solution, first navigate to the solution directory, before running 
```bash
geco-postprocessing-data
```
if GECo is installed, or 
```bash
python3 ~/geco/bin/geco-postprocessing-data
```
otherwise. 
This will run the postprocessing script to generate additional scalar data from the computed solution. 
Several of the postprocessing scripts rely on this data, so it is a good idea to run this script first. 
Other postprocessing scripts can be used similarly. 
Documentation on the various postrpocessing routines can be found in [Postprocessing](./postprocessing.md).

## New Ansatz functions
GECo is a convenient tool for exploring solutions and their properties corresponding to new distribution ansatz functions. 
To write and use a new ansatz, the following steps should be taken: 
- Define your ansatz as a subclass of the `EVAnsatz` or `VPAnsatz` classes, as in `geco/cppcode/EV-E-Polytropic-L-Polytropic.h`. 
  The ansatz should contain member functions `ansatz`, `init_parameters`, and `read_parameters`.
  The model must be defined in a file with the pattern `cppcode/EV-[MY-MODEL-NAME].h` or `cppcode/VP-MY-MODEL-NAME].h` as appropriate.
- Add your ansatz to the list of `PYBIND11_MODULES` in either `geco/cppcod/EVBindings.h` or `geco/cppcod/VPBindings.h`, following the examples already present. 
  Make sure to `#include` the header file corresponding to your new ansatz class.
  This provides Python bindings.
- Finally use your model in a `run_my_model.py` file as in `model = MaterialModel("MY_MODEL-NAME")`. 
  Reference existing demos for more complete information.