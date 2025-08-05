# FrontUQ-2024-proceedings-likelihood-ratios
Manuscript and experiments for the article titled "An Investigation into the Distribution of Ratios of Particle Solver-based Likelihoods", submitted to the proceedings of FrontUQ 2024.

## Running the code
The repository contains two scripts. The first one `setup_script.sh` builds the software and runs a small scale version of the paper experiments to verify that everything is working correctly. The second one `run_full_simulation.sh` does the same, but runs the full simulation correspondinding to the figures from the paper. This simulation will take some time (potentially in the order of an hour or more on a laptop).

## Dependencies
The C++ code for generating the results relies on Eigen, however the build system downloads this automatically assuming you have a working internet connection. The build is based on CMake and Make. If you do not have these on your system, they can be installed using your favourite package manager. The Python code for plottig requires PyPlot. If you do not have this library installed, you can do so using `pip3 install pyplot`. Otherwise, please reach out with any issues to emil.loevbak@kit.edu.
