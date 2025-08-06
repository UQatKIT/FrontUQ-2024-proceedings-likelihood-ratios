# FrontUQ-2024-proceedings-likelihood-ratios
Code and experiments for the article titled "An Investigation into the Distribution of Ratios of Particle Solver-based Likelihoods", submitted to the proceedings of FrontUQ 2024.

## Running the code
The repository contains two scripts. The first one `setup_script.sh` builds the software and runs a small scale version of the paper experiments to verify that everything is working correctly. The second one `run_full_simulation.sh` does the same, but runs the full simulation correspondinding to the figures from the paper. This simulation will take some time (potentially in the order of an hour or more on a laptop).

## Dependencies
The build is based on CMake and Make. If you do not have these on your system, they can be installed using your favourite package manager. The Python code for plottig requires PyPlot. If you do not have this library installed, you can do so using `pip3 install pyplot`. Otherwise, please reach out with any issues to emil.loevbak@kit.edu.

## Supported systems
The code should work out of the box on regular Linux distributions and has been tested on Ubuntu 24.04 using both GCC and Clang. Testing of the code on other systems is currently ongoing progress. If it does not work on your system, please open an issue.

## Note for reviewers
The figures in the submitted manuscript were generated with the code as it apperas in commit `d72a8dff170ffd04efb96e952ba67c5cad599e40`. This version of the code had a dependency on Eigen for the computation of the reference solution with finite differences (not the Monte Carlo solution). This dependency caused compilation issues on Mac. To improve reproducibility, we have since removed this dependency, substituting our own solver for the matrix system in the diffusion equation. We have since re-run the experiments with the latest commit, reproducing the same qualitative figures. Even though there are no qualitative changes in the results, we intend to update the figures upon revision for the sake of correctness. In the mean time, we leave this note here.
