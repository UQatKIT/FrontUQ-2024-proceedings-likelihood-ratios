# FrontUQ-2024-proceedings-likelihood-ratios
Code and experiments for the article titled "An Investigation into the Distribution of Ratios of Particle Solver-based Likelihoods", submitted to the proceedings of FrontUQ 2024.

## Running the code
The repository contains two scripts. The first one `setup_script.sh` builds the software and runs a small scale version of the paper experiments to verify that everything is working correctly. The second one `run_full_simulation.sh` does the same, but runs the full simulation correspondinding to the figures from the paper. This simulation will take some time (potentially in the order of an hour or more on a laptop).

## Dependencies
The build is based on CMake and Make. If you do not have these on your system, they can be installed using your favourite package manager. The Python code for plottig requires PyPlot. If you do not have this library installed, you can do so using `pip3 install pyplot`. Otherwise, please reach out with any issues to emil.loevbak@kit.edu.

## Supported systems
The results in the publication were produced on Ubuntu 24.04 using GCC 13.3.0. The code has also been verified to compile on Clang 18.1.3 on the same system. On Mac, there are known build issues with Apple Clang (accessed under the alias gcc under the default configuration) due to the lack of OpenMP support in the compiler. These issues can be avoided by using Homebrew GCC version 15 before running the tests (see below). Testing of the code on other systems is currently ongoing progress. If it does not work on your system, please open an issue.

### Compiling under Homebrew GCC
Follow these steps to use GCC 15. If you do not have version 15 available, another version will likely work by changing the numbers accordingly.
1. Ensure you have Homebrew GCC installed: `brew install gcc@15`
2. Set the compiler environment variables using `export CC=gcc-15` and `export CXX=g++-15`
3. Run the setup script as specified above.

## Note for reviewers
The figures in the submitted manuscript were generated with the code as it apperas in commit `d72a8dff170ffd04efb96e952ba67c5cad599e40`. That version of the code had a dependency on Eigen for the computation of the reference solution with finite differences (not the Monte Carlo solution). This dependency caused compilation issues on Mac. To improve reproducibility, we have since removed this dependency, substituting our own solver for the matrix system in the diffusion equation. We have since re-run the experiments with the latest commit, reproducing the same qualitative figures. Even though there are no qualitative changes in the results, we intend to update the figures upon revision for the sake of correctness. The arXiv paper [2508.05303](https://www.arxiv.org/abs/2508.05303) already contains the updated figures. Here the biggest difference can be seen on the rightmost datapoints in Figure 1c, i.e., computing the ratio of two high-variance estimates based on few Monte Carlo particles. This behavior is not unexpected. Until we have revised the paper, we leave this note here for the sake of transparency.
