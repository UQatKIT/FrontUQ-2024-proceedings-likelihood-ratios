#!/bin/bash

# Make build directory if it doesn't exist
mkdir -p build
# Change to build directory and build the project
cd build
cmake ../code
make

echo "Build complete. Now running a small test program..."
cd ..
mkdir -p output
for measurement_noise in 0.1 1; do
    for particle_count in 1000 10000; do
        echo "Running with measurement noise: $measurement_noise and particle count: $particle_count"
        echo "0.1 to 0.1"
        build/frontuq_likelihood_ratios 0.1 0.1 0.1 $measurement_noise $particle_count 100
        echo "0.1 to 0.08"
        build/frontuq_likelihood_ratios 0.1 0.1 0.08 $measurement_noise $particle_count 100
        echo "0.08 to 0.1"
        build/frontuq_likelihood_ratios 0.1 0.08 0.1 $measurement_noise $particle_count 100
    done
done

echo "Generating plots for the small test program..."
python3 code/plotLikelihoodRatiosTest.py
echo "Plots written as CSV files in the output directory."
echo "To run the full-scale simulation, run the script 'run_full_simulation.sh'. This will likely take several hours to complete."
echo "The code makes use of OpenMP, so it is reccomended to set the environment variable OMP_NUM_THREADS to the number of threads you want to use."
