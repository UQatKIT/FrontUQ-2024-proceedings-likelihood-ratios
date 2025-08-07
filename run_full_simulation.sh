#!/bin/bash
set -e

# Make build directory if it doesn't exist
mkdir -p build
# Change to build directory and build the project
cd build
cmake ../code
make

echo "Build complete. Now running simulation..."
cd ..
mkdir -p output
for measurement_noise in 0.01 0.025 0.05 0.1 0.25 0.5 1; do
    for particle_count in 100 1000 10000 100000 1000000; do
        echo "Running with measurement noise: $measurement_noise and particle count: $particle_count"
        echo "0.1 to 0.1"
        build/frontuq_likelihood_ratios 0.1 0.1 0.1 $measurement_noise $particle_count 1000
        echo "0.1 to 0.08"
        build/frontuq_likelihood_ratios 0.1 0.1 0.08 $measurement_noise $particle_count 1000
        echo "0.08 to 0.1"
        build/frontuq_likelihood_ratios 0.1 0.08 0.1 $measurement_noise $particle_count 1000
    done
done

echo "The full-scale program has completed successfully."
echo "Removing the build directory..."
rm -rf build

echo ""
echo "You can plot the results using the provided Python script 'code/plotLikelihoodRatiosFullScale.py' from this directory."