# This file is part of FrontUQ-2024-proceedings-likelihood-ratios.
# Copyright (C) 2025 Emil Loevbak emil.loevbak@kit.edu
#
# FrontUQ-2024-proceedings-likelihood-ratios is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FrontUQ-2024-proceedings-likelihood-ratios is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with FrontUQ-2024-proceedings-likelihood-ratios.  If not, see <https://www.gnu.org/licenses/>.

import os
import re
import numpy as np
from itertools import product

import matplotlib.pyplot as plt

# Directory containing the output files
output_dir = "output"

true_vals = ["0.1"] 
current_vals = ["0.1", "0.08"]
new_vals =  ["0.1", "0.08"]
sigma_vals = ["0.025", "0.1"]
particles_vals = ["100", "1000", "10000"]
samples_vals = ["100"]

for true, current, new, samples in product(
    true_vals, current_vals, new_vals, samples_vals
):
    # Skip the case where new_vals=current_vals=0.08
    if current == "0.08" and new == "0.08":
        continue

    # Create figures for each quantity at every iteration of the outermost loop
    figs, axes = plt.subplots(4, 1, figsize=(10, 20), sharex=True)
    if not isinstance(axes, np.ndarray):
        axes = [axes]
    lines_labels = [[] for _ in range(4)]  # To keep track of labels for legend

    for sigma in sigma_vals:
        particles_list = []
        mean_current = []
        mean_proposal = []
        mean_ratio = []
        mean_acceptance = []
        variance_current = []
        variance_proposal = []
        variance_ratio = []
        variance_acceptance = []
        mean_solver_variance = []
        min_solver_variance = []
        max_solver_variance = []

        for particles in particles_vals:
            filename = f"likelihood_results_true_{true}_current_{current}_new_{new}_sigma_{sigma}_particles_{particles}_samples_{samples}.csv"
            filepath = os.path.join(output_dir, filename)

            if not os.path.exists(filepath):
                print(f"Missing file: {filepath}")

            with open(filepath, "r") as f:
                mean_line = f.readline().strip()
                var_line = f.readline().strip()
                mean_vals = list(map(float, mean_line.split(",")))
                var_vals = list(map(float, var_line.split(",")))
                mean_current.append(mean_vals[0])
                mean_proposal.append(mean_vals[1])
                mean_ratio.append(mean_vals[2])
                mean_acceptance.append(mean_vals[3])
                variance_current.append(var_vals[0])
                variance_proposal.append(var_vals[1])
                variance_ratio.append(var_vals[2])
                variance_acceptance.append(var_vals[3])

                # Read variance values from file
                output_variance_line = f.readline().strip()
                output_variance_vals = list(map(float, output_variance_line.split(",")))
                mean_solver_variance.append(output_variance_vals[0])
                min_solver_variance.append(output_variance_vals[1])
                max_solver_variance.append(output_variance_vals[2])

        # Plotting
        quantity_names = ["Current likelihood", "Proposal likelihood", "Likelihood ratio", "Acceptance rate"]

        particles_arr = np.array(list(map(int, particles_vals)))
        means_arr = [
            np.array(mean_current),
            np.array(mean_proposal),
            np.array(mean_ratio),
            np.array(mean_acceptance),
        ]

        label = f"sigma={sigma}"
        for i in range(4):
            axes[i].plot(np.sqrt(max_solver_variance), means_arr[i], marker='o', label=label)
            axes[i].set_ylabel(f"Mean of {quantity_names[i]}")
            axes[i].set_xscale("log")
            if i == 3:
                axes[i].set_ylim(0, 1)
            else:
                axes[i].set_yscale("log")
            axes[i].grid(True)
            lines_labels[i].append(label)

        # Set xlabel only for the last subplot
        axes[-1].set_xlabel("max_solver_sigma")

        # Add legends after all lines are added
        for i in range(4):
            axes[i].legend()

    plt.suptitle(
        f"True={true}, Current={current}, Proposal={new}"
    )
    plt.show()
    # Write out raw plot data as CSV per subplot
    for i, quantity in enumerate(quantity_names):
        csv_filename = (
            f"plotdata_{quantity.replace(' ', '_').lower()}_true_{true}_current_{current}_proposal_{new}_samples_{samples}.csv"
        )
        csv_filepath = os.path.join(output_dir, csv_filename)
        with open(csv_filepath, "w") as csvfile:
            # Write header
            csvfile.write("particles," + ",".join([f"sigma={s}" for s in sigma_vals]) + "\n")
            # Write data rows
            for idx, particles in enumerate(particles_vals):
                row = [particles]
                for sigma in sigma_vals:
                    # Find the index of sigma in sigma_vals
                    sigma_idx = sigma_vals.index(sigma)
                    # Get the corresponding means_arr for this subplot and sigma
                    # means_arr[i] is a list of arrays, one per sigma
                    # But in this code, means_arr[i] is only for the last sigma in the loop
                    # So we need to collect all means for all sigmas
                    # Instead, let's build a 2D array for each subplot
                    pass  # We'll build this below

    # Build 2D arrays for each subplot and sigma, and collect max_solver_variance
    means_by_sigma = [[[] for _ in sigma_vals] for _ in range(4)]
    max_solver_variances = [[[] for _ in sigma_vals] for _ in range(len(particles_vals))]
    for sigma_idx, sigma in enumerate(sigma_vals):
        for p_idx, particles in enumerate(particles_vals):
            filename = f"likelihood_results_true_{true}_current_{current}_new_{new}_sigma_{sigma}_particles_{particles}_samples_{samples}.csv"
            filepath = os.path.join(output_dir, filename)
            if not os.path.exists(filepath):
                means = [np.nan] * 4
                max_var = np.nan
            else:
                with open(filepath, "r") as f:
                    mean_line = f.readline().strip()
                    mean_vals = list(map(float, mean_line.split(",")))
                    means = mean_vals
                    _ = f.readline()  # skip variance line
                    output_variance_line = f.readline().strip()
                    output_variance_vals = list(map(float, output_variance_line.split(",")))
                    max_var = output_variance_vals[2]
            for i in range(4):
                means_by_sigma[i][sigma_idx].append(means[i])
            # Store max_solver_variance for this particles/sigma
            if len(max_solver_variances[p_idx]) < len(sigma_vals):
                max_solver_variances[p_idx].append(max_var)
            else:
                max_solver_variances[p_idx][sigma_idx] = max_var

    # Now write the CSVs
    for i, quantity in enumerate(quantity_names):
        csv_filename = (
            f"plotdata_{quantity.replace(' ', '_').lower()}_true_{true}_current_{current}_proposal_{new}_samples_{samples}.csv"
        )
        csv_filepath = os.path.join(output_dir, csv_filename)
        with open(csv_filepath, "w") as csvfile:
            # Write header: particles, then for each sigma: mean, max_solver_variance
            header = ["max_solver_sigma"]
            for s in sigma_vals:
                header.append(f"meansigma{s}")
                header.append(f"maxsolverstdsigma{s}")
            csvfile.write(",".join(header) + "\n")
            for idx, particles in enumerate(particles_vals):
                row = [particles]
                for sigma_idx in range(len(sigma_vals)):
                    row.append(str(means_by_sigma[i][sigma_idx][idx]))
                    row.append(str(np.sqrt(max_solver_variances[idx][sigma_idx])))
                csvfile.write(",".join(row) + "\n")

    # Write CSV with ratio of expectations of likelihoods
    ratio_csv_filename = (
        f"plotdata_ratio_expectations_true_{true}_current_{current}_proposal_{new}_samples_{samples}.csv"
    )
    ratio_csv_filepath = os.path.join(output_dir, ratio_csv_filename)
    with open(ratio_csv_filepath, "w") as csvfile:
        # Write header: particles, then for each sigma: ratio of means
        header = ["max_solver_sigma"]
        for s in sigma_vals:
            header.append(f"ratiosigma{s}")
            header.append(f"maxsolverstdsigma{s}")
        csvfile.write(",".join(header) + "\n")
        for idx, particles in enumerate(particles_vals):
            row = [particles]
            for sigma_idx in range(len(sigma_vals)):
                if means_by_sigma[0][sigma_idx][idx] != 0:
                    ratio = means_by_sigma[1][sigma_idx][idx] / means_by_sigma[0][sigma_idx][idx]
                else:
                    ratio = 0.0
                row.append(str(ratio))
                row.append(str(np.sqrt(max_solver_variances[idx][sigma_idx])))
            csvfile.write(",".join(row) + "\n")

print("Plots written as CSV files in the output directory.")