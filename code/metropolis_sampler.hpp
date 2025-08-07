/*
 * This file is part of FrontUQ-2024-proceedings-likelihood-ratios.
 * Copyright (C) 2025 Emil Loevbak emil.loevbak@kit.edu
 *
 * FrontUQ-2024-proceedings-likelihood-ratios is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FrontUQ-2024-proceedings-likelihood-ratios is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FrontUQ-2024-proceedings-likelihood-ratios.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <random>
#include <vector>
#include <memory>
#include <algorithm>
#include <omp.h>

namespace samplers
{
    template <typename Solver>
    class MetropolisSampler
    {
    private:
        std::shared_ptr<Solver> solver;
        std::vector<double> output;
        double input;
        std::vector<double> mu;
        double sigma;
        double minInput;
        double maxInput;

        std::default_random_engine generator;
        std::uniform_real_distribution<double> uniformDistribution;
        std::vector<double> workVector;

        double generateProposal(double input)
        {
            std::lognormal_distribution<double> logNormalDistribution(std::log(input), 0.25);
            auto proposal = logNormalDistribution(generator);
            // Ensure the proposal is within bounds
            while (proposal < minInput || proposal > maxInput)
            {
                proposal = generateProposal(input);
            }
            return proposal;
        }

    public:
        double proposalRatio(double input, double proposal)
        {
            return proposal / input;
        }
        double likelihoodFunction(std::vector<double> x)
        {
            std::transform(x.begin(), x.end(), mu.begin(), workVector.begin(),
                           std::minus<double>());
            double differenceNorm = std::accumulate(workVector.begin(), workVector.end(), 0.0, [](double a, double b)
                                                    { return a + b * b; });
            return 1.0 / std::pow(sigma / std::sqrt(2 * M_PI), x.size()) * std::exp(-0.5 * differenceNorm / std::pow(sigma, 2.0));
        }
        MetropolisSampler(std::shared_ptr<Solver> solver, double initialInput, std::vector<double> mu, double sigma, double minInput, double maxInput)
            : solver(solver), input(initialInput), mu(mu), sigma(sigma), minInput(minInput), maxInput(maxInput)
        {
            generator = std::default_random_engine(std::random_device{}());
            uniformDistribution = std::uniform_real_distribution<double>(0.0, 1.0);
            output = solver->solve(input);
            workVector.resize(mu.size());
        }

        MetropolisSampler(const MetropolisSampler &other)
            : solver(std::make_shared<Solver>(*other.solver)), input(other.input), mu(other.mu), sigma(other.sigma)
        {
            generator = std::default_random_engine(std::random_device{}());
            uniformDistribution = std::uniform_real_distribution<double>(0.0, 1.0);
            output = solver->solve(input);
            workVector.resize(mu.size());
        }

        double sampleGIMH()
        {
            bool accepted;
            return sampleGIMH(accepted);
        }
        double sampleMCWM()
        {
            bool accepted;
            return sampleMCWM(accepted);
        }
        double sampleGIMH(bool &accepted)
        {
            // Generate new proposal sample
            double proposal = generateProposal(input);

            // Solve the model with the new sample
            auto proposalOutput = solver->solve(proposal);

            // Evaluate the likelihood
            double likelihood = proposalRatio(input, proposal) * likelihoodFunction(proposalOutput) / likelihoodFunction(output);

            // Accept or reject the proposal
            double u = uniformDistribution(generator);
            if (u < likelihood)
            {
                input = proposal;
                output = proposalOutput;
                accepted = true;
            }
            else
            {
                accepted = false;
            }

            // Return the sample
            return input;
        }
        double sampleMCWM(bool &accepted)
        {
            double proposal = generateProposal(input);
            auto proposalOutput = solver->solve(proposal);
            auto outputRecomputation = solver->solve(input);
            double likelihood = proposalRatio(input, proposal) * likelihoodFunction(proposalOutput) / likelihoodFunction(outputRecomputation);
            double u = uniformDistribution(generator);
            if (u < likelihood)
            {
                input = proposal;
                output = proposalOutput;
                accepted = true;
            }
            else
            {
                accepted = false;
            }

            // Return the sample
            return input;
        }
    };
} // namespace samplers