/**
 * A rudementary MCMC sampler using Metropolis-Hastings.
 * @author: Emil Løvbak
 */

#pragma once

#include <Eigen/Dense>
#include <random>

namespace samplers
{
    template <typename Solver>
    class MetropolisSampler
    {
    private:
        std::shared_ptr<Solver> solver;
        Eigen::VectorXd output;
        double input;
        Eigen::VectorXd mu;
        double sigma;
        double minInput;
        double maxInput;

        std::default_random_engine generator;
        std::uniform_real_distribution<double> uniformDistribution;

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
        double likelihoodFunction(Eigen::VectorXd x)
        {
            return 1.0 / std::pow(sigma / std::sqrt(2 * M_PI), x.size()) * std::exp(-0.5 * std::pow((x - mu).norm() / sigma, 2.0));
        }
        MetropolisSampler(std::shared_ptr<Solver> solver, double initialInput, Eigen::VectorXd mu, double sigma, double minInput, double maxInput)
            : solver(solver), input(initialInput), mu(mu), sigma(sigma), minInput(minInput), maxInput(maxInput)
        {
            generator = std::default_random_engine(std::random_device{}());
            uniformDistribution = std::uniform_real_distribution<double>(0.0, 1.0);
            output = solver->solve(input);
        }

        MetropolisSampler(const MetropolisSampler &other)
            : solver(std::make_shared<Solver>(*other.solver)), input(other.input), mu(other.mu), sigma(other.sigma)
        {
            generator = std::default_random_engine(std::random_device{}());
            uniformDistribution = std::uniform_real_distribution<double>(0.0, 1.0);
            output = solver->solve(input);
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