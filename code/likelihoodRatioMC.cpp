#include "heat_equation_grid.hpp"
#include "heat_equation_mc.hpp"
#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include "metropolis_sampler.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

int main(int argc, char **argv)
{
    if (argc != 7)
    {
        std::cerr << "Usage: " << argv[0] << " <true_state> <current_state> <new_state> <measurement_sigma> <number_of_particles> <number_of_samples>" << std::endl;
        return 1;
    }
    double true_state;
    double current_state;
    double new_state;
    double measurement_sigma;
    size_t number_of_particles;
    size_t number_of_samples;

    try
    {
        true_state = std::stod(argv[1]);
        current_state = std::stod(argv[2]);
        new_state = std::stod(argv[3]);
        measurement_sigma = std::stod(argv[4]);
        number_of_particles = std::stoul(argv[5]);
        number_of_samples = std::stoul(argv[6]);
    }
    catch (std::exception &e)
    {
        std::cerr << "Error in parsing command line argument: " << e.what() << std::endl;
        return 1;
    }
    assert(true_state >= 0.0);
    assert(current_state >= 0.0);
    assert(new_state >= 0.0);
    assert(measurement_sigma > 0.0);
    assert(number_of_particles > 0);
    assert(number_of_samples > 0);

    double domain_length = 10.0;
    double minInput = 0.0001;
    double maxInput = 1000.0;
    size_t number_of_cells = 100;
    double dt = 0.1;
    double end_time = 10.0;
    std::default_random_engine generator(std::random_device{}());

    solvers::HeatEquationGrid solver(domain_length, number_of_cells, dt, end_time);
    auto true_solution = solver.solve(true_state);

    std::cout << "True solution: ";
    for (const auto &val : true_solution)
    {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::vector<double> sum(4, 0.0);
    std::vector<double> sum_squared(4, 0.0);
    std::vector<double> output_sum(number_of_cells, 0.0);
    std::vector<double> output_sum_squared(number_of_cells, 0.0);

    std::shared_ptr<solvers::HeatEquationMonteCarlo> solver_mc = std::make_shared<solvers::HeatEquationMonteCarlo>(domain_length, number_of_cells, dt, end_time, number_of_particles);

    for (size_t i = 0; i < number_of_samples; ++i)
    {
        // Each thread needs its own generator and distribution
        std::normal_distribution<double> distribution(0.0, measurement_sigma);

        std::vector<double> noise(true_solution.size());
        for (size_t j = 0; j < true_solution.size(); ++j)
        {
            noise[j] = distribution(generator);
        }
        std::vector<double> perturbed_solution(true_solution.size());
        std::transform(true_solution.begin(), true_solution.end(), noise.begin(), perturbed_solution.begin(),
                       std::plus<double>());
        samplers::MetropolisSampler sampler_mc(solver_mc, true_state, perturbed_solution, measurement_sigma, minInput, maxInput);

        auto current_output = solver_mc->solve(current_state);
        double current_likelihood = sampler_mc.likelihoodFunction(current_output);
        auto new_output = solver_mc->solve(new_state);
        double new_likelihood = sampler_mc.likelihoodFunction(new_output);

        sum[0] += current_likelihood;
        sum[1] += new_likelihood;
        sum[2] += new_likelihood / current_likelihood;
        sum[3] += std::min(new_likelihood / current_likelihood, 1.0);
        sum_squared[0] += current_likelihood * current_likelihood;
        sum_squared[1] += new_likelihood * new_likelihood;
        sum_squared[2] += (new_likelihood / current_likelihood) * (new_likelihood / current_likelihood);
        sum_squared[3] += std::min(new_likelihood / current_likelihood, 1.0) * std::min(new_likelihood / current_likelihood, 1.0);

        std::transform(new_output.begin(), new_output.end(), output_sum.begin(), output_sum.begin(), std::plus<double>());
        std::transform(new_output.begin(), new_output.end(), output_sum_squared.begin(), output_sum_squared.begin(), [](double value, double acc)
                       { return acc + value * value; });
    }

    std::cout << "Output sum: ";
    for (const auto &val : output_sum)
    {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    std::cout << "Output sum squared: ";
    for (const auto &val : output_sum_squared)
    {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    std::cout << "Sum: ";
    for (const auto &val : sum)
    {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    std::cout << "Sum squared: ";
    for (const auto &val : sum_squared)
    {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::vector<double> average(4);
    std::transform(sum.begin(), sum.end(), average.begin(), [number_of_samples](double value)
                   { return value / number_of_samples; });
    std::vector<double> variance(4);
    std::transform(sum_squared.begin(), sum_squared.end(), average.begin(), variance.begin(),
                   [number_of_samples](double value, double avg)
                   { return (value / number_of_samples) - (avg * avg); });
    std::cout << "Average values: " << average[0] << ", " << average[1] << ", " << average[2] << ", " << average[3] << std::endl;
    std::cout << "Variance values: " << variance[0] << ", " << variance[1] << ", " << variance[2] << ", " << variance[3] << std::endl;
    std::vector<double> output_variance(number_of_cells);
    std::transform(output_sum.begin(), output_sum.end(), output_sum_squared.begin(), output_variance.begin(),
                   [number_of_samples](double sum, double sum_sq)
                   { return (sum_sq / number_of_samples) - (sum / number_of_samples) * (sum / number_of_samples); });
    double output_variance_mean = std::accumulate(output_variance.begin(), output_variance.end(), 0.0) / output_variance.size();
    double output_variance_min = *std::min_element(output_variance.begin(), output_variance.end());
    double output_variance_max = *std::max_element(output_variance.begin(), output_variance.end());
    std::cout << "Output variance statistics: Mean=" << output_variance_mean << ", Min=" << output_variance_min << ", Max=" << output_variance_max << std::endl;

    std::ostringstream filename;
    filename << "output/likelihood_results_true_" << true_state
             << "_current_" << current_state
             << "_new_" << new_state
             << "_sigma_" << measurement_sigma
             << "_particles_" << number_of_particles
             << "_samples_" << number_of_samples
             << ".csv";
    std::ofstream outfile(filename.str());
    if (outfile.is_open())
    {
        for (int i = 0; i < 4; ++i)
        {
            outfile << average[i] << (i < 3 ? ", " : "\n");
        }
        for (int i = 0; i < 4; ++i)
        {
            outfile << variance[i] << (i < 3 ? ", " : "\n");
        }
        outfile << output_variance_mean << ", " << output_variance_min << ", " << output_variance_max << "\n";
        outfile.close();
    }
    else
    {
        std::cerr << "Failed to open output file: " << filename.str() << std::endl;
    }
    return 0;
}