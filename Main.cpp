#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <numeric>

// Function to simulate asset price paths using the Black-Scholes model
std::vector<std::vector<double>> simulatePricePaths(double S0, double mu, double sigma, double T, int numSteps, int numPaths) {
    std::vector<std::vector<double>> pricePaths(numPaths, std::vector<double>(numSteps + 1, S0));
    double dt = T / numSteps;
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<> distribution(0.0, 1.0);

    for (int i = 0; i < numPaths; ++i) {
        for (int j = 1; j <= numSteps; ++j) {
            double dW = distribution(generator) * sqrt(dt);
            pricePaths[i][j] = pricePaths[i][j-1] * exp((mu - 0.5 * sigma * sigma) * dt + sigma * dW);
        }
    }
    return pricePaths;
}

// Function to calculate the realized variance within a corridor
double calculateRealizedVariance(const std::vector<std::vector<double>>& pricePaths, double lowerBound, double upperBound) {
    std::vector<double> returns;
    for (const auto& path : pricePaths) {
        for (size_t i = 1; i < path.size(); ++i) {
            if (path[i-1] >= lowerBound && path[i-1] <= upperBound) {
                double ret = log(path[i] / path[i-1]);
                returns.push_back(ret);
            }
        }
    }
    double mean = std::accumulate(returns.begin(), returns.end(), 0.0) / returns.size();
    double sq_sum = std::inner_product(returns.begin(), returns.end(), returns.begin(), 0.0);
    return sq_sum / returns.size() - mean * mean;
}

// Main function to price the corridor variance swap using Monte Carlo simulation
int main() {
    // Parameters
    double S0 = 100.0;     // Initial asset price
    double mu = 0.1;       // Drift
    double sigma = 0.2;    // Volatility
    double T = 1.0;        // Time to maturity in years
    int numSteps = 252;    // Number of time steps (e.g., daily steps)
    int numPaths = 10000;  // Number of Monte Carlo paths
    double lowerBound = 90.0; // Lower bound of the corridor
    double upperBound = 110.0; // Upper bound of the corridor
    double strikeVariance = 0.04; // Strike variance

    // Simulate price paths
    auto pricePaths = simulatePricePaths(S0, mu, sigma, T, numSteps, numPaths);

    // Calculate the realized variance within the corridor
    double realizedVariance = calculateRealizedVariance(pricePaths, lowerBound, upperBound);

    // Calculate the payoff
    double notional = 1000000; // Example notional amount
    double payoff = (realizedVariance - strikeVariance) * notional;

    // Output the result
    std::cout << "Realized Variance within Corridor: " << realizedVariance << std::endl;
    std::cout << "Payoff: " << payoff << std::endl;

    return 0;
}
