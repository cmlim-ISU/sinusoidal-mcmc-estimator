# sinusoidal-mcmc-estimator
Estimation of Frequency and Phase of Sinusoidal Signal using MCMC

# Sinusoidal Signal Parameter Estimation using MCMC

This repository contains code and documentation for estimating the frequency and phase of a sinusoidal signal using Markov Chain Monte Carlo (MCMC) methods with Stan.

## Project Overview

This project implements Bayesian inference to estimate parameters of a sinusoidal signal:

```
y = α * sin(ω * t + γ) + ε
```

where:
- α (alpha): Amplitude of the sine wave
- ω (omega): Angular frequency
- γ (gamma): Phase shift
- ε (epsilon): Random noise (error)

## Files in this Repository

- `CaiMin_544_data.csv`: Input data containing time points and observations
- `project_544-00.stan`: Stan model specification for Bayesian inference
- `projectR.R`: R script for data preparation, MCMC sampling, and result visualization
- `project.rds`: Saved Stan fit object containing MCMC results

## Model Description

The Stan model implements:

1. **Data Section**: Defines inputs including observations and time points
2. **Parameters Section**: Defines model parameters (alpha, omega, gamma, sigma, epsilon)
3. **Transformed Parameters**: Computes the expected sine wave values
4. **Model Section**: Specifies priors and likelihood functions
5. **Generated Quantities**: Creates posterior predictive samples

### Parameter Priors

- α (alpha): Uniform(0, 10)
- ω (omega): Uniform(0, 10)
- γ (gamma): Uniform(-π, π)
- σ (sigma): Normal(0, 1)

## Results

The MCMC sampling produces posterior distributions for all parameters. The model effectively recovers the true sine wave from noisy observations, as demonstrated in the visualization code.

## Usage

1. Ensure you have R installed with the following packages:
   - rstan
   - bayesplot
   - rstanarm
   - ggplot2
   - mcmcse
   - dplyr
   - plyr
   - corrplot
   - PerformanceAnalytics
   - GGally
   - viridis

2. Run the R script:
   ```R
   source("projectR.R")
   ```

3. The script will:
   - Load and prepare the data
   - Run the MCMC sampling using the Stan model
   - Generate visualizations of the posterior distributions
   - Create diagnostic plots (trace plots, correlation plots)
   - Calculate summary statistics and credible intervals

## Visualization Outputs

The code generates multiple visualizations:
- Parameter trace plots to assess convergence
- Posterior distribution plots
- Model fit visualizations comparing observed data to inferred sine wave
- Correlation plots between parameters
- Posterior predictive checks

## References

This project demonstrates Bayesian parameter estimation using MCMC for sinusoidal signal analysis, useful in applications like signal processing, time series analysis, and noise filtering.

