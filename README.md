# Estimation of heterogeneous agent models using macro and micro data

Matlab code for full-information Bayesian inference in heterogeneous agent models using both (i) macro time series data and (ii) repeated cross sections of micro data

**Reference:**
[Liu, Laura](https://laurayuliu.com/), and [Mikkel Plagborg-Møller](https://scholar.princeton.edu/mikkelpm) (2020), "Full-Information Estimation of Heterogeneous Agent Models Using Macro and Micro Data", https://scholar.princeton.edu/mikkelpm/het_agents

**Acknowledgements:**
We build on the excellent Dynare code kindly made available by [Thomas Winberry](http://www.thomaswinberry.com/research/index.html) (see also [Winberry, QE 2018](https://qeconomics.org/ojs/index.php/qe/article/view/617))

**Requirements:**
[Dynare](https://www.dynare.org/) version 4.6.1 or higher

Tested in: Matlab R2020a on Windows 10 PC (64-bit) with Dynare 4.6.1

## Contents

**[program](program):** Matlab routines
- [run_mcmc_hh.m](program/run_mcmc_hh.m): simulate and estimate heterogeneous household model
- [run_mcmc_firm.m](program/run_mcmc_firm.m): simulate and estimate heterogeneous firm model
- [plot_mcmc.m](program/plot_mcmc.m): plot estimation output
- [run_likelihood_hh.m](program/run_likelihood_hh.m): compute likelihood functions (for various observables) in heterogeneous household model
- [plot_likelihood.m](program/plot_likelihood.m): plot likelihood functions

**[program/functions](program/functions):** general functions for MCMC, likelihood evaluation, simulations, and plotting
- [likelihood/loglike_compute.m](program/functions/likelihood/loglike_compute.m): main function for numerically unbiased likelihood estimate

**[program/hh_model](program/hh_model):** files specific to the heterogeneous household model
- [dynare](program/hh_model/dynare): sub-folder with Dynare model files adapted from Winberry (2018)
- [auxiliary_functions/likelihood/likelihood_micro.m](program/hh_model/auxiliary_functions/likelihood/likelihood_micro.m): micro likelihood function

**[program/firm_model](program/firm_model):** files specific to the heterogeneous firm model
- [dynare](program/firm_model/dynare): sub-folder with Dynare model files adapted from Winberry (2018)
- [auxiliary_functions/likelihood/likelihood_micro.m](program/firm_model/auxiliary_functions/likelihood/likelihood_micro.m): micro likelihood function
