# IDSA

This repository contains the scripts and datasets used to run the Interval Dynamic Survival Analysis (IDSA) method as introduced in our paper "Monitoring SARS-COV-2 Transmission and Prevalence under Repeated Testing"

The contents are as follows:

1. IDSA_closed.R - simulates an epidemic in a single population, fits the IDSA closed model, and generates plots of the model fit and R_t. This requires the plot_ly package. The default settings correspond to those used to produce fig 4 in the paper.
2. IDSA_open.R - simulates an epidemic in two coupled populations, fits the IDSA open model, and generates plots of the model fit and R_t. This requires the plot_ly package. The default settings correspond to those used to produce fig 5 in the paper.
3. IDSA_OSU - Fits the IDSA closed model to data from Ohio State University collected during the Fall 2020 semester and generates plots of the model fit and R_t. This requires the plot_ly package and the data file test_mat_OSU.csv
4. test_mat_OSU_partial.csv - This is partial de-identified data from Ohio State University collected during the Fall 2020 semester. A subset of 10000 individuals was selected uniformly at random from the total Fall 2020 on-campus population. For each unique individual, it gives the first time (if ever) they were known to be positive to SARS-COV-2 and the most recent time before that (if one exists) they were known to be negative for SARS-COV-2. IDSA_OSU reads in data from this file.
