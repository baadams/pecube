# Pecube-D-inrusion
This repository contains the codes used to produce a suite codes used to analyze output from Pecube-D-intrusion. These codes were written by Byron A. Adams. 

*Directory contents*

bring_the_heat.m - 
This MATLAB function, bring_the_heat, processes output from a single Pecube-D 1D thermal/exhumation simulation and turns it into diagnostic plots and a temperature-depth movie. It reads the main Pecube.in file plus Pecube output files for temperature histories, exhumation rates, predicted ages, and observed cooling ages. This is a post-processing and visualization script for checking how a Pecube model evolves thermally, how fast material is exhumed, how isotherms move, and whether the model’s predicted cooling ages are consistent with observed thermochronology.

In broad terms, the code does six things:
It parses optional inputs such as target_depth, iso_temp, movie frame rate, whether to save figures, and video format.
It reads the Pecube model setup from Pecube.in, including model dimensions, number of time steps, model depth, vertical nodes, thermal conductivity, and surface temperature.
It enters the output folder and extracts exhumation rates, temperature-depth profiles, sample time-temperature histories, and predicted thermochronologic ages.
It reconstructs a sample’s depth history through time using the exhumation-rate history, and compares predicted thermochronologic ages against observed ages from the CSV input file.
It calculates thermal metrics through time, especially the shallow geothermal gradient to a chosen target_depth, the depth of a chosen isotherm such as 70 °C, and a cumulative heat/energy proxy based on conductivity, geothermal gradient, and model area.
It generates figures showing temperature-depth evolution, exhumation rate, isotherm depth, geothermal gradient, cumulative advected energy, and the sample’s time-temperature-depth path. It also writes a movie of the evolving geotherm and optionally saves figures as .eps, .png, and .fig.
  
Well_fit.m - 
This MATLAB function, well_fit, is a post-processing tool for Pecube-D Monte Carlo results. Its purpose is to identify not just the single “best” exhumation history, but the set of exhumation histories that are statistically consistent with the observed cooling ages, and optionally with an independent estimate of total exhumed thickness. This code asks: “Which Pecube exhumation histories reproduce all of my observed thermochronology within uncertainty, and which of those also make geological sense given independent exhumation-depth constraints?” It is therefore a filtering and uncertainty-visualization tool, rather than a simple best-fit picker.

In broad terms, the code does the following:

It reads optional inputs, including the chi-squared tolerance, integration time for calculating total exhumation, percentage uncertainties, an independent exhumation-depth constraint, and plotting/saving options.
It reads Pecube.in to extract the Monte Carlo time slices, then reads the observed cooling ages and uncertainties from the ages_XXXX.csv file.
It loads the full Monte Carlo output files: ages_all_XXXX.txt, containing predicted cooling ages, and rates_all_XXXX.txt, containing the exhumation-rate histories that generated those ages.
It calculates a chi-squared misfit for each predicted age against each observed age. Any simulation that exceeds the specified tolerance for any chronometer is rejected.
If an integration_time is supplied, it integrates each surviving exhumation-rate history to calculate total exhumed thickness. If cull_data is set to 'y', it further filters the results to keep only histories consistent with the independent exhumation-depth constraint.
It then plots the distributions of acceptable predicted ages against the observed age constraints, the distribution of acceptable exhumed thicknesses if relevant, and the mean acceptable exhumation-rate history with uncertainty.
If depth culling is enabled, it also plots rejected versus accepted exhumation histories, helping show how the independent depth constraint narrows the possible thermal/exhumation histories.
  
# Attribution
If you use or modify these codes and use the results in a publication, please cite the corresponding publication: Cooper, F.J. Adams, B.A., S.I.R. Dahlström, T.A. Ehlers, M.C. van Soest, K.V. Hodges, B.R. Jicha, B.S. Singer, J. Cortes Yañez, R.J. Perkins h (2026). Reassessing long-term exhumation rates in magmatic terranes. Earth and Planetary Science Letters, 691, doi.org/10.1016/j.epsl.2026.120217.

# Error Reporting and Feature Request
If you encounter a bug or have a suggestion for a new feature / improvement the preferred method of communication is to use the 'Issues' function built into GitHub. You can also email Byron [byron.adams 'at' ucl.ac.uk]. If you encounter an issue that you know how to fix and are comfortable with how git works, please feel free to fork the code and submit a pull request with fixes and improvements.
