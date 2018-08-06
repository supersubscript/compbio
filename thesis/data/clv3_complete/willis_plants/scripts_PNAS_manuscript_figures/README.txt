Before running each script, check the variables at the beginning of the script to see whether any need modifying. 

The config (config_2015_02_16_plantX.py) files contain :
-paths to data extracted from segmentation and tracking
-times of data aquisition, x-y-z resolution of segmentation, and stretching factor (already factored into data extraction)
This file should be compiled as the data is segmented and tracked. Set the "savePath" variable to path for saving (some) figures

The parameter (indices_parameters_2015_02_16_plantX.py) files contain :
-the 3D centre of the SAM at first and last time-point, selected manually
-lists of indices of L1 cells within certain radii of SAM centre at first and last time-point
-the 3D coordinates of the L1 centre, for all time-points
-the mean volumes among L1 cells < 30 um from the SAM centre
This file should be compiled by manually identifying the SAM center at the first and last timepoint, and applying the following scripts:
1. identify_L1_L2_L3_indices_and_avVolumes_in_circular_region.py 
	to extract lists of indices at specified time-points within specified distances from centre, which should then be copied into parameter file
2. trackBackwards_meanVolumeEvolutionAndCentreL1.py, which produces L1centre and meanVolumes, to be copied into parameter file.


Select a particular plant, then radius from centre for which cells to include, then run script:
area_vs_vols_vs_dist_within_fixed_range_newAreas.py, updating config and parameters_indices files and the radius accordingly
to generate figures for 
-division rates
-time evolution violin plots
-cell size distributions
-num. neighbours vs. cell volumes
-distance from centre vs. size scatter plots
-log(volume) vs. log(surface area metric) scatter plots
-phase of light/dark cycle vs. cell volume violin plot


Select a particular plant and radius from centre (usually 30 micron), then run script
divisions_within_fixed_range.py, updating config and parameters_indices files and the radius accordingly
(make sure indicesL1Init in parameters_indices file is for radius exceeding selected radius)
to produce data files of completed cell cycles within 30 micron radius of centre
Div_t0_new_X.dat, each line of file corresponds to a cell birth with format  (j0, m0, d0)
Div_t1_new_X.dat, the same line is for the corresponding cell division (j1,m1,d11, d12)
where cell is born between time[j0] with mother index m0 and time[j0+1] with daughter index d0, and divides between time[j1] and time[j1+1] with mother index m1 and daughter indices d11, d12


Select a radius from centre and a region (match the tag "X" in the Div_t0_new_X.dat file), and what data to include (large cells at birth, cells born in the morning, etc) modifying the "zone" variable accordingly, then run script:
divisionLawPooledPlants.py
to generate figures for pooled plant data on completed cell cycles
-birth size vs. inter-division time (referred to as cell cycle time)
-birth size vs. division size
-birth size vs. increment
-dist. from centre vs. inter-division time
-volumetric birth asymmetry vs. birth size
-time vs. inter-division time


Select a radius from centre and a region (again, match the tag "X" in the Div_t0_new_X.dat file), and choose a multiple of the standard deviation beyond which outliers are excluded, then run script
divisionLawAllPlants.py
to generate figures similar to divisionLawPooledPlants.py but for individual plants

Select a radius from centre and a "zone", also standard deviations for excluding outliers, then run script
sisterStatsPooledPlants.py
to generate figures for pooled plant data on sisters where both sisters complete a cell cycle
-mother's size vs. birth asymmetry
-distance from centre vs. birth, division, and inter-division time asymmetry 
-time vs. birth, division, and inter-division time asymmetry
-birth asymmetry vs. volume asymmetry among sister cell cycles included in data
-birth asymmetry vs. volume compared with non-sister neighbour volume (V - V^ns)/(V + V^ns) where V^ns is mean of volume of non-sister neighbours
-sister vs. sister: size at division, birth, inter-division time
-birth asymmetry (alpha_b) vs. division asymmetry (alpha_d) with outliers removed
-birth asymmetry (alpha_b) vs. increment asymmetry (alpha_delta) with outliers removed
-total vol. sisters at birth vs. total vol. sisters at division
-histogram of time of day at division and birth 


Select a radius from centre and a "zone" (make sure the corresponding results folder exists), then run script
sisterGrowthStatsPooledPlants.py
to generate figures for pooled plant data on sisters where both sisters complete a cell cycle
-phase of day vs. per unit size growth rate
-phase of cell cycle vs. absolute growth rate
-phase of cell cycle vs. per unit size growth rate
-phase of cell cycle vs. per unit size growth rate for symmetric/symmetric divisions, for cells born large/small
-birth size vs. exponential growth rate for symmetric/asymmetric divisions
-normalized difference between cell and non-sister neighbours vs exponential growth rate for symmetric/asymmetric divisions
-birth asymmetry (alpha_b) vs exponential growth rate for intermediate/extreme birth sizes
-birth asymmetry (alpha_b) vs exponential growth rate for similar/different sizes relative to non-sister neighbours
-large sister/small sister size time evolution over cell cycle
-birth asymmetry vs. exponential growth rate averaged over corresponding cell cycle
 

Select a radius from centre and a minimum value for z-coordinate, only cells within this radius and above this z-coordinate will be included in the analysis.
Check that the 'label' matches the plant in imported config and indices_parameters files, then run script
mother_daughter_growth_rates_within_fixed_range.py
to generate
-2D heat map of relative (or exponential) growth rate normalized by relative growth rate of neighbours of mother vs. daughter cell
-scatter plot of relative (or exponential) growth rate normalized by relative growth rate of neighbours of mother vs. daughter cell

Run script nuclei_birth_division_stats_newArea.py to generate data on inter-division times, nuclear volumes, cell volumes for all completed cell cycles within 'radius_out' and with z-coordinates greater than 'zMin'. Enter 'radius_in' and 'region_in' to import corresponding tracking and cell indices data on complete cell cycles; 'radius_in' should be greater than 'radius_out'. 

In script, all_res_birth_div_nuclei.py, select 'radius' and minimum z-coordinate 'zThres'; update all paths as appropriate. Generates
-time since birth vs. nuclear volume
-time since birth vs. ratio of nuclear:cell volume

















