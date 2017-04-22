# Species distribution modelling for Callipogon relictus

This repository contains scripts and data needed to reproduce the species distribution models published in 
Kim S, de Medeiros BAS, et al. 2017. Journal of Biogeography (submitted)

Contents of this folder:

1. relictus_distribution_final.txt
File with locality data for Callipogon relictus

2. SDM_Callipogon.R

Runs all models and projections using biomod2, in an utm map projection

3. plot_projections_lonlat.R

Transform the results of previous scripts back to latlon to enable plotting pretty maps

4. plots_expansion_loss.R

Plots graphs of area expansion and loss (figure 5)

5. figure_ensemble_projection.R

Plots figure 4 (ensemble projection to present and one future scenario)

6. supplementary_plots.R

Plots all supplementary figures

7. maxent/

Directory with maxent program

8. islands_exclude/

Directory with kml file containg polygon to exclude islands.

## Notes

Initially, we excluded Japan and other islands from figure 4 and supplementary maps, bu later decided to include them. Files 5 and 6, with suffix w_islands, are modifications of these scripts to include islands in the plots.

*.sh files are slurm scripts to run some of the scripts above in Harvard's Odyssey Cluster. The others were run in interactive sessions.
