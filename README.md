# Diatom Tolerance Metrics to Identify Total Phosphorus as a Candidate Cause of Aquatic Life Impairment in Connecticut Freshwater Streams

Biological tolerance metrics developed by combining responses of individual diatom species along the observed phosphorus gradient in State of Connecticut, USA.  Individual diatom species responses to TP were examined using generalized additive models and curve classification (Yuan 2004,2006). These metrics were found to discriminate well between different levels of ambient phosphorus concentrations and had a greater response to phosphorus than alternative ecological gradients previously identified as affecting variation in diatom species composition, pH and water temperature.  

* R ver. 3.3.1, 
* gam, ggplot2, grid, plyr, reshape2

## DATA

**SPP_040617:** 
Species occurrence data for each sample used in the analysis <br> <br>
**TP:** 
Total Phosphorus data for each sample used in the analysis <br> <br>
**DiatomMetrics_040517:** 
Combined species responses for each sample with Total Phosphorus groups Low (L) & High (H) <br>
H - % Relative Abundance Tolerant Species <br>
L - % Relative Abundance Sensitive Species <br>
R - TP Index<br> <br>
**DiatomMetricspH:** pH data with L / H groups <br> <br>
**DiatomMetricsJulyTemp:** 
Temperature data with L / H groups 

## ANALYSIS

#### GAM_DiatomSppTP_032816.R:  
Calculates tolerance values for each diatom species
#### GAMRun040517 (RWorkspace):
Tolerance values analysis with results used in paper.  Results could vary slightly due to bootstrap analysis with additional runs.
#### DiscriminationEfficiency_Diatom_LH.R:
Calculates discrimination efficiency for metrics and creates boxplots.<br><br>

*Data Visualization Tool in Development.  Current DRAFT version at: https://ctriverresearch.shinyapps.io/TPAQLapp*
