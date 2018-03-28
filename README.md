# Diatom Tolerance Metrics to Identify Total Phosphorus as a Candidate Cause of Aquatic Life Impairment in Connecticut Freshwater Streams

Biological tolerance metrics developed by combining responses of individual diatom species along the observed phosphorus gradient in State of Connecticut, USA.  Individual diatom species responses to TP were examined using generalized additive models and curve classification (Yuan 2004,2006). These metrics were found to discriminate well between different levels of ambient phosphorus concentrations and had a greater response to phosphorus than alternative ecological gradients previously identified as affecting variation in diatom species composition, pH and water temperature.  

* R ver. 3.3.1, 
* gam, ggplot2, grid, plyr, reshape2

## DATA

**SPP:** <br>
Species occurrence data for each sample used to develop species models with GAMs <br> <br>
**SPP_RelAbund:**<br>
Species relative abundance data for each sample used to develop species tolerances with indicator species analysis and weighted averaging<br> <br>
**TP:** <br>
Total Phosphorus data for each sample used in the analysis <br> <br>
**Diatom Metrics Data:** <br>
Combined species responses for each sample with Total Phosphorus groups Low (L) & High (H) <br>
H - % Relative Abundance Tolerant Species <br>
L - % Relative Abundance Sensitive Species <br>
R - TP Index<br> <br>
DiatomMetrics_GAM: Diatom metrics calculated using GAM tolerance values for sites the calibration dataset<br>
DiatomMetrics_ANSP: Diatom metrics calculated using National tolerance values for sites the calibration dataset<br>
DiatomMetrics_ANSPRegional:  Diatom metrics calculated using Regional tolerance values for sites the calibration dataset<br>
DiatomMetrics_IndVal:  Diatom metrics calculated using IndVal/WA tolerance values for sites the calibration dataset<br>
DiatomMetrics_TESTGAM:  Diatom metrics calculated using GAM tolerance values for sites the test dataset<br>
DiatomMetrics_TESTANSP: Diatom metrics calculated using National tolerance values for sites the test dataset<br>
DiatomMetrics_TESTANSPRegional: Diatom metrics calculated using Regional tolerance values for sites the test dataset<br>
DiatomMetrics_TESTCHEM:  pH and Chloride data used to test alternative gradients<br>
DiatomMetrics_TESTJTEMP:  Temperature data used to test alternative gradients<br> <br>

## SCRIPTS

#### GAM: <br> 
Method to develop tolerance values for each diatom species using GAMs and curve shape classification
#### INDValWA:<br>
Method to develop tolerance values for each diatom species using indicator species analysis and weighted averaging
#### DiscriminationEfficiency_Diatom_LH_10102017.R:<br>
Calculates discrimination efficiency for metrics and creates boxplots.<br><br>
