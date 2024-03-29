# Diatom Tolerance Metrics to Identify Total Phosphorus as a Candidate Cause of Aquatic Life Impairment in Connecticut, USA Freshwater Streams

Mary E. Becker, Timothy J. Becker, Christopher J. Bellucci

https://doi.org/10.1016/j.ecolind.2018.05.046

## ABSTRACT

Anthropogenic phosphorus is a major driver of cultural eutrophication in rivers and streams, leading to numerous water quality impairments, including detrimental shifts in biological communities. Phosphorus has not been identified as a cause of aquatic life impairment in the State of Connecticut (CT), USA, rivers and streams. That is because phosphorus effects on aquatic life are complex, varying spatially and temporally, and often have indirect effects on biological communities typically used for water quality assessment, such as macroinvertebrates and fish. Biological tolerance metrics can be useful in identifying biological impairments due to pollutants that do not lend themselves to conventional toxicity testing, like phosphorus, by providing a measure of the sensitivity of aquatic organisms to anthropogenic disturbance over time. Diatom species tolerances to phosphorus have been derived at national and regional scales in the USA, but not specifically for CT. National scale studies often have the advantage of utilizing larger datasets to derive tolerances over a wide range of environmental conditions, however, developing tolerances specific to a region or for CT may better capture localized conditions. Our study aims to identify diatom species tolerance value metrics suited to aiding aquatic life assessments in CT. We developed diatom tolerance metrics using two different methods that combined responses of individual diatom species along the observed phosphorus gradient using data collected in CT. We then compared the existing national and regional diatom tolerance metrics to the CT tolerance metrics. Our results found the best performing metrics for use in CT aquatic life assessments were CT tolerance values derived using a generalized additive modeling approach. These metrics were CT specific, discriminated well between high and low levels of phosphorus concentrations, and had a greater response to phosphorus than alternative ecological gradients (chloride, pH, and water temperature) that also affect diatom species composition. These results show that diatom tolerance metrics for phosphorus can be effectively used in a weight of evidence approach to identify phosphorus as a cause of aquatic life impairment in CT. All source code and data for this project is freely available and open source at: https://github.com/marybecker/DiatomTPMetrics.

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
