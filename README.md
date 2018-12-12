SSM-SeaScallops: R code for simulating and fitting the alternative state space model of Yin et al. (2018)
---------------------------------------------------------------------------------------------------------

These separate R scripts provide functions to simulate data according to the alternative state space model (SSM) introduced in Yin et al. (2018) for the assessment of the Bay of Fundy sea scallops fishery (Scallop Production Area 4 SPA4 as labelled by Fisheries and Oceans Canada). These scripts also provide functions to fit the model according to a Laplace-approximated maximum likelihood estimation method which relies on the R package Template Model Builder (TMB). Details about the model and the comparison to the reference model can be found in Yin et al. (2018).

Updates can be found at https://github.com/williamaeberhard/SSM-SeaScallops.

Any requests/comments/bug reports should be sent to william.aeberhard@gmail.com.

### Contents

Files contained in this repository:

* SS2v4_main.r: an R script loading all required libraries, compiling the TMB C++ template, simulating data according to the model, fitting the model, and outputting some figure as in Yin et al. (2018);
* SS2v4simul.r: an R script that creates the SS2v4simul function which simulates data according to the SSM;
* SS2v4.cpp: a TMB C++ template file coding the joint negative log-likelihood of the SSM;
* this README file.

### Version History

This is SSM-SeaScallops version 0.1. This is the initial release.

### References

Yin, Y., Aeberhard, W. H., Smith, S. J., and Mills Flemming, J. (2018). Identifiable state-space models: A case study of the Bay of Fundy sea scallop fishery. Canadian Journal of Statistics, In press.

