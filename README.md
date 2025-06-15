# TransformationsInHierarchicalModels
This repository contains the materials for the paper "A biased parameter transformation in hierarchical models"


## Complete reproducibility
To completely reproduce the exact results from scratch (which will take pretty long), do the following:
- remove all `.RData` files, 
- install JAGS (version 4.3.1)
- open the RProject, install the `renv` package, and call `renv::activate()` and `renv::restore()`
- set in line 8 of `Main_Script.R`: `REDOALLANALYSIS <- TRUE` and run the script

## Reproduce results without doing MCMC sampling again
Install all packages (or use `renv::restore()`), and source the `Main_Script.R`. 

## Structure of the project
- Folder *Rieskamp_2008_data* contains the raw data (behavioral data and gambles) from Rieskamp (2008)
....
 
### References

Nilsson, H. & Rieskamp, J. & Wagenmakers, E.-J. (2011). Hierarchical Bayesian parameter estimation for cumulative prospect theory. Journal of Mathematical Psychology. 55. 84-93. [10.1016/j.jmp.2010.08.006](https://doi.org/10.1016/j.jmp.2010.08.006). 

Rieskamp, J. (2008). The probabilistic nature of preferential choice. Journal of Experimental Psychology: Learning, Memory, and Cognition, 34(6), 1446–1465. [10.1037/a0013646](https://doi.org/10.1037/a0013646).

Pachur, T., Mata, R., & Hertwig, R. (2017). Who dares, who errs? disentangling cognitive and motivational roots of age differences in decisions under risk. Psychological Science, 28 (4), 504–518. [10.1177/0956797616687729](https://doi.org/10.1177/0956797616687729)

*Add Preprint Here*

## Contact

For comments, remarks, and questions please contact: [sebastian.hellmann\@tum.de](mailto:sebastian.hellmann@tum.de){.email}
