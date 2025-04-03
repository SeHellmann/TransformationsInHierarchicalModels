# TransformationsInHierarchicalModels
This repository contains the materials for the paper "A biased parameter transformation in hierarchical models"


## Complete reproducibility
To completely reproduce the exact results from sctrach (which will take pretty long), do the follwing:
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

Rieskamp, J. (2008). The probabilistic nature of preferential choice. Journal of Experimental Psychology: Learning, Memory, and Cognition, 34(6), 1446–1465. [10.1037/a0013646](https://doi.org/10.1037/a0013646).

Nilsson, H. & Rieskamp, J. & Wagenmakers, E.-J. (2011). Hierarchical Bayesian parameter estimation for cumulative prospect theory. Journal of Mathematical Psychology. 55. 84-93. [10.1016/j.jmp.2010.08.006](https://doi.org/10.1016/j.jmp.2010.08.006). 

*Add Preprint Here*

## Contact details

For comments, remarks, and questions please contact: [sebastian.hellmann\@tum.de](mailto:sebastian.hellmann@tum.de){.email}
