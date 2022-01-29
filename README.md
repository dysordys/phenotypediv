## Code and data for the manuscript "The evolution of trait variance creates a tension between species diversity and functional diversity"

### System requirements

* Operating systems: CentOS 7, Ubuntu 19.10, Ubuntu 20.04, MacOS 10.13, Windows 10.
* Software dependencies: R (has been tested with R 4.1.0).
* Required R packages:
  - `deSolve`: Integrating differential equations.
  - `mvtnorm`: Multivariate normal distributions.
  - `tidyverse`: Efficient data manipulation and plotting.
* Required non-standard hardware: none.
* Typical installation time on a normal desktop computer: no appreciable time if R is already installed. Otherwise, it is the installation time of R and the above three packages.


### Repository contents

The directory `Galapagos_snails` has the files:

* `snaildata.csv`: the raw snail data, in comma-separated value format
* `readme.txt`: short description of the data file `snaildata.csv`
* `vegzonetotals.csv`: the number of distinct host plant species found on the humid and arid zones of each island, in comma-separated value format
* `snail_analysis.R`: code to obtain and plot species- and functional diversity in the snail dataset
* `snail_distance_dep.R`: code to obtain and plot species- and functional diversity in the snail dataset, with distance-dependent diversity metric

The directory `Model_simulation` has the files:

* `model.R`: code that integrates the model equations and displays their results, for a single parameterization
* `solve_eqs.R`: functions used for numerically solving the model equations, and to re-format the final results in a tidy table
* `plot.R`: functions for plotting the results
* `diversity.R`: functions to obtain species- and functional diversity
* `run_replicates.R`: script to help directly replicate the numerical results

The directory `No_variance_evoluton_null_model` has the same files as `Model_simulation`. Everything is identical, *except* that trait variances do not evolve: the right-hand sides of the differential equations for the entries of the genetic covariance matrix are all set to zero. (Note: this implementation is not optimized for efficiency. Instead, it emphasizes that the sole difference between this null model and the original model is the lack of trait variance evolution.)

The directory `No_interaction_null_model` has the files:

* `diversity.R`: functions to obtain species- and functional diversity
* `no_interaction.R`: script to perform all null model simulations, as well as obtain species- and functional diversity calculations for each

The directory `General_1D_model` is like `Model_simulations`, but it only works for one-dimensional trait spaces. In return, it implements a fully general interaction structure between phenotypes, instead of restricting analysis to Gaussian interaction kernels. The code here implements two alternative interaction structures: one that is asymmetric, and one that is nonlocal (hierarchical competition with a competition-mortality tradeoff).
