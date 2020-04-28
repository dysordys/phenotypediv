## Code and data for the manuscript "The evolution of trait variance creates a tension between species diversity and functional diversity"

The folder `Galapagos_snails` has the files:

* `snaildata.csv`: the raw snail data, in comma-separated value format
* `readme.txt`: short description of the data file `snaildata.csv`
* `vegzonetotals.csv`: the number of distinct host plant species found on the humid and arid zones of each island, in comma-separated value format
* `snail_analysis.R`: code to obtain and plot species- and functional diversity in the snail dataset
* `snail_distance_dep.R`: code to obtain and plot species- and functional diversity in the snail dataset, with distance-dependent diversity metric

The folder `Model_simulation` has the files:

* `QGR_model.R`: code that integrates the model equations and displays their results
* `solve_eqs.R`: functions used for numerically solving the model equations, and to re-format the final results in a tidy table
* `plot.R`: functions for plotting the results
* `diversity.R`: functions to obtain species- and functional diversity
* `run_replicates.R`: script to help directly replicate our numerical results
