
require(deSolve)   ## for solving differential equations
require(mvtnorm)   ## for multivariate normal distributions
require(tidyverse) ## for efficient data manipulation and plotting

source("solve_eqs.R") ## functions for integrating the ODEs
source("plot.R") ## functions for plotting results
source("diversity.R") ## functions for calculating functional diversity

set.seed(57310) ## set random seed for reproducibility
S <- 3 ## number of species
L <- 1 ## number of trait dimensions
w <- 0.15 ## competition width
W <- diag(rep(w, L)^2, nrow=L, ncol=L)/2
E <- array(0, c(S, L, L)) ## species' environmental trait covariances
for (i in 1:S) E[i,,] <- diag(runif(L, 0.005, 0.008)^2, nrow=L, ncol=L)
bshape <- "b_quadratic" ## shape of intrinsic growth function
params <- list(W=W, E=E, bshape=bshape) ## list of parameters

ninit <- rep(1, S) ## initial species densities
minit <- matrix(runif(S*L, -0.5, 0.5), S, L) ## initial species trait means
Ginit <- array(0, c(S, L, L)) ## initial genetic covariances
for (i in 1:S) {
    O <- qr.Q(qr(matrix(rnorm(L^2), L, L))) ## create random orthogonal matrix
    B <- diag(runif(L, 0.01, 0.05)^2, nrow=L, ncol=L) ## random PC axis lengths
    Ginit[i,,] <- O%*%B%*%t(O) ## merge decomposition
}
ic <- c(ninit, minit, Ginit) ## merge initial conditions in one vector

## solve ODEs for 1e10 time units, with the parameters defined above
dat <- run_dynamics(ic, seq(0, 1e10, by=1e7), params) %>%
    organize_results(params) ## put results in tidy table

## in 1D and 2D, plot the final trait distributions
if (L==1) plot_snapshot(dat %>% filter(time==max(time)), limits=c(-0.6, 0.6, NA))
if (L==2) plot_snapshot_2D(dat %>% filter(time==max(time)))

## obtain data for the final time step
final <- dat %>%
    filter(time==max(time)) %>% ## only final community state
    arrange(-n) %>% ## order rows in decreasing order of population density
    mutate(n=ifelse(n<=1e-6, 0, n)) ## species below threshold density are extinct

## calculate species diversity
spdiv <- species_diversity(final$n, 2)
## calculate functional diversity of order 2 (inverse Simpson index):
funcdiv <- functional_diversity(final, 2, seq(-1, 1, l=101))

## print results on screen
print(final) ## final community state
print(spdiv) ## species diversity
print(funcdiv) ## functional diversity
