
require(deSolve)   ## for solving differential equations
require(mvtnorm)   ## for multivariate normal distributions
require(tidyverse) ## for efficient data manipulation and plotting

source("solve_eqs.R") ## functions for integrating the ODEs
source("plot.R")      ## functions for plotting results
source("diversity.R") ## functions for obtaining species- and functional diversity

Ss <- 2:10 ## initial numbers of species
ws <- c(0.1, 0.15) ## competition widths
replicates <- 1:1 ## replicates; increase number after the colon to do more
bshapes <- c("b_quadratic") ## shape(s) of intrinsic growth function
Ls <- 1:1 ## trait dimension(s)

## set up factorial numerical experiment, varying the initial number of species,
## number of trait dimensions, the shape of the intrinsic growth function,
## the competition width, and replicate ID
experiment <- expand_grid(S=Ss, L=Ls, bshape=bshapes, w=ws, repl=replicates) %>%
    mutate(spdiv=0, funcdiv=0) ## add columns for future diversity calculations

## do numerical experiment for each row
for (r in 1:nrow(experiment)) {
    set.seed(57300+10*experiment$repl) ## set random seed (for reproducibility)
    S <- experiment$S[r] ## initial number of species
    L <- experiment$L[r] ## number of trait dimensions
    bshape <- experiment$bshape[r] ## shape of intrinsic growth function
    w <- experiment$w[r] ## competition width
    W <- diag(rep(w, L)^2, nrow=L, ncol=L)/2 ## competition matrix
    E <- array(0, c(S, L, L)) ## species' environmental trait covariances
    for (i in 1:S) E[i,,] <- diag(runif(L, 0.005, 0.008)^2, nrow=L, ncol=L)
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
    dat <- run_dynamics(ic, seq(0, 1e10, l=1001), params) %>%
        organize_results(params) %>% ## put results in tidy table
        filter(time==max(time)) %>% ## keep only final state
        mutate(n=ifelse(n<=1e-6, 0, n)) ## species below threshold are extinct
    ## species diversity and functional diversity:
    experiment$spdiv[r] <- species_diversity(dat$n, q=2)
    experiment$funcdiv[r] <- functional_diversity(dat, q=2, seq(-1, 1, l=101))
    ## track numerical experiment number on screen
    cat(paste(r, "/", nrow(experiment), "\n"))
}

## plot results
ggplot(experiment) +
    aes(x=spdiv, y=funcdiv, colour=factor(w)) +
    geom_point() +
    scale_x_continuous(name="species diversity") +
    scale_y_log10(name="functional diversity") +
    scale_colour_manual(values=cpal, name="competition width")
