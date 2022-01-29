require(deSolve)
require(tidyverse)


# Apply smoothed step function to array. Values less than 0 are set to 0; values
# between 0 and 1 are set to 10*n^3-15*n^4+6*x^5; values larger than 1 are
# set to 1. The function aids numerical integration of ODEs.
# Input:
# - n: vector or arbitrary array of values
# - a: cutoff threshold
# Output:
# - array of values with smoothed step function applied to them
cutoff <- function(n) {
  ifelse(n<1, (1*(n>0))*(n*n*n*(10+n*(-15+6*n))), 1)
}

# Per capita growth rates
# Input:
# - n: vector of abundances
# - m: vector of trait means
# - s: vector of total phenotypic standard deviations
# - z: vector with sampling points along the trait axis
# - b: intrinsic growth as a function of trait value
# - a: interaction kernel as a function of trait value difference
# Output:
# - vector of per capita growth rates at each z
pgr <- function(n, m, s, z, b, a) {
  S <- length(n) # number of species
  Z <- length(z) # number of grid points along the trait axis
  sump <- rep(0, Z) # allocate memory for sum_j n_j p_j(z)
  for (j in 1:S) sump <- sump+n[j]*dnorm(z, m[j], s[j])
  zdiff <- outer(z, z, FUN = `-`) # difference matrix of trait positions
  dz <- z[2]-z[1] # difference between two adjacent trait positions
  an <- (a(zdiff) %*% sump) * dz # int a(z,z') sump(z') dz'
  matrix(c(b(z)-an), S, Z, byrow = TRUE)
}

# Right hand side of ODE system
# Input:
# - time: moment of time
# - state: state variables (n, m) coerced into a single vector
# - pars: list of parameters, with the following elements:
#         $b: a function, giving phenotypes' intrinsic growth rates
#         $a: the competition kernel, as a function of trait difference
#         $E: species' environmental phenotypic variances
#         $dz: trait axis grid size
# Output:
# - time derivatives dn/dt and dm/dt, coerced into a vector and then a list
eqs <- function(time, state, pars) {
  S <- length(pars$E) # number of species
  n <- state[1:S] # species abundances
  m <- state[(S+1):(2 * S)] # species trait means
  G <- state[(2*S+1):(3 * S)] # species genetic variances
  P <- G + pars$E # total phenotypic variances = genetic + environmental
  s <- sqrt(P) # total phenotypic standard deviations
  h2 <- G / P # heritabilities
  limits <- range(c(m-4*s, m+4*s)) # limits of integration
  z <- seq(limits[1], limits[2], by = pars$dz) # trait axis
  r <- pgr(n, m, s, z, pars$b, pars$a) # per capita growth rates
  dndt <- rep(0, S) # allocate memory for abundance equations
  dmdt <- rep(0, S) # allocate memory for trait mean equations
  dGdt <- rep(0, S) # allocate memory for trait variance equations
  for (j in 1:S) {
    p <- dnorm(z, m[j], s[j]) # trait distribution of species j
    dndt[j] <- n[j]*as.numeric(p%*%r[j,])*pars$dz*cutoff(n[j]/(1e-6))
    dmdt[j] <- h2[j]*as.numeric(p%*%(r[j,]*(z-m[j])))*pars$dz
    dGdt[j] <- (h2[j]^2/2)*as.numeric(p%*%(r[j,]*((z-m[j])^2-P[j])))*
      pars$dz*cutoff(G[j]/(1e-10))
  }
  list(c(dndt, dmdt, dGdt) * pars$dz)
}

# Organize simulation results into tidy tibble
# Input:
# - sol: output produced by the function ode
# - pars: list of parameters, with the following elements:
#         $b: a function, giving phenotypes' intrinsic growth rates
#         $a: the competition kernel, as a function of trait difference
#         $E: species' environmental phenotypic variances
#         $dz: trait axis grid size
# Output:
# - a tibble w/ columns time, species, n (abundance), m (trait mean), G, E, and s
organize_results <- function(sol, pars) {
  S <- length(pars$E) # number of species
  as_tibble(as.data.frame(sol)) %>% # convert to tibble
    rename_with(~paste0("n_", 1:S), 1+1:S) %>% # density column names: n_1, ..., n_S
    rename_with(~paste0("m_", 1:S), S+1+1:S) %>% # trait mean columns: m_1, ..., m_S
    rename_with(~paste0("G_", 1:S), 2*S+1+1:S) %>% # genetic covars.:  G_1, ..., G_S
    pivot_longer(cols=!"time", names_to="variable", values_to="val") %>%
    separate(variable, c("type", "species"), sep="_", convert=TRUE) %>%
    pivot_wider(names_from=type, values_from=val) %>%
    mutate(E=pars$E[species], s=sqrt(G+E))
}
