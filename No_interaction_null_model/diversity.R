# obtain, from the table returned by "organize_results()" defined in
# "solve_eqs.R" at a given moment of time, the densities, trait means,
# and trait covariances of each species
# Input:
# - snap: table returned by "organize_results()" defined in "solve_eqs.R"
# Output:
# - a named list with the following entries:
#   * n: an S-component vector; n[i] is the density of species i
#   * m: an SxL matrix; m[i,k] is the kth trait mean of species i
#   * P: an SxLxL array; P[i,k,l] is the (k,l)th covar matrix entry of species i
extract_dynamical_variables <- function(snap) {
  S <- nrow(snap) # number of species
  # solve for number of trait axes, given the number of columns in snap
  L <- select(snap, matches("(m\\d)")) %>% ncol()
  n <- rep(0, S) # species densities
  m <- matrix(0, S, L) # trait means
  P <- array(0, c(S, L, L)) # phenotypic covariances
  for (i in 1:S) {
    n[i] <- snap$n[i]
    m[i,] <- select(snap, matches(paste0("m[1-", L, "]"))) %>%
      slice(i) %>%
      as.numeric()
    Gvec <- select(snap, matches(paste0("G[1-", L, "][1-", L, "]"))) %>%
      slice(i) %>%
      as.numeric()
    Evec <- select(snap, matches(paste0("E[1-", L, "][1-", L, "]"))) %>%
      slice(i) %>%
      as.numeric()
    P[i,,] <- matrix(Gvec + Evec, L, L, byrow=TRUE) # reshape into matrix
    # symmetrize matrix, to eliminate asymmetries from tiny numerical errors
    P[i,,] <- (P[i,,] + t(P[i,,]))/2
  }
  list(n=n, m=m, P=P)
}

# obtain species diversity (Hill number) of order q
# Input:
# - n: vector of species abundances (does not have to be relative)
# - q: order of diversity index (q = 2 is the inverse Simpson index)
# Output:
# - the species diversity of the community
species_diversity <- function(n, q) {
  p <- n/sum(n) # convert abundances to relative frequencies
  sum(p^q)^(1/(1-q))
}

# calculate functional diversity over trait space, given species' densities,
# trait means, and trait covariance matrices
# Input:
# - snap: table returned by "organize_results()" defined in "solve_eqs.R",
#         filtered for one given moment of time
# - q: order of diversity index (q = 2 is the inverse Simpson index)
# - taxis: trait axis, e.g. seq(from, to, length.out), used for each trait axis
# Output:
# - the functional diversity of the community
functional_diversity <- function(snap, q, taxis) {
  vars <- extract_dynamical_variables(snap)
  n <- vars$n # species densities
  m <- vars$m # trait means
  P <- vars$P # trait covariance matrices
  S <- dim(m)[1] # number of species
  L <- dim(m)[2] # number of trait dimensions
  # construct an expression which performs expand.grid on L trait axes
  create_grid <- "expand.grid(taxis" # start of instruction
  # if L>1, do expand.grid L times on taxis
  if (L>1) for (i in 1:(L-1)) create_grid <- paste0(create_grid, ", taxis")
  create_grid <- paste0(create_grid, ")") # closing parenthesis
  grid <- eval(parse(text=create_grid)) %>% as.matrix %>% unname # create grid
  s <- rep(0, nrow(grid)) # vector to store diversity data per grid cell
  # for each species, add up the multivariate normals at each grid cell:
  for (i in 1:length(n)) {
    s <- s + n[i]*dmvnorm(x=grid, mean=m[i,], sigma=matrix(P[i,,], L, L))
  }
  s <- s/sum(s) # to get diversity values, the entries of s should sum to 1
  res <- taxis[2]-taxis[1] # grid resolution of trait space
  (res^L)*(sum(s^q))^(1/(1-q)) # normalized functional diversity of order q
}

# calculate distance-dependent functional diversity
# Input:
# - snap: table returned by "organize_results()" defined in "solve_eqs.R",
#         filtered for one given moment of time
# - q: order of diversity index (q = 2 for inverse Simpson index)
# - xi: decay coefficient weighting phenotype similarity
# - taxis: trait axis, e.g. seq(from, to, length.out), used for each trait axis
# Output:
# - the functional diversity of the community
functional_diversity_distance_dep <- function(snap, q, xi, taxis) {
  vars <- extract_dynamical_variables(snap)
  n <- vars$n # species densities
  m <- vars$m # trait means
  P <- vars$P # trait covariance matrices
  S <- dim(m)[1] # number of species
  L <- dim(m)[2] # number of trait dimensions
  # construct an expression which performs expand.grid on L trait axes
  create_grid <- "expand.grid(taxis" # start of instruction
  # if L>1, do expand.grid L times on taxis
  if (L>1) for (i in 1:(L-1)) create_grid <- paste0(create_grid, ", taxis")
  create_grid <- paste0(create_grid, ")") # closing parenthesis
  grid <- eval(parse(text=create_grid)) %>% as.matrix %>% unname # create grid
  s <- rep(0, nrow(grid))
  for (i in 1:length(n)) {
    s <- s + n[i]*dmvnorm(x=grid, mean=m[i,], sigma=matrix(P[i,,], L, L))
  }
  s <- s/sum(s)
  Z <- exp(-as.matrix(dist(grid, diag=TRUE, upper=TRUE))/xi)
  (s%*%(Z%*%s)^(q-1))^(1/(1-q))
}
