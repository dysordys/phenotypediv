# obtain b_i, g_i, and Q_i, assuming the intrinsic growth function is quadratic
# Input:
# - m: an SxL matrix; m[i,k] is the kth trait mean of species i
# - P: an SxLxL array; P[i,k,l] is the (k,l)th covar matrix entry of species i
# Output:
# - a list with three named entries:
#   * b: vector with S entries
#   * g: an SxL matrix
#   * Q: an SxLxL array
b_quadratic <- function(m, P) {
  S <- dim(m)[1]
  L <- dim(m)[2]
  b <- rep(0, S)
  g <- matrix(0, S, L)
  Q <- array(0, c(S, L, L))
  for (i in 1:S) {
    trP <- as.numeric(ifelse(L==1, P[i,,], sum(diag(P[i,,])))) # trace(P_i)
    mm <- sum(m[i,]^2) # m_i %*% m_i
    b[i] <- 1-(trP+mm)/0.5^2
    g[i,] <- -2*m[i,]/0.5^2
    Q[i,,] <- -2*diag(L)/0.5^2
  }
  list(b=b, g=g, Q=Q)
}

# obtain b_i, g_i, and Q_i, assuming the intrinsic growth function is quartic
# Input:
# - m: an SxL matrix; m[i,k] is the kth trait mean of species i
# - P: an SxLxL array; P[i,k,l] is the (k,l)th covar matrix entry of species i
# Output:
# - a list with three named entries:
#   * b: vector with S entries
#   * g: an SxL matrix
#   * Q: an SxLxL array
b_quartic <- function(m, P) {
  S <- dim(m)[1]
  L <- dim(m)[2]
  b <- rep(0, S)
  g <- matrix(0, S, L)
  Q <- array(0, c(S, L, L))
  for (i in 1:S) {
    trP <- as.numeric(ifelse(L==1, P[i,,], sum(diag(P[i,,])))) # trace(P_i)
    trP2 <- as.numeric(sum(diag(P[i,,]%*%P[i,,]))) # trace(P_i %*% P_i)
    mm <- sum(m[i,]^2) # m_i %*% m_i
    b[i] <- 1-(trP^2+2*trP2+2*trP*mm+4*m[i,]%*%P[i,,]%*%m[i,]+mm^2)/0.5^4
    g[i,] <- -4*(trP*m[i,]+mm*m[i,]+2*P[i,,]%*%m[i,])/0.5^4
    Q[i,,] <- -4*(trP*diag(L)+mm*diag(L)+2*(m[i,]%o%m[i,])+2*P[i,,])/0.5^4
  }
  list(b=b, g=g, Q=Q)
}

# apply smoothed step function to array: values less than 0 are set to 0;
# values between 0 and 1 are set to 10*n^3-15*n^4+6*n^5; and values larger
# than 1 are set to 1
# Input:
# - n: vector or arbitrary array of values
# Output:
# - array of values with smoothed step function applied to them
cutoff <- function(n) {
  ifelse(n<1, (1*(n>0))*(n*n*n*(10+n*(-15+6*n))), 1)
}

# right-hand side of governing equations
# Input:
# - time: moment of time (this is here for compatibility reasons; the
#         system of equations does not actually depend on time explicitly)
# - state: the dynamical state variables (densities, trait means, and trait
#          covariance matrices), packaged into a single vector
# - pars: list of parameters, with the following elements:
#         $W: an LxL matrix, the covariance matrix of the competition kernel
#         $E: an SxLxL array; E[i,k,l] is the (k,l)th entry of species i's
#             environmental covariance matrix
#         $bshape: the name of the function (as a string) to be executed for
#                  calculating b and g, which depend on the shape of the
#                  intrinsic growth function. For instance, it can be
#                  "b_quadratic" or "b_quartic"
# Output:
# - time derivatives of the state variables coerced into a single vector
eqs <- function(time, state, pars) {
  S <- dim(pars$E)[1] # number of species
  L <- dim(pars$E)[2] # number of trait dimensions
  n <- state[1:S] # population densities
  m <- matrix(state[(S+1):(S*(L+1))], S, L) # trait means
  G <- array(state[(S*(L+1)+1):(S*(L^2+L+1))], c(S, L, L)) # genetic covars
  P <- G + pars$E # total phenotypic covar matrices
  alpha <- matrix(0, S, S)
  beta <- array(0, c(S, S, L))
  for (i in 1:S) {
    for (j in 1:S) {
      A <- solve(P[i,,]+P[j,,]+pars$W) # (P_i + P_j + W)^(-1)
      dm <- m[j,]-m[i,] # difference of trait means
      alpha[i,j] <- sqrt(det(pars$W)*det(A))*exp(-dm%*%A%*%dm/2)
      beta[i,j,] <- (A%*%dm)*alpha[i,j]
    }
  }
  alphan <- alpha%*%n
  betan <- matrix(0, S, L)
  for (k in 1:L) betan[,k] <- beta[,,k]%*%n
  p <- match.fun(pars$bshape)(m, P)
  dndt <- (p$b-alphan)*n*cutoff(n/(1e-8))
  dmdt <- matrix(0, S, L)
  dGdt <- array(0, c(S, L, L))
  for (i in 1:S) dmdt[i,] <- G[i,,]%*%(p$g[i,]-betan[i,])
  list(c(dndt, dmdt, dGdt))
}

# put the matrix containing the solution of the equations in tidy format
# Input:
# - dat: dynamical data returned by the function "run_dynamics()" below
# - pars: list of parameters, with the following elements:
#         $W: an LxL matrix, the covariance matrix of the competition kernel
#         $E: an SxLxL array; E[i,k,l] is the (k,l)th entry of species i's
#             environmental covariance matrix
#         $bshape: the name of the function (as a string) to be executed for
#                  calculating b and g, which depend on the shape of the
#                  intrinsic growth function. For instance, it can be
#                  "b_quadratic" or "b_quartic"
# Output:
# - a tibble, organized as follows: each row is a species at a given moment
#   of time, and the columns are:
#   * time: the moment of time
#   * species: the ID of the species in the given row (a positive integer)
#   * G11, G12, ...: entries of the genetic covariance matrix
#   * m1, m2, ...: entries of the mean trait vector
#   * n: population density
#   * bshape: shape of the intrinsic growth function
#   * E11, E12, ...: entries of the environmental covariance matrix
#   * W11, W12, ...: entries of the competition kernel's covariance matrix
#   * Sinit: the number of species at the initial moment of time
organize_results <- function(dat, pars) {
  S <- dim(pars$E)[1] # number of species
  L <- dim(pars$E)[2] # number of trait dimensions
  names(dat)[1] <- "time"
  ind <- 2
  for (i in 1:S) {
    names(dat)[ind] <- paste0("n_", i)
    ind <- ind + 1
  }
  for (k in 1:L) {
    for (i in 1:S) {
      names(dat)[ind] <- paste0("m", k, "_", i)
      ind <- ind + 1
    }
  }
  for (l in 1:L) {
    for (k in 1:L) {
      for (i in 1:S) {
        names(dat)[ind] <- paste0("G", k, l, "_", i)
        ind <- ind + 1
      }
    }
  }
  dat <- dat %>%
    gather("variable", "val", 2:ncol(dat)) %>%
    separate(variable, c("type", "species"), sep="_") %>%
    mutate(species=as.integer(species)) %>%
    spread(type, val) %>%
    mutate(bshape=pars$bshape)
  for (k in 1:L) {
    for (l in 1:L) {
      dat[paste0("E", k, l)] <- pars$E[dat$species, k, l]
    }
  }
  for (k in 1:L) {
    for (l in 1:L) {
      dat[paste0("W", k, l)]=pars$W[k, l]
    }
  }
  mutate(dat, Sinit=as.integer(S))
}

# solve the ODEs of the eco-evolutionary model
# Input:
# - state: the initial conditions, coerced into a single vector
# - tseq: the time sequence for which results should be returned
# - pars: list of parameters, with the following elements:
#         $W: an LxL matrix, the covariance matrix of the competition kernel
#         $E: an SxLxL array; E[i,k,l] is the (k,l)th entry of species i's
#             environmental covariance matrix
#         $bshape: the name of the function (as a string) to be executed for
#                  calculating b and g, which depend on the shape of the
#                  intrinsic growth function. For instance, it can be
#                  "b_quadratic" or "b_quartic"
# Output:
# - the solution converted to a tibble
run_dynamics <- function(state, tseq, pars) {
  ode(func=eqs, y=state, parms=pars, times=tseq) %>% # solve ODEs
    as.data.frame() %>% # convert to data frame (necessary for next step)
    as_tibble() # convert to tibble
}
