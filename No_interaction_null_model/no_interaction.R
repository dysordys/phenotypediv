library(tidyverse)
library(mvtnorm)

source("./diversity.R")


b_quadratic <- function(m, P) {
  S <- dim(m)[1]
  L <- dim(m)[2]
  b <- rep(0, S)
  for (i in 1:S) {
    trP <- as.numeric(ifelse(L==1, P[i,,], sum(diag(P[i,,])))) # trace(P_i)
    mm <- sum(m[i,]^2) # m_i %*% m_i
    b[i] <- 1-(trP+mm)/0.5^2
  }
  return(b)
}

b_quartic <- function(m, P) {
  S <- dim(m)[1]
  L <- dim(m)[2]
  b <- rep(0, S)
  for (i in 1:S) {
    trP <- as.numeric(ifelse(L==1, P[i,,], sum(diag(P[i,,])))) # trace(P_i)
    trP2 <- as.numeric(sum(diag(P[i,,]%*%P[i,,]))) # trace(P_i %*% P_i)
    mm <- sum(m[i,]^2) # m_i %*% m_i
    b[i] <- 1-(trP^2+2*trP2+2*trP*mm+4*m[i,]%*%P[i,,]%*%m[i,]+mm^2)/0.5^4
  }
  return(b)
}

organize_results <- function(n, m, pars) {
  S <- dim(m)[1]
  L <- dim(m)[2]
  dat <- pivot_wider(enframe(c(n, m)))
  ind <- 1
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
  dat <- dat %>%
    gather("variable", "val", 1:ncol(dat)) %>%
    separate(variable, c("type", "species"), sep="_") %>%
    mutate(species=as.integer(species)) %>%
    spread(type, val) %>%
    mutate(bshape=pars$bshape)
  for (k in 1:L) {
    for (l in 1:L) {
      dat[paste0("G", k, l)] <- pars$G[dat$species, k, l]
    }
  }
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
  dat <- dat %>% mutate(Sinit=as.integer(S))
  return(dat)
}

eqs_nocomp <- function(S, L, bshape, Glower, Gupper, Elower, Eupper, w, rseed) {
  set.seed(rseed) # set random seed
  E <- array(0, c(S, L, L)) # species' environmental trait covariances
  for (i in 1:S) E[i,,] <- diag(runif(L, Elower, Eupper)^2, L, L)
  W <- diag(rep(w, L)^2, L, L)/2 # competition width matrix
  m <- matrix(runif(S*L, -0.5, 0.5), S, L) # initial species trait means
  G <- array(0, c(S, L, L)) # initial genetic covariances
  for (i in 1:S) {
    O <- qr.Q(qr(matrix(rnorm(L^2), L, L))) # create random orthogonal matrix
    B <- diag(runif(L, Glower, Gupper)^2, L, L) # random lengths of principal axes
    G[i,,] <- O%*%B%*%t(O) # merge decomposition into the genetic covariance
  }
  n <- match.fun(bshape)(m, E + G)
  params <- list(W=W, E=E, G=G, bshape=bshape, w=w) # list of parameters
  organize_results(n, m, params) %>% # tidy up results
    mutate(Glower=Glower, Gupper=Gupper, Elower=Elower, Eupper=Eupper,
           rseed=as.integer(rseed)) %>% # add Glower,Gupper,Elower,Eupper,rseed
    unite(Gscen, c(Glower, Gupper), sep="_") %>% # Glower and Gupper into Gscen
    unite(Escen, c(Elower, Eupper), sep="_") # Elower and Eupper into Escen
}


# this will run all replicates, which might take too long on a single machine;
# one way to control for this is to only analyze a manageable subset of the
# rows of the table at a time
expand_grid(
  S=2:25,
  L=1:3,
  bshape=c("b_quadratic", "b_quartic"),
  Gscen=c("0.01_0.05", "0.05_0.1"),
  Escen=c("0.005_0.008", "0.015_0.018", "0.025_0.028"),
  w=c(0.1, 0.15, 0.25, 0.3, 0.4, 0.45),
  rseed=as.integer(100*(0:9)+57313)
) %>%
  filter((L==1 & w %in% c(0.1, 0.15)) |
           (L==2 & w %in% c(0.25, 0.3)) |
           (L==3 & w %in% c(0.4, 0.45))) %>%
  separate(Gscen, c("Glower", "Gupper"), sep="_") %>%
  separate(Escen, c("Elower", "Eupper"), sep="_") %>%
  mutate(across(c(Glower, Gupper, Elower, Eupper), as.numeric)) %>%
  mutate(sol=pmap(., eqs_nocomp)) %>%
  mutate(spdiv=map_dbl(sol, ~species_diversity(.$n, 2))) %>%
  mutate(funcdiv=map_dbl(sol, functional_diversity, 2, seq(-1, 1, l=101)))
