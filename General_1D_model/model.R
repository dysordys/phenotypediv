source("./solve_eqs.R")
source("./plot.R")


set.seed(54325) # set random seed (for reproducibility)

S <- 3 # number of species

params <- list(
  b=function(z) 1-exp(-z-0.5),
  a=function(z) (tanh(z/0.15)+1)/2,
  E=runif(S, 0.005, 0.008)^2, # environmental variances
  dz=0.005
)

ninit <- rep(0.1, S)
minit <- runif(S, -0.5, 0.5)
Ginit <- runif(S, 0.05, 0.1)^2
ic <- c(ninit, minit, Ginit)
tseq <- seq(0, 1e6, l=1001)

sol <- ode(func=eqs, y=ic, parms=params, times=tseq) %>%
  organize_results(params)

plot_snapshot(filter(sol, time==max(time)), c(NA, NA, NA), 1001, params$b)
