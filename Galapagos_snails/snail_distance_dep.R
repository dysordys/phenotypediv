require(tidyverse) # for efficient data manipulation and plotting
require(mvtnorm) # for multivariate normal distributions


# parameters controlling distance-dependent functional diversity metric:
q  <- 2 # Hill exponent; set to 2 for Rao's quadratic entropy
xi <- 1 # distance parameter of the diversity metric


# calculate functional diversity over 2D trait space
# Input:
# - n: vector of species abundances, with as many entries as species
# - m1: m1[i] is the 1st entry of species i's mean trait vector
# - m2: m2[i] is the 2nd entry of species i's mean trait vector
# - P11: P11[i] is the (1,1)th entry of species i's covariance matrix
# - P12, P21, P22: as above, but for the corresponding matrix entries
# - q: Hill exponent (=2 for inverse Simpson index)
# - xi: decay coefficient weighting phenotype similarity
# - taxis: trait axis, i.e. seq(from, to, length.out), for both traits
# Output:
# - the functional diversity of the community
get_div <- function(n, m1, m2, P11, P12, P21, P22, q, xi, taxis) {
  m <- matrix(c(m1, m2), length(n), 2)
  P <- array(c(P11, P21, P12, P22), c(length(n), 2, 2))
  ztab <- unname(as.matrix(expand.grid(t1=taxis, t2=taxis)))
  p <- rep(0, nrow(ztab))
  for (i in 1:length(n)) p <- p+n[i]*dmvnorm(x=ztab, mean=m[i,], sigma=P[i,,])
  p <- p/sum(p)
  Z <- exp(-as.matrix(dist(ztab, diag=TRUE, upper=TRUE))/xi)
  func_div <- (p%*%(Z%*%p)^(q-1))^(1/(1-q))
  return(func_div)
}

# maximum likelihood estimates of mean and covariance of binormal distribution
# Input:
# - dat: matrix, where dat[i,k] is individual i's trait component k
# Output:
# - a vector with 6 entries: the two entries of the mean vector (m1, m2) and
#   the four entries of the covariance matrix (P11, P12, P21, P22)
mle_binorm <- function(dat) {
  N <- nrow(dat) # number of data points
  me <- unname(colSums(dat)/N) # fitted mean
  Pe <- c(0, 0)%o%c(0, 0) # initialize covariance matrix
  for (i in 1:N) Pe <- Pe + (dat[i,]-me)%o%(dat[i,]-me) # obtain cov matrix
  Pe <- Pe/(N-1) # normalize by degrees of freedom minus 1 for unbiased estimate
  return(c(me, Pe)) # return fitted mean and covariance, coerced into a vector
}


# load data with number of plants in humid/arid zones per island
snailveg <- read_csv("vegzonetotals.csv")

# load morphological data
snail <- read_csv("snaildata.csv") %>%
  # trait 1: centroid size; trait 2: PC1 of shape (explains 80% of variance)
  rename(t1="Centroidsize", t2="shape1") %>%
  # standardize trait data by subtracting mean & dividing by std dev
  mutate(t1=(t1-mean(t1))/sd(t1), t2=(t2-mean(t2))/sd(t2)) %>%
  # remove three outlier satellite islands
  filter(!(island %in% c("CH", "ED", "GA"))) %>%
  # throw out achatellinus species, which has only 2 sampled individuals
  filter(!(species=="achatellinus")) %>%
  # replace NAs with 0s in vegetation zone columns
  replace_na(list(vegzoneHumid=0, vegzoneTransition=0, vegzoneArid=0)) %>%
  # convert numerical vegetation zone info into "humid" / "arid"
  mutate(zone=case_when(
    vegzoneHumid==1 ~ "humid",
    vegzoneArid==1 ~ "arid",
    TRUE ~ "none")
  ) %>%
  filter(zone %in% c("humid", "arid")) %>%
  # keep only those columns that will be needed
  select(species, island, zone, t1, t2)


# lower and upper limits of data in both trait directions, used to define taxis
lim <- snail %>% select(t1, t2) %>% range(na.rm=TRUE)
# increasing length and/or resolution of trait axis does not change the results
taxis <- seq(3*lim[1], 3*lim[2], l=51)

snail_diversity <- snail %>%
  # obtain number of individuals per species, plus trait mean and covar matrix
  group_by(species, island, zone) %>%
  summarise(n=n(), fit=list(mle_binorm(as.matrix(tibble(t1=t1, t2=t2)))),
            .groups="drop") %>%
  # new columns for mean trait (m1,m2) and covar matrix (P11,P12,P21,P22)
  mutate(
    m1=map_dbl(fit, `[`(1)),
    m2=map_dbl(fit, `[`(2)),
    P11=map_dbl(fit, `[`(3)),
    P12=map_dbl(fit, `[`(4)),
    P21=map_dbl(fit, `[`(5)),
    P22=map_dbl(fit, `[`(6)),
  )

# obtain species- and functional diversity measures
snail_diversity <- snail_diversity %>%
  group_by(island, zone) %>%
  # for each island and zone, obtain:
  summarise(S=length(unique(species)), # species richness
            ntot=sum(n), # number of sampled individuals in local community
            # species diversity (change "2" below to change type of index)
            spdiv=sum((n/sum(n))^2)^(1/(1-2)),
            # functional diversity
            funcdiv=get_div(n, m1, m2, P11, P12, P12, P22, q, xi, taxis),
            .groups="drop") %>%
  # add data on plant species richness in each zone on each island
  mutate(plants=ifelse(zone=="humid",
                       snailveg$HumidTotal[match(island, snailveg$island)],
                       snailveg$AridTotal[match(island, snailveg$island)]))

# plot diversity data
# x-axis: species diversity per host plant species
# y-axis: functional diversity
snail_diversity %>%
  # normalize species diversity by number of available host plant species:
  mutate(spdiv=spdiv/plants) %>%
  # begin plotting:
  ggplot() +
  aes(x=spdiv, y=funcdiv, label=island, colour=zone) +
  geom_text(fontface="bold") +
  scale_x_continuous("species diversity per host plant species") +
  scale_y_continuous("functional diversity") +
  scale_colour_manual(values=c("#E69F00", "#0072B2")) +
  theme_bw()
