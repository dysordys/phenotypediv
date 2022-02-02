require(tidyverse) # for efficient data manipulation and plotting
require(mvtnorm) # for multivariate normal distributions


# change these parameters to adjust species- and functional diversity metrics
q_sp <- 2 # order of species diversity index (2 corresponds to inverse Simpson)
q_fn <- 2 # order of functional diversity index
res <- 501 # Resolution of trait axis binning; increase for more bins


# Calculate functional diversity over 2D trait space
# Input:
# - n: vector of species abundances, with as many entries as species
# - m: list of mean trait vectors; m[[i]][k] is species i's kth vector entry
# - P: list of trait covariance matrices; P[[i]][k,l] is species i's (k,l)th entry
# - q: order of diversity index
# - taxis: trait axis, i.e. seq(from, to, l), used for both traits
# Output:
# - A single number, the community's functional diversity index
diversity <- function(n, m, P, q, taxis) {
  ztab <- unname(as.matrix(expand.grid(taxis, taxis)))
  s <- rep(0, nrow(ztab))
  for (i in 1:length(n)) s <- s+n[i]*dmvnorm(x=ztab, mean=m[[i]], sigma=P[[i]])
  s <- s/sum(s)
  sum(s^q)^(1/(1-q))/length(s)
}

# Obtain discretized trait axis (which will be used along both trait directions)
# Input:
# - snail_data: data frame for whole snail assemblage, with columns t1 and t2
# - resolution: number of bins to divide trait axis into
# Output:
# - A sequence going from three times the smallest to three times the largest
#   recorded trait value (after normalization), with the number of entries
#   equal to resolution
trait_axis <- function(snail_data, resolution) {
  lim <- range(c(snail_data$t1, snail_data$t2), na.rm=TRUE)
  seq(3*lim[1], 3*lim[2], l=resolution)
}

# Maximum likelihood estimate of mean and covariance of binormal distribution
# Input:
# - dat: tibble with 3 columns: species, t1, and t2 (for the two traits)
# Output:
# - A tibble with columns species, data, m, and P. The latter 3 are nested;
#   m and P are lists of vectors and matrices, respectively
mle_binorm <- function(dat) {
  # Sub-function to extract sample mean trait
  # Input: a tibble with 3 columns: species, t1, and t2 (but "species" is unique)
  # Output: a vector with 2 entries, the mean t1 and mean t2
  sample_mean <- function(data) c(mean(data$t1), mean(data$t2))
  # Sub-function to extract sample trait covariance
  # Input: a tibble with 3 columns: species, t1, and t2 (but "species" is unique)
  # Output: a 2x2 matrix, the estimated trait covariance
  sample_cov <- function(data) {
    tab <- as.matrix(select(data, t1, t2))
    N <- nrow(tab)
    m <- sample_mean(data)
    P <- c(0, 0)%o%c(0, 0)
    for (i in 1:N) P <- P + (tab[i,]-m)%o%(tab[i,]-m)
    unname(P/(N-1))
  }
  dat %>%
    group_by(species) %>%
    nest() %>%
    ungroup() %>%
    mutate(m=map(data, sample_mean), # Calculate mean for each species,
           P=map(data, sample_cov))  # and the covariance for each species
}


snail <- read_csv("snaildata.csv") %>%
  # Trait 1: centroid size; trait 2: PC1 of shape (explains >80% of variance):
  rename(t1="Centroidsize", t2="shape1") %>%
  # Standardize trait data by subtracting mean & dividing by std dev:
  mutate(t1=(t1-mean(t1))/sd(t1), t2=(t2-mean(t2))/sd(t2)) %>%
  # Remove three outlier satellite islands:
  filter(!(island %in% c("CH", "ED", "GA"))) %>%
  # Remove the species N. achatellinus, which has only 2 sampled individuals:
  filter(!(species=="achatellinus")) %>%
  # Replace NAs with 0s in vegetation zone columns:
  replace_na(list(vegzoneHumid=0, vegzoneTransition=0, vegzoneArid=0)) %>%
  # Convert numerical vegetation zone labels into "humid" / "arid":
  mutate(habitat=case_when(
    vegzoneHumid==1 ~ "humid",
    vegzoneArid==1 ~ "arid",
    TRUE ~ "none")
  ) %>%
  filter(habitat %in% c("humid", "arid")) %>%
  select(island, habitat, species, t1, t2) # Relevant columns only

snail %>%
  group_by(island, habitat) %>%
  nest() %>%
  ungroup() %>%
  left_join(read_csv("vegzonetotals.csv"),
            by="island") %>% # Add host plant species richness data
  rowwise() %>% # Choose plant richness from appropriate habitat:
  mutate(plants=if_else(habitat=="arid", AridTotal, HumidTotal)) %>%
  ungroup() %>%
  select(-AridTotal, -HumidTotal) %>%
  mutate(fit=map(data, mle_binorm), # Calculate mean & covariance for all species
         n=map(data, ~count(.x, species)$n), # Abundance of each species
         S=map_int(n, length), # Species richness per island
         m=map(fit, ~.x$m), # Put means in separate column
         P=map(fit, ~.x$P), # Put covariances in separate column
         spdiv=map_dbl(n, ~sum((.x/sum(.x))^q_sp)^(1/(1-q_sp))), # Sp. diversity
         funcdiv=pmap_dbl(list(n, m, P), diversity, # Func. diversity
                          q=q_fn, taxis=trait_axis(snail, res))) %>%
  mutate(spdiv=spdiv/plants) %>% # Calculate species diversity per host plant spp.
  ggplot(aes(x=spdiv, y=funcdiv, label=island, colour=habitat)) + # Plot results
  geom_text(fontface="bold", size=4) +
  scale_x_continuous("species diversity per host plant species") +
  scale_y_continuous("functional diversity") +
  scale_colour_manual(values=c("#E69F00", "#0072B2")) +
  theme_bw() +
  theme(panel.grid=element_blank())
