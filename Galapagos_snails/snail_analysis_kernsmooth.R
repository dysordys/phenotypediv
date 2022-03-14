library(tidyverse) # For efficient data manipulation and plotting
library(mvtnorm) # For multivariate normal distributions
library(ks) # For multivariate kernel density estimation


# Change these parameters to adjust species- and functional diversity metrics
q_sp <- 2 # Order of species diversity index (2 corresponds to inverse Simpson)
q_fn <- 2 # Order of functional diversity index
res <- 501 # Number of bins to divide trait axis into; increase for better resolution


# Calculate functional diversity over 2D trait space
# Input:
# - pts: a data frame with two columns t1 and t2, which hold the coordinates
#        of the trait values of each individual (rows)
# - q: order of diversity index
# - grid_size: number of equally-sized square bins to discretize trait space into
# Output:
# - A single number, the community's functional diversity index
diversity_kernsmooth <- function(pts, q, grid_size) {
  s <- kde(x=pts, H=Hpi(pts), binned=FALSE, xmin=c(-3, -3), xmax=c(6, 6),
           gridsize=grid_size)$estimate
  s <- s/sum(s)
  sum(s^q)^(1/(1-q))/length(s)
}


# Load morphological data
snail_kernsmooth <- read_csv("snaildata.csv") %>%
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

snail_kernsmooth %>%
  nest(data = c(species, t1, t2)) %>%
  left_join(read_csv("vegzonetotals.csv"), by="island") %>% # Host plant sp. richnesses
  rowwise() %>% # Choose plant richness from appropriate habitat:
  mutate(plants=if_else(habitat=="arid", AridTotal, HumidTotal)) %>%
  ungroup() %>%
  select(-AridTotal, -HumidTotal) %>%
  mutate(n=map(data, ~count(.x, species)$n), # Abundance of each species
         S=map_int(n, length), # Species richness per island
         spdiv=map_dbl(n, ~sum((.x/sum(.x))^q_sp)^(1/(1-q_sp))), # Species diversity
         funcdiv=map_dbl(data, ~diversity_kernsmooth(select(.x, t1, t2), # Func. div.
                                                     q=q_fn, grid_size=res))) %>%
  mutate(spdiv=spdiv/plants) %>%
  ggplot() +
  aes(x=spdiv, y=funcdiv, label=island, colour=habitat) +
  geom_text(fontface="bold", size=4) +
  scale_x_continuous("species diversity per host plant species") +
  scale_y_continuous("functional diversity") +
  scale_colour_manual(values=c("#E69F00", "#0072B2")) +
  theme_bw() +
  theme(panel.grid=element_blank())
