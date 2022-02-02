source("./solve_eqs.R")
source("./diversity.R")

set.seed(24137) # set random seed (for reproducibility)

species_diversity <- function(n, q) {
  p <- n/sum(n) # convert abundances to relative abundances
  sum(p^q)^(1/(1-q))
}

functional_diversity <- function(n, m, s, q = 2, delta = 0.001) {
  limits <- range(c(m - 4 * s, m + 4 * s)) # limits of integration
  taxis <- seq(limits[1], limits[2], by = delta) # trait axis
  p <- n / sum(n) # convert abundances to relative abundances
  dens <- rep(0, length(taxis)) # vector to store diversity data per grid cell
  for (i in 1:length(n)) dens <- dens + p[i]*dnorm(x=taxis, mean=m[i], sd=s[i])
  sum(dens^q)^(1/(1-q))/length(dens) # normalized functional diversity of order q
}


solution_tab <-
  crossing(S=2:13,
           replicate=1:10,
           b=c(exp_n=function(z) 1-exp(-z-0.5),
               quad_n=function(z) 1-4*z^2),
           a=c(hier_n=function(z) (tanh(z/0.1)+1)/2,
               gamma_n=function(z) 10*dgamma(200*(z+0.05), 4, 0.5)),
           E=c("low", "medium", "high")) %>%
  filter(paste(names(b), names(a)) %in% c("exp_n hier_n", "quad_n gamma_n")) %>%
  rowwise() %>%
  mutate(E = case_when(
    E=="low" ~ list(runif(S, 0.005, 0.008)^2),
    E=="medium" ~ list(runif(S, 0.015, 0.018)^2),
    E=="high" ~ list(runif(S, 0.025, 0.028)^2)
  )) %>%
  mutate(ic=list(c(rep(0.1, S), runif(S, -0.5, 0.5), runif(S, 0.05, 0.1)^2))) %>%
  mutate(params=list(list(b=b, a=a, E=E, dz=0.005))) %>%
  ungroup() %>%
  mutate(a=names(a), b=names(b)) %>%
  mutate(sol=map2(ic, params,
                  ~ode(y=.x, times=seq(0, 1e10, l=1001), func=eqs, parms=.y) %>%
                    organize_results(.y))) %>%
  select(S, replicate, sol, b, a) %>%
  unnest(sol)


# Plot results
solution_tab %>%
  unite(type, c(b, a), sep="_") %>%
  mutate(w=if_else(str_detect(type, "_n"), "narrow", "wide")) %>%
  mutate(type=if_else(type=="quad_gamma" | type=="quad_n_gamma_n",
                      "asymmetric~competition", "hierarchical~competition")) %>%
  group_by(S, replicate, type) %>%
  mutate(E = case_when(
    pmin(E)>0.005^2 & pmax(E)<0.008^2 ~ "E[i]:~low",
    pmin(E)>0.015^2 & pmax(E)<0.018^2 ~ "E[i]:~medium",
    pmin(E)>0.025^2 & pmax(E)<0.028^2 ~ "E[i]:~high",
  )) %>%
  ungroup() %>%
  group_by(S, replicate, E, type, w) %>%
  summarise(`species diversity`=mean(S),
            spdiv=species_diversity(n, 2),
            `functional diversity`=functional_diversity(n, m, s, 2, 0.001)) %>%
  ungroup() %>%
  ggplot(aes(x=`species diversity`, y=`functional diversity`, colour=w)) +
  geom_point(alpha=0.6) +
  facet_grid(type~E, scales="free_y", labeller=label_parsed) +
  scale_colour_manual(name="competition width", values=c("#E69F00", "#0072B2")) +
  guides(colour=guide_legend(override.aes=list(alpha=1), nrow=1)) +
  theme_bw() +
  theme(panel.grid=element_blank(), legend.position="bottom")
