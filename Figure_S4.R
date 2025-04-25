### PURPOSE OF THIS SCRIPT
## Reproduce Figure S4 of EZbakR-suite paper.
##
## Workflow:
## 1) Simulate NR-seq data with simple nuclear -> cytoplasmic model.
## 2) Estimate parameters of model with 4 strategies: 1) no read count
##    modeling; 2) library size normalization; 3) EZbakR normalization;
##    4) Spike-in normalization (i.e., provide true scale factors to EZbakR)
## 3) Make plots assessing accuracy of parameter estimates
##
## Required data: None
##
## Estimated total runtime: ~1 hour



# Load dependencies ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(EZbakR)
library(MASS)


# Directory to save figures and data to
savedir <- getwd()


# Source: https://slowkow.com/notes/ggplot2-color-by-density/
# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}



# Figure S3 --------------------------------------------------------------------
# Simulate data from 0 -> N -> C -> 0 model and fit with either fraction new
# only, library size normalized, or EZbakR normalized data

##### SIMULATE DATA

# Number of features to simulate
nfeatures <- 1000

# graph
graph <- matrix(c(0, 1, 0,
                  0, 0, 2,
                  3, 0, 0),
                nrow = 3,
                ncol = 3,
                byrow = TRUE)

colnames(graph) <- c("0", "N", "C")
rownames(graph) <- colnames(graph)

# formula list
total_list <- list(GF ~ C + N)
nuc_list <- list(GF ~ N)
cyt_list <- list(GF ~ C)

formula_list <- list(sampleA = total_list,
                     sampleB = total_list,
                     sampleC = total_list,
                     sampleD = total_list,
                     sampleE = nuc_list,
                     sampleF = nuc_list,
                     sampleG = nuc_list,
                     sampleH = nuc_list,
                     sampleI = cyt_list,
                     sampleJ = cyt_list,
                     sampleK = cyt_list,
                     sampleL = cyt_list)

# metadf
metadf <- dplyr::tibble(sample = paste0("sample", LETTERS[1:length(formula_list)]),
                        compartment = rep(c("total", "nuclear", "cytoplasm"),
                                          each = 4),
                        tl = rep(c(3, 3,
                                   1, 1), times = 3))

# means of log of parameters
log_means <- c(1, -0.3, -2)

# population sds on log scale of parameters
log_sds <- rep(0.4, times = max(graph))

# Unassigned indicator
unassigned_name <- "__no_feature"

# Sequencing depth
seqdepth <- nfeatures * 2500

# Negative binomial dispersion parameter (size in `rnbinom()`)
dispersion <- 1000

# Logit(fn) replicate variability (homoskedastic for now)
lfn_sd <- 0.2



simdata <- SimulateDynamics(nfeatures=nfeatures,
                            graph = graph,
                            metadf = metadf,
                            formula_list = formula_list,
                            log_means = log_means,
                            log_sds = log_sds,
                            unassigned_name = unassigned_name,
                            seqdepth = seqdepth,
                            dispersion = dispersion,
                            lfn_sd = lfn_sd)


metadf <- metadf

ezbdo <- EZbakRData(simdata$cB, metadf)

ezbdo <- EstimateFractions(ezbdo)


ezbdo <- AverageAndRegularize(ezbdo,
                              formula_mean = ~tl:compartment - 1,
                              type = "fractions",
                              parameter = "logit_fraction_highTC")


### Estimate parameters

ezbdo_ez <- EZDynamics(ezbdo,
                       graph = graph,
                       sub_features = "GF",
                       grouping_features = "GF",
                       sample_feature = "compartment",
                       modeled_to_measured = list(
                         total = total_list,
                         nuclear = nuc_list,
                         cytoplasm = cyt_list
                       ))


### Assess accuracy
gt <- simdata$ground_truth$parameter_truth

# Change names to reflect old naming convention used to make
# the plots originally
dynfit <- ezbdo_ez$dynamics$dynamics1 %>%
  dplyr::rename(
    k1 = logk1,
    k2 = logk2,
    k3 = logk3,
    k1_se = se_logk1,
    k2_se = se_logk2,
    k3_se = se_logk3
  ) %>%
  dplyr::filter(!((k1 > 9.9 | k1 < -9.9) |
                    (k2 > 9.9 | k2 < -9.9) |
                    (k3 > 9.9 | k3 < -9.9)))

compare <- dplyr::inner_join(dynfit, gt %>% dplyr::rename(GF = feature),
                             by = "GF")


true_scale_factor <- mean(exp(compare$k1[compare$k1 < 9.9]) / compare$true_k1[compare$k1 < 9.9])


gPk1 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k1),
    y = log(exp(k1)/true_scale_factor),
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k1),
             y = log(exp(k1)/true_scale_factor),
             color = density)) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true ksyn)") +
  ylab("log(estimated ksyn)") +
  geom_abline(slope =1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.5,
              linetype = 'dotted') +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10)) +
  theme(legend.position = "none")




gPk2 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k2),
    y = k2,
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k2),
             y = k2,
             color = density)) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true kexp)") +
  ylab("log(estimated kexp)") +
  geom_abline(slope =1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.5,
              linetype = 'dotted') +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


gPk3 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k3),
    y = k3,
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k3),
             y = k3,
             color = density)) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true kdeg)") +
  ylab("log(estimated kdeg)") +
  geom_abline(slope =1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.5,
              linetype = 'dotted') +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")

gPk1
gPk2
gPk3

### Save figures
setwd(savedir)
ggsave(
  filename = "EZnormalization_ksyn_accuracy.pdf",
  plot = gPk1,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "EZnormalization_kexp_accuracy.pdf",
  plot = gPk2,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "EZnormalization_kdeg_accuracy.pdf",
  plot = gPk3,
  width = 2,
  height = 1.67
)


#### FRACTION NEW ONLY ESTIMATES


ezbdo_un <- EZDynamics(ezbdo,
                       graph = graph,
                       sub_features = "GF",
                       grouping_features = "GF",
                       sample_feature = "compartment",
                       modeled_to_measured = list(
                         total = total_list,
                         nuclear = nuc_list,
                         cytoplasm = cyt_list
                       ),
                       use_coverage = FALSE)


### Assess accuracy
gt <- simdata$ground_truth$parameter_truth

# Currently, EZbakR removes synthesis parameter
# if normalization impossible; future versions
# will change this behavior though to perform
# post-hoc normalization.
dynfit <- ezbdo_un$dynamics$dynamics1 %>%
  dplyr::rename(
    k3 = logk2,
    k2 = logk1,
    se_k3 = se_logk2,
    se_k2 = se_logk1,
  ) %>%
  dplyr::filter(!((k2 > 9.9 | k2 < -9.9) |
                    (k3 > 9.9 | k3 < -9.9)))

compare <- dplyr::inner_join(dynfit, gt %>% dplyr::rename(GF = feature),
                             by = "GF")



gPk2_un <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k2),
    y = k2,
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k2),
             y = k2,
             color = density)) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true kexp)") +
  ylab("log(estimated kexp)") +
  geom_abline(slope =1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.5,
              linetype = 'dotted') +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


gPk3_un <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k3),
    y = k3,
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k3),
             y = k3,
             color = density)) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true kdeg)") +
  ylab("log(estimated kdeg)") +
  geom_abline(slope =1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.5,
              linetype = 'dotted') +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")

setwd(savedir)

ggsave(
  filename = "Unnormalized_kexp_accuracy.pdf",
  plot = gPk2_un,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "Unnormalized_kdeg_accuracy.pdf",
  plot = gPk3_un,
  width = 2,
  height = 1.67
)

#### LIBRARY SIZE NORMALIZED READ COUNTS


scale_df_unscaled <- tibble(
  compartment = unique(metadf$compartment)
) %>%
  mutate(scale = 1)

ezbdo_ls <- EZDynamics(ezbdo,
                       graph = graph,
                       sub_features = "GF",
                       grouping_features = "GF",
                       sample_feature = "compartment",
                       scale_factors = scale_df_unscaled,
                       modeled_to_measured = list(
                         total = total_list,
                         nuclear = nuc_list,
                         cytoplasm = cyt_list
                       ))


### Assess accuracy
gt <- simdata$ground_truth$parameter_truth

dynfit <- ezbdo_ls$dynamics$dynamics1 %>%
  dplyr::rename(
    k1 = logk1,
    k2 = logk2,
    k3 = logk3,
    k1_se = se_logk1,
    k2_se = se_logk2,
    k3_se = se_logk3
  ) %>%
  dplyr::filter(!((k1 > 9.9 | k1 < -9.9) |
                    (k2 > 9.9 | k2 < -9.9) |
                    (k3 > 9.9 | k3 < -9.9)))


compare <- dplyr::inner_join(dynfit, gt %>% dplyr::rename(GF = feature),
                             by = "GF")


gPk1_ls <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k1),
    y = log(exp(k1)/true_scale_factor),
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k1),
             y = log(exp(k1)/true_scale_factor),
             color = density)) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true ksyn)") +
  ylab("log(estimated ksyn)") +
  geom_abline(slope =1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.5,
              linetype = 'dotted') +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")




gPk2_ls <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k2),
    y = k2,
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k2),
             y = k2,
             color = density)) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true kexp)") +
  ylab("log(estimated kexp)") +
  geom_abline(slope =1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.5,
              linetype = 'dotted') +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


gPk3_ls <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k3),
    y = k3,
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k3),
             y = k3,
             color = density)) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true kdeg)") +
  ylab("log(estimated kdeg)") +
  geom_abline(slope =1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.5,
              linetype = 'dotted') +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")

setwd(savedir)
ggsave(
  filename = "LibrarySizeNorm_ksyn_accuracy.pdf",
  plot = gPk1_ls,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "LibrarySizeNorm_kexp_accuracy.pdf",
  plot = gPk2_ls,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "LibrarySizeNorm_kdeg_accuracy.pdf",
  plot = gPk3_ls,
  width = 2,
  height = 1.67
)


setwd(savedir)
saveRDS(ezbdo_ez,
        "NucCyto_EZnormalized_fit.rds")
saveRDS(ezbdo_un,
        "NucCyto_Unnormalized_fit.rds")
saveRDS(ezbdo_ls,
        "NucCyto_LibrarySizeNorm_fit.rds")



#### Spike-ins ####
# For simplicity, will just provide true scale factors to EZbakR
# to simulate best case spike-in strategy.


metadf <- metadf

ezbdo <- EZbakRData(simdata$cB, metadf)

ezbdo <- EstimateFractions(ezbdo)


ezbdo <- AverageAndRegularize(ezbdo,
                              formula_mean = ~tl:compartment - 1,
                              type = "fractions",
                              parameter = "logit_fraction_highTC")

### Estimate scale factors

get_solutions <- function(parameter_ests, graph, tl){

  # Parameters are on log-scale for ease of optimization
  param_extend <- c(0, exp(parameter_ests))
  param_graph <- matrix(param_extend[graph + 1],
                        nrow = nrow(graph),
                        ncol = ncol(graph),
                        byrow = FALSE)

  # A is same size as graph minus the "0" row and column
  A <- matrix(0,
              nrow = nrow(graph) - 1,
              ncol = ncol(graph) - 1)

  rownames(A) <- rownames(graph[-1,])
  colnames(A) <- rownames(graph[-1,])



  zero_index <- which(colnames(graph) == "0")

  # Left as an exercise to the reader
  # Just kidding, I'll document this somewhere.
  # TO-DO: Should really compile documentation for all
  # non-obvious mathematical results
  diag(A) <- -rowSums(param_graph[-zero_index,])
  A <- A + t(param_graph[-zero_index,-zero_index])


  ### Step 2: infer general solution

  Rss <- solve(a = A,
               b = -param_graph[zero_index,-zero_index])


  ev <- eigen(A)

  lambda <- ev$values
  V<- ev$vectors
  cs <- solve(V, -Rss)


  ### Step 3: Infer data for actual measured species
  exp_lambda <- exp(lambda*tl)

  scaled_eigenvectors <- V %*% diag(exp_lambda*cs)

  result_vector <- rowSums(scaled_eigenvectors) + Rss

  names(result_vector) <- rownames(A)

  fns <- result_vector/Rss
  return(list(fns = fns,
              ss = Rss))
}


true_params <- simdata$ground_truth$parameter_truth

Ns <- rep(0, times = nrow(true_params))
Cs <- Ns

for(i in 1:nrow(true_params)){

  ss <- get_solutions(log(c(simdata$ground_truth$parameter_truth$true_k1[i],
                            simdata$ground_truth$parameter_truth$true_k2[i],
                            simdata$ground_truth$parameter_truth$true_k3[i])),
                      graph,
                      1)$ss

  Ns[i] <- ss[1]
  Cs[i] <- ss[2]

}

scale <- sum(Ns) / (sum(Cs) + sum(Ns))

scale_df <- tibble(scale = c(1, (1-scale)/scale, 1 + ((1-scale)/scale)),
                   compartment = c("nuclear", "cytoplasm", "total"))
scale_df

### Estimate parameters

ezbdo_ez <- EZDynamics(ezbdo,
                       graph = graph,
                       sub_features = "GF",
                       grouping_features = "GF",
                       sample_feature = "compartment",
                       scale_factors = scale_df,
                       modeled_to_measured = list(
                         total = total_list,
                         nuclear = nuc_list,
                         cytoplasm = cyt_list
                       ))


### Assess accuracy
gt <- simdata$ground_truth$parameter_truth

dynfit <- ezbdo_ez$dynamics$dynamics1 %>%
  dplyr::filter(!((logk1 > 9.9 | logk1 < -9.9) |
                    (logk2 > 9.9 | logk2 < -9.9) |
                    (logk3 > 9.9 | logk3 < -9.9)))

compare <- dplyr::inner_join(dynfit, gt %>% dplyr::rename(GF = feature),
                             by = "GF")

true_scale_factor <- mean(exp(compare$logk1[compare$logk1 < 9.9]) / compare$true_k1[compare$logk1 < 9.9])


gPk1 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k1),
    y = log(exp(logk1)/true_scale_factor),
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k1),
             y = log(exp(logk1)/true_scale_factor),
             color = density)) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true ksyn)") +
  ylab("log(estimated ksyn)") +
  geom_abline(slope =1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.5,
              linetype = 'dotted') +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10)) +
  theme(legend.position = "none")




gPk2 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k2),
    y = logk2,
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k2),
             y = logk2,
             color = density)) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true kexp)") +
  ylab("log(estimated kexp)") +
  geom_abline(slope =1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.5,
              linetype = 'dotted') +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


gPk3 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k3),
    y = logk3,
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k3),
             y = logk3,
             color = density)) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true kdeg)") +
  ylab("log(estimated kdeg)") +
  geom_abline(slope =1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.5,
              linetype = 'dotted') +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")



### Save figures
setwd(savedir)
ggsave(
  filename = "SpikeInNormalization_ksyn_accuracy.pdf",
  plot = gPk1,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "SpikeInNormalization_kexp_accuracy.pdf",
  plot = gPk2,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "SpikeInNormalization_kdeg_accuracy.pdf",
  plot = gPk3,
  width = 2,
  height = 1.67
)

