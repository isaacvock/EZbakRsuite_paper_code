### PURPOSE OF THIS SCRIPT
## Reproduce Figure S4 of EZbakR-suite paper.
##
## Workflow:
## 1) Simulate NR-seq data with simple nuclear -> cytoplasmic model.
## 2) Estimate parameters of model with 3 strategies: 1) no read count
##    modeling; 2) library size normalization; 3) EZbakR normalization
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


savedir <- getwd()

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



simdata <- SimulateDynamics(nfeatures,
                            graph,
                            metadf,
                            formula_list,
                            log_means,
                            log_sds,
                            unassigned_name,
                            seqdepth,
                            dispersion,
                            lfn_sd)


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

dynfit <- ezbdo_ez$dynamics$dynamics1 %>%
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
                       scale_factors = scale_df,
                       modeled_to_measured = list(
                         total = total_list,
                         nuclear = nuc_list,
                         cytoplasm = cyt_list
                       ),
                       use_coverage = FALSE)


### Assess accuracy
gt <- simdata$ground_truth$parameter_truth

dynfit <- ezbdo_un$dynamics$dynamics1 %>%
  dplyr::rename(
    k3 = k2,
    k2 = k1,
    se_k3 = se_k2,
    se_k2 = se_k1,
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


scale_df_unscaled <- scale_df %>%
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
             y = log(exp(k1)/scale_factor),
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

