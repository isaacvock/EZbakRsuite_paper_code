### PURPOSE OF THIS SCRIPT
## Reproduce Figure 5 of EZbakR-suite paper.
##
## Workflow:
## 1) Simulate pre-RNA processing dynamics
## 2) Analyze simulated data with EZbakR
## 3) Make plots visualzing parameter estimate accuracy
## 4) Simulate nuclear -> cytoplasmic export with nuclear decay
## 5) Analyze simulated data with EZbakR
## 6) Make plots visualizing parameter estimat accuracy
##
## Required data: None
##
## Estimated total runtime: ~1 hour


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(EZbakR)
library(MASS)


# Path to save figures to
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


# Figure 5A --------------------------------------------------------------------


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

colnames(graph) <- c("0", "P", "M")
rownames(graph) <- colnames(graph)

# formula list
total_list <- list(GF ~ P,
                   XF ~ M)

formula_list <- list(sampleA = total_list,
                     sampleB = total_list,
                     sampleC = total_list,
                     sampleD = total_list)

# metadf
metadf <- dplyr::tibble(sample = c('sampleA', 'sampleB',
                                   'sampleC', 'sampleD'),
                        compartment = c('total', 'total',
                                        'total', 'total'),
                        tl = c(3, 3,
                               1, 1))

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


cB <- simdata$cB %>%
  dplyr::mutate(feature = dplyr::case_when(
    GF == "__no_feature" ~ XF,
    .default = GF
  ))

metadf <- metadf

ezbdo <- EZbakRData(cB, metadf)

ezbdo <- EstimateFractions(ezbdo)


ezbdo <- AverageAndRegularize(ezbdo,
                              formula_mean = ~tl - 1,
                              type = "fractions",
                              parameter = "logit_fraction_highTC")

### Estimate parameters

ezbdo <- EZDynamics(ezbdo,
                    graph = graph,
                    sub_features = c("GF", "XF"),
                    grouping_features = "feature",
                    modeled_to_measured = total_list)


### Assess accuracy
gt <- simdata$ground_truth$parameter_truth

dynfit <- ezbdo$dynamics$dynamics1

compare <- dplyr::inner_join(dynfit, gt,
                             by = "feature")


scale_factor <- mean(exp(compare$k1) / compare$true_k1)


gPk1 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k1),
    y = log(exp(k1)/scale_factor),
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k1),
             y = log(exp(k1)/scale_factor),
             color = density)) +
  geom_point(size = 0.3) +
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



gPk2 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k2),
    y = k2,
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k2),
             y = k2,
             color = density)) +
  geom_point(size = 0.3) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true kp)") +
  ylab("log(estimated kp)") +
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
  geom_point(size = 0.3) +
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
ggsave(filename = "PtoM_ksyn_accuracy.pdf",
       plot = gPk1,
       width = 2,
       height = 1.67)
ggsave(filename = "PtoM_kp_accuracy.pdf",
       plot = gPk2,
       width = 2,
       height = 1.67)
ggsave(filename = "PtoM_kdeg_accuracy.pdf",
       plot = gPk3,
       width = 2,
       height = 1.67)



# Figure 5B --------------------------------------------------------------------

##### SIMULATE DATA

# Number of features to simulate
nfeatures <- 1000

# graph
graph <- matrix(c(0, 1, 0,
                  2, 0, 3,
                  4, 0, 0),
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
log_means <- c(1, -1, -0.3, -2)

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

ezbdo <- EZDynamics(ezbdo,
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

dynfit <- ezbdo$dynamics$dynamics1 %>%
  dplyr::filter(!((k1 > 9.9 | k1 < -9.9) |
                    (k2 > 9.9 | k2 < -9.9) |
                    (k3 > 9.9 | k3 < -9.9) |
                    (k4 > 9.9 | k4 < -9.9)))

compare <- dplyr::inner_join(dynfit, gt %>% dplyr::rename(GF = feature),
                             by = "GF")


scale_factor <- mean(exp(compare$k1[compare$k1 < 9.9]) / compare$true_k1[compare$k1 < 9.9])


gPk1 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k1),
    y = log(exp(k1)/scale_factor),
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
  xlab("log(true kNdeg)") +
  ylab("log(estimated kNdeg)") +
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

gPk2se <- compare %>%
  arrange(se_k2) %>%
  ggplot(aes(x = log(true_k2),
             y = k2,
             color = se_k2)) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c(option = "plasma") +
  xlab("log(true kNdeg)") +
  ylab("log(estimated kNdeg)") +
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


gPk4 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k4),
    y = k4,
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k4),
             y = k4,
             color = density)) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true kCdeg)") +
  ylab("log(estimated kCdeg)") +
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
gPk2se
gPk3
gPk4



gPkpdeg_explore <- compare %>%
  arrange(se_k2) %>%
  ggplot(aes(x = log(true_k2),
             y = log(true_k3),
             color = se_k2)) +
  geom_point(size=0.3) +
  theme_classic() +
  scale_color_viridis_c(option = "plasma") +
  xlab("log(true kNdeg)") +
  ylab("log(true kexp)") +
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

gkpdeg_color <- compare %>%
  ggplot(aes(x = log(true_k2),
             y = log(true_k3),
             color = se_k2)) +
  geom_point(size=0.3) +
  theme_classic() +
  scale_color_viridis_c(option = "plasma") +
  xlab("log(true kpdeg)") +
  ylab("log(true kp)") +
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
    legend.title=element_text(size=10))


### Save figures
setwd(savedir)
ggsave(filename = "ksyn_accuracy.pdf",
       plot = gPk1,
       width = 2,
       height = 1.67)
ggsave(filename = "kexp_accuracy.pdf",
       plot = gPk3,
       width = 2,
       height = 1.67)
ggsave(filename = "kCdeg_accuracy.pdf",
       plot = gPk4,
       width = 2,
       height = 1.67)
ggsave(filename = "kNdeg_accuracy.pdf",
       plot = gPk2,
       width = 2,
       height = 1.67)
ggsave(filename = "kNdeg_accuracy_secolor.pdf",
       plot = gPk2se,
       width = 2,
       height = 1.67)
ggsave(filename = "kNdeg_vs_kexp_secolor.pdf",
       plot = gPkpdeg_explore,
       width = 2,
       height = 1.67)
ggsave(filename = "secolor_bar.pdf",
       plot = gkpdeg_color,
       width = 2,
       height = 1.67)
