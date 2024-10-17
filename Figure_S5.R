### PURPOSE OF THIS SCRIPT
## Reproduce Figure S5 of EZbakR-suite paper.
##
## Workflow:
## 1) Simulate NR-seq data with combined pre-RNA/nuclear/cytoplasmic model
## 2) Analyze data with EZbakR
## 3) Make plots assessing accuracy of parameter estimates
##
## Required data: None
##
## Estimated total runtime: ~30 minutes


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



# Figure S4 --------------------------------------------------------------------

##### SIMULATE DATA

# Number of features to simulate
nfeatures <- 1000

# graph
graph <- matrix(c(0, 1, 0, 0, 0,
                  0, 0, 2, 3, 0,
                  0, 0, 0, 0, 4,
                  0, 0, 0, 0, 5,
                  6, 0, 0, 0, 0),
                nrow = 5,
                ncol = 5,
                byrow = TRUE)

colnames(graph) <- c("0", "NP", "NM", "CP","CM")
rownames(graph) <- colnames(graph)

# formula list
total_list <- list(GF ~ NP + CP,
                   XF ~ NM + CM)
nuc_list <- list(GF ~ NP,
                 XF ~ NM)
cyt_list <- list(GF ~ CP,
                 XF ~ CM)

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
log_means <- c(1, seq(from = -0.3, to = -2, length.out = max(graph) - 1))

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

cB <- simdata$cB %>%
  dplyr::mutate(feature = dplyr::case_when(
    GF == "__no_feature" ~ XF,
    .default = GF
  ))

ezbdo <- EZbakRData(cB, metadf)

ezbdo <- EstimateFractions(ezbdo)


ezbdo <- AverageAndRegularize(ezbdo,
                              formula_mean = ~tl:compartment - 1,
                              type = "fractions",
                              parameter = "logit_fraction_highTC",
                              min_reads = 1)

### Estimate parameters

ezbdo <- EZDynamics(ezbdo,
                    graph = graph,
                    sub_features = c("GF", "XF"),
                    grouping_features = "feature",
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
dynfit <- ezbdo$dynamics$dynamics1 %>%
  dplyr::rename(
    k1 = logk1,
    k2 = logk2,
    k3 = logk3,
    k4 = logk4,
    k5 = logk5,
    k6 = logk6,
    k1_se = se_logk1,
    k2_se = se_logk2,
    k3_se = se_logk3,
    k4_se = se_logk4,
    k5_se = se_logk5,
    k6_se = se_logk6
  )


compare <- dplyr::inner_join(dynfit, gt,
                             by = "feature")


scale_factor <- mean(exp(compare$k1) / compare$true_k1)

point_size <- 0.4
gPk1 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k1),
    y = log(exp(k1)/scale_factor),
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k1),
             y = log(exp(k1)/scale_factor),
             color = density)) +
  geom_point(size=point_size) +
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
  geom_point(size=point_size) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true kNPtoNM)") +
  ylab("log(estimated kNPtoNM)") +
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
  geom_point(size=point_size) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true kNPtoCP)") +
  ylab("log(estimated kNPtoCP)") +
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
  geom_point(size=point_size) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true kNMtoCM)") +
  ylab("log(estimated kNMtoCM)") +
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


gPk5 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k5),
    y = k5,
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k5),
             y = k5,
             color = density)) +
  geom_point(size=point_size) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true kCPtoCM)") +
  ylab("log(estimated kCPtoCM)") +
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


gPk6 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k6),
    y = k6,
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k6),
             y = k6,
             color = density)) +
  geom_point(size=point_size) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(true kCMdeg)") +
  ylab("log(estimated kCMdeg)") +
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
gPk4
gPk5
gPk6

setwd(savedir)
saveRDS(ezbdo,
        "NucCyto_PandM_EZnormalized_fit.rds")


setwd(savedir)
ggsave(
  filename = "CytoNuc_PandM_ksyn_accuracy.pdf",
  plot = gPk1,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "CytoNuc_PandM_kNPtoNM_accuracy.pdf",
  plot = gPk2,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "CytoNuc_PandM_kNPtoCP_accuracy.pdf",
  plot = gPk3,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "CytoNuc_PandM_kNMtoCM_accuracy.pdf",
  plot = gPk4,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "CytoNuc_PandM_kCPtoCM_accuracy.pdf",
  plot = gPk5,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "CytoNuc_PandM_kdeg_accuracy.pdf",
  plot = gPk6,
  width = 2,
  height = 1.67
)

