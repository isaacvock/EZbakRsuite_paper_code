### PURPOSE OF THIS SCRIPT
## Reproduce Figure S1 of EZbakR-suite paper.
##
## Workflow:
## 1) Simulate NR-seq data with 3 different plabeled's
## 2) Analyze simulated data with EZbakR, and two naive cutoff-based strategies
## 3) Make plots visualzing fraction labeled estimate acurracies
##
## Required data: None
##
## Estimated total runtime: ~1 minute


# Load dependencies ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(EZbakR)
library(MASS)
library(data.table)


# Path to save files and figures to
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

# How bad is cutoff strategy? --------------------------------------------------
#' Going to test the following combos:
#' 1) plabeled = 0.01; cutoff = 0
#' 2) plabeled = 0.025; cutoff = 0
#' 3) plabeled = 0.05; cutoff = 0
#' 4) plabeled = 0.01; cutoff = 1
#' 5) plabeled = 0.025; cutoff = 1
#' 6) plabeled = 0.05; cutoff = 1


### Simulate the three datasets (three plabeled's)

# Number of features to simulate
NF <- 3000


simdata1 <- SimulateOneRep(NF,
                           pnew = 0.01)
simdata2 <- SimulateOneRep(NF,
                           pnew = 0.05)
simdata3 <- SimulateOneRep(NF,
                           pnew = 0.025)

### plabeled = 0.01; cutoff = 0

cB1 <- simdata1$cB
gt1 <- simdata1$ground_truth


naive_ests1 <- cB1[,.(fraction_new = sum(n[TC > 0])/sum(n),
                      reads = sum(n)),
                   by = .(feature)]

compare1 <- naive_ests1 %>%
  inner_join(gt1, by = "feature")


g0.01_0 <-  compare1 %>%
  mutate(density = get_density(y = fraction_new,
                               x = true_fraction_highTC,
                               n = 200)) %>%
  ggplot(aes(x = true_fraction_highTC,
             y = fraction_new,
             color = density)) +
  geom_point(size = 0.3) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'darkred',
              linetype = "dotted",
              linewidth = 0.5) +
  scale_color_viridis_c() +
  theme_classic() +
  xlab("True fn") +
  ylab("Estimated fn") +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")



### plabeled = 0.05; cutoff = 0

cB2 <- simdata2$cB
gt2 <- simdata2$ground_truth


naive_ests2 <- cB2[,.(fraction_new = sum(n[TC > 0])/sum(n),
                      reads = sum(n)),
                   by = .(feature)]

compare0.05_0 <- naive_ests2 %>%
  inner_join(gt2, by = "feature")


g0.05_0 <-  compare0.05_0 %>%
  mutate(density = get_density(y = fraction_new,
                               x = true_fraction_highTC,
                               n = 200)) %>%
  ggplot(aes(x = true_fraction_highTC,
             y = fraction_new,
             color = density)) +
  geom_point(size = 0.3) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'darkred',
              linetype = "dotted",
              linewidth = 0.5) +
  scale_color_viridis_c() +
  theme_classic() +
  xlab("True fn") +
  ylab("Estimated fn") +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


g0.05_0


### plabeled = 0.01; cutoff = 1

cB1 <- simdata1$cB
gt1 <- simdata1$ground_truth


naive_ests1 <- cB1[,.(fraction_new = sum(n[TC > 1])/sum(n),
                      reads = sum(n)),
                   by = .(feature)]

compare0.01_1 <- naive_ests1 %>%
  inner_join(gt1, by = "feature")


g0.01_1 <-  compare0.01_1 %>%
  mutate(density = get_density(y = fraction_new,
                               x = true_fraction_highTC,
                               n = 200)) %>%
  ggplot(aes(x = true_fraction_highTC,
             y = fraction_new,
             color = density)) +
  geom_point(size = 0.3) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'darkred',
              linetype = "dotted",
              linewidth = 0.5) +
  scale_color_viridis_c() +
  theme_classic() +
  xlab("True fn") +
  ylab("Estimated fn") +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


g0.01_1

### plabeled = 0.05; cutoff = 1

cB2 <- simdata2$cB
gt2 <- simdata2$ground_truth


naive_ests2 <- cB2[,.(fraction_new = sum(n[TC > 1])/sum(n),
                      reads = sum(n)),
                   by = .(feature)]

compare0.05_1 <- naive_ests2 %>%
  inner_join(gt2, by = "feature")


g0.05_1 <-  compare0.05_1 %>%
  mutate(density = get_density(y = fraction_new,
                               x = true_fraction_highTC,
                               n = 200)) %>%
  ggplot(aes(x = true_fraction_highTC,
             y = fraction_new,
             color = density)) +
  geom_point(size = 0.3) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'darkred',
              linetype = "dotted",
              linewidth = 0.5) +
  scale_color_viridis_c() +
  theme_classic() +
  xlab("True fn") +
  ylab("Estimated fn") +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


g0.05_1



### plabeled = 0.01; mixture modeling

cB1 <- simdata1$cB
gt1 <- simdata1$ground_truth
ezbdo1 <- EZbakRData(cB1,
                     metadf = tibble(sample = 'sampleA',
                                     tl = 2))
ezbdo1 <- EstimateFractions(ezbdo1)

mmests1 <- ezbdo1$fractions$feature

compare_mm1 <- mmests1 %>%
  inner_join(gt1, by = "feature")


g0.01_mm <-  compare_mm1 %>%
  mutate(density = get_density(y = fraction_highTC,
                               x = true_fraction_highTC,
                               n = 200)) %>%
  ggplot(aes(x = true_fraction_highTC,
             y = fraction_highTC,
             color = density)) +
  geom_point(size = 0.3) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'darkred',
              linetype = "dotted",
              linewidth = 0.5) +
  scale_color_viridis_c() +
  theme_classic() +
  xlab("True fn") +
  ylab("Estimated fn") +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


g0.01_mm

### plabeled = 0.05; cutoff = 0

cB2 <- simdata2$cB
gt2 <- simdata2$ground_truth
ezbdo2 <- EZbakRData(cB2,
                     metadf = tibble(sample = 'sampleA',
                                     tl = 2))
ezbdo2 <- EstimateFractions(ezbdo2)

mmests2 <- ezbdo2$fractions$feature

compare_mm2 <- mmests2 %>%
  inner_join(gt2, by = "feature")


g0.05_mm <-  compare_mm2 %>%
  mutate(density = get_density(y = fraction_highTC,
                               x = true_fraction_highTC,
                               n = 200)) %>%
  ggplot(aes(x = true_fraction_highTC,
             y = fraction_highTC,
             color = density)) +
  geom_point(size = 0.3) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'darkred',
              linetype = "dotted",
              linewidth = 0.5) +
  scale_color_viridis_c() +
  theme_classic() +
  xlab("True fn") +
  ylab("Estimated fn") +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


g0.05_mm



### plabeled = 0.025; cutoff = 0

cB4 <- simdata3$cB
gt4 <- simdata3$ground_truth


naive_ests4 <- cB4[,.(fraction_new = sum(n[TC > 0])/sum(n),
                      reads = sum(n)),
                   by = .(feature)]

compare0.025_0 <- naive_ests4 %>%
  inner_join(gt4, by = "feature")


g0.025_0 <-  compare0.025_0 %>%
  mutate(density = get_density(y = fraction_new,
                               x = true_fraction_highTC,
                               n = 200)) %>%
  ggplot(aes(x = true_fraction_highTC,
             y = fraction_new,
             color = density)) +
  geom_point(size = 0.3) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'darkred',
              linetype = "dotted",
              linewidth = 0.5) +
  scale_color_viridis_c() +
  theme_classic() +
  xlab("True fn") +
  ylab("Estimated fn") +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


g0.025_0



### plabeled = 0.025; cutoff = 1

cB4 <- simdata3$cB
gt4 <- simdata3$ground_truth


naive_ests4 <- cB4[,.(fraction_new = sum(n[TC > 1])/sum(n),
                      reads = sum(n)),
                   by = .(feature)]

compare0.025_1 <- naive_ests4 %>%
  inner_join(gt4, by = "feature")


g0.025_1 <-  compare0.025_1 %>%
  mutate(density = get_density(y = fraction_new,
                               x = true_fraction_highTC,
                               n = 200)) %>%
  ggplot(aes(x = true_fraction_highTC,
             y = fraction_new,
             color = density)) +
  geom_point(size = 0.3) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'darkred',
              linetype = "dotted",
              linewidth = 0.5) +
  scale_color_viridis_c() +
  theme_classic() +
  xlab("True fn") +
  ylab("Estimated fn") +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


g0.025_1


### plabeled = 0.025; mixture modeling

cB4 <- simdata3$cB
gt4 <- simdata3$ground_truth
ezbdo4 <- EZbakRData(cB4,
                     metadf = tibble(sample = 'sampleA',
                                     tl = 2))
ezbdo4 <- EstimateFractions(ezbdo4)

mmests4 <- ezbdo4$fractions$feature

compare_mm4 <- mmests4 %>%
  inner_join(gt4, by = "feature")


g0.025_mm <-  compare_mm4 %>%
  mutate(density = get_density(y = fraction_highTC,
                               x = true_fraction_highTC,
                               n = 200)) %>%
  ggplot(aes(x = true_fraction_highTC,
             y = fraction_highTC,
             color = density)) +
  geom_point(size = 0.3) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'darkred',
              linetype = "dotted",
              linewidth = 0.5) +
  scale_color_viridis_c() +
  theme_classic() +
  xlab("True fn") +
  ylab("Estimated fn") +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


g0.025_mm




### Save figures

setwd(savedir)
ggsave(filename = "plabeled_0.01_cutoff_0.pdf",
       plot = g0.01_0,
       width = 2,
       height = 1.67)
ggsave(filename = "plabeled_0.025_cutoff_0.pdf",
       plot = g0.025_0,
       width = 2,
       height = 1.67)
ggsave(filename = "plabeled_0.05_cutoff_0.pdf",
       plot = g0.05_0,
       width = 2,
       height = 1.67)


ggsave(filename = "plabeled_0.01_cutoff_1.pdf",
       plot = g0.01_1,
       width = 2,
       height = 1.67)
ggsave(filename = "plabeled_0.025_cutoff_1.pdf",
       plot = g0.025_1,
       width = 2,
       height = 1.67)
ggsave(filename = "plabeled_0.05_cutoff_1.pdf",
       plot = g0.05_1,
       width = 2,
       height = 1.67)



ggsave(filename = "plabeled_0.01_mm.pdf",
       plot = g0.01_mm,
       width = 2,
       height = 1.67)
ggsave(filename = "plabeled_0.025_mm.pdf",
       plot = g0.025_mm,
       width = 2,
       height = 1.67)
ggsave(filename = "plabeled_0.05_mm.pdf",
       plot = g0.05_mm,
       width = 2,
       height = 1.67)




