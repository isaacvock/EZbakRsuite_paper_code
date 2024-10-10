### PURPOSE OF THIS SCRIPT
## Reproduce Figure 3 of EZbakR-suite paper.
##
## Workflow:
## 1) Simulate TILAC data
## 2) Analyze simulated data with EZbakR
## 3) Make plots assessing fraction estimate accuracy
##
## Required data: None
##
## Estimated total runtime: < 10 minutes


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(EZbakR)


# Path to save figures to
savedir <- getwd()


# Make Figure 3 ----------------------------------------------------------------

# Takes several minutes
simdata <- VectSimulateMultiLabel(5000,
                                  populations = c("TC", "GA"),
                                  fraction_design = tilac_fraction_design,
                                  phighs = setNames(c(0.05, 0.05), c("TC", "GA")))


metadf <- tibble(sample = "sampleA",
                 tl_TC = 1,
                 tl_GA = 1)

ezbdo <- EZbakRData(simdata$cB,
                    metadf)

# Takes several minutes
ezbdo <- EstimateFractions(ezbdo,
                           fraction_design = tilac_fraction_design)


fractions <- EZget(ezbdo, type = "fractions")

gt <- simdata$ground_truth
compare <- inner_join(fractions, gt,
                      by = 'feature')


g1 <- compare %>%
  arrange(n) %>%
  ggplot(aes(x = true_fraction_highTC_lowGA,
             y = fraction_highTC_lowGA,
             color = log10(n))) +
  geom_point(size = 0.2) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("truth") +
  ylab("estimate") +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.3,
              linetype = 'dotted') +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


g2 <- compare %>%
  arrange(n) %>%
  ggplot(aes(x = true_fraction_lowTC_highGA,
             y = fraction_lowTC_highGA,
             color = log10(n))) +
  geom_point(size = 0.2) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("truth") +
  ylab("estimate") +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.3,
              linetype = 'dotted') +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")



g3 <- compare %>%
  arrange(n) %>%
  ggplot(aes(x = true_fraction_lowTC_lowGA,
             y = fraction_lowTC_lowGA,
             color = log10(n))) +
  geom_point(size = 0.2) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("truth") +
  ylab("estimate") +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.3,
              linetype = 'dotted') +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")

# Save figures
setwd(savedir)
ggsave(filename = "fraction_highTC_lowGA_acc.pdf",
       plot = g1,
       width = 2,
       height = 1.67)
ggsave(filename = "fraction_lowTC_highGA_acc.pdf",
       plot = g2,
       width = 2,
       height = 1.67)
ggsave(filename = "fraction_lowTC_lowGA_acc.pdf",
       plot = g3,
       width = 2,
       height = 1.67)

