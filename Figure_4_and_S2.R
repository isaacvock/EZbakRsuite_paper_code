### PURPOSE OF THIS SCRIPT
## Reproduce Figures 4 and S2 of EZbakR-suite paper.
##
## Workflow:
## 1) Simulate feature-to-feature plabeled variance
## 2) Analyze simulated data with EZbakR
## 3) Make plots assessing estimate accuracy
## 4) Analyze real dataset with EZbakR (from Ietswaart et al., 2024)
## 5) Make plots assessing feature-to-feature plabeled variance
## 6) Make S2 plot comparing hierarchical and non-hierarchical estimates
##    on real data
##
## Required data:
## 1) Total RNA cB file from Ietswaart et al., 2024 data
## 2) Annotation GTF used for cB creation
##
## Estimated total runtime: ~30 minutes


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(rtracklayer)
library(MASS)
library(EZbakR)
library(data.table)
library(tidyr)

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



# Path to save figures to
savedir <- getwd()

# Path to cB file
cB_path <- "G:/Shared drives/Matthew_Simon/IWV/EZbakR_paper/Data/Subcellular_TLseq/total_K562_ensembl/cB/cB.csv.gz"

# Path to annotation GTF file
gtf_path <- "G:/Shared drives/Matthew_Simon/IWV/Genomes/Human_HISAT3N/Homo_sapiens.GRCh38.104.chr_chr.gtf"

# Analyze simulated data -------------------------------------------------------

### Simulate

simdata <- SimulateOneRep(5000,
                          feature_pnew = TRUE,
                          logit_pnew_mean = -3.25,
                          logit_pnew_sd = 0.3)


### Fit

metadf <- tibble(sample = 'sampleA',
                 tl = 2)

ezbdo <- EZbakRData(simdata$cB, metadf)

ezbdo <- EstimateFractions(ezbdo,
                           strategy = "hierarchical")


### Assess

gt <- simdata$ground_truth
fractions <- EZget(ezbdo, type = "fractions")
mutrates <- ezbdo$mutation_rates$feature_TC

compare <- fractions %>%
  inner_join(gt, by = 'feature')

compare_m <- mutrates %>%
  inner_join(gt,
             by = c("sample", "feature")) %>%
  inner_join(fractions %>%
               dplyr::select(feature, n),
             by = "feature")

gf <- compare %>%
  arrange(n) %>%
  ggplot(aes(x = true_fraction_highTC,
             y = fraction_highTC,
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
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


gm <- compare_m %>%
  arrange(n) %>%
  ggplot(aes(x = true_pnew,
             y = pnew,
             color = log10(n))) +
  geom_point(size = 0.2) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("truth") +
  ylab("estimate") +
  coord_cartesian(xlim = c(0.01, 0.1),
                  ylim = c(0.01, 0.1)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.3,
              linetype = 'dotted') +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")



gp <- ggplot(compare_m, aes(x = true_pnew)) +
  geom_histogram(color = "black",
                 fill = "gray40") +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("Simulated phighs") +
  ylab("Count") +
  coord_cartesian(xlim = c(0.01, 0.1)) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


gf
gm
gp


# Save figures
setwd(savedir)
ggsave(filename = "hier_fraction_acc.pdf",
       plot = gf,
       width = 2,
       height = 1.67)
ggsave(filename = "hier_phigh_acc.pdf",
       plot = gm,
       width = 2,
       height = 1.67)
ggsave(filename = "sim_phigh_dist.pdf",
       plot = gp,
       width = 2,
       height = 1.67)


# Analyze real data ------------------------------------------------------------

### Fit

cB <- fread(cB_path)

metadf <- tibble(sample = unique(cB$sample),
                 tl = c(0, 0, 60, 60,
                        120, 120),
                 compartment = 'total')
ezbdo <- EZbakRData(cB, metadf)

start <- Sys.time()
ezbdo_h <- EstimateFractions(ezbdo,
                             pold_from_nolabel = TRUE,
                             features = "XF",
                             strategy = "hierarchical",
                             pnew_prior_sd_max = 0.15)
end <- Sys.time()
print(end - start)
# 3.32 minutes

start <- Sys.time()
ezbdo_nh <- EstimateFractions(ezbdo,
                              pold_from_nolabel = TRUE,
                              features = "XF",
                              overwrite = FALSE)
end <- Sys.time()
print(end - start)
# 1.08 minutes

setwd(savedir)
saveRDS(ezbdo_h,
        file = "Hierarchical_MM_fit_totalRNA.rds")
saveRDS(ezbdo_nh,
        file = "Non_hierarchical_MM_fit_totalRNA.rds")


### Explore fit

mutrates <- ezbdo_h$mutation_rates$XF_TC %>%
  dplyr::filter(sample == "K562_total_WT_120min_rep1")

fractions <- ezbdo_h$fractions$XF %>%
  dplyr::filter(n >= 50)

fractions <- fractions %>%
  dplyr::select(sample, XF, fraction_highTC, n)

mutrates <- mutrates %>%
  inner_join(fractions,
             by = c('XF', 'sample'))

gps <- mutrates %>%
  filter(n > 100) %>%
  mutate(density = get_density(
    x = EZbakR:::logit(fraction_highTC),
    y = pnew,
    n = 200
  )) %>%
  ggplot(aes(x = EZbakR:::logit(fraction_highTC),
             y = pnew,
             color = density)) +
  geom_point(size = 0.25) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("logit(fraction new)") +
  ylab("K562 gene-specific phigh") +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=10), #change font size of axis text
    axis.title=element_text(size=12), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=10), #change font size of legend text
    legend.title=element_text(size=12)) + #change font size of legend title
  theme(legend.position = "none")
gps


gtf <- rtracklayer::import(
  gtf_path
)
gtf <- as_tibble(gtf)

mito_df <- gtf %>%
  dplyr::filter(type == "gene") %>%
  dplyr::select(seqnames, gene_id, gene_name) %>%
  dplyr::distinct() %>%
  dplyr::mutate(chrMT = case_when(
    seqnames == "chrMT" ~ TRUE,
    .default = FALSE
  ),
  XF = gene_id,
  MTlocalized = case_when(
    seqnames == "chrMT" ~ TRUE,
    grepl("^MT", gene_name) ~ TRUE,
    .default = FALSE
  )) %>%
  dplyr::select(XF, chrMT, MTlocalized)


gmt <- mutrates %>%
  filter(n > 100) %>%
  dplyr::inner_join(mito_df, by = "XF") %>%
  dplyr::arrange(chrMT) %>%
  dplyr::mutate(logit_fraction_highTC = EZbakR:::logit(fraction_highTC)) %>%
  ggplot(aes(x = logit_fraction_highTC,
             y = pnew,
             size = chrMT,
             color = chrMT)) +
  geom_point() +
  scale_color_manual(values = c("gray70", "darkred"))  +
  scale_size_manual(values = c(0.25, 0.75)) +
  theme_classic() +
  xlab("logit(fraction new)") +
  ylab("K562 gene-specific phigh") +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=10), #change font size of axis text
    axis.title=element_text(size=12), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=10), #change font size of legend text
    legend.title=element_text(size=12)) + #change font size of legend title
  theme(legend.position = "none")


gmt
gps

setwd(savedir)
ggsave(filename = "real_pnews_reprocessed.pdf",
       plot = gps,
       width = 3,
       height = 1.67*1.5)
ggsave(filename = "real_pnews_MTlabel_reprocessed.pdf",
       plot = gmt,
       width = 3,
       height = 1.67*1.5)




### Supplemental figures


fractions_h <- EZget(ezbdo_h, type = "fractions") %>%
  dplyr::select(sample, XF, fraction_highTC, n) %>%
  dplyr::rename(hier_fraction = fraction_highTC)

fractions_nh <- EZget(ezbdo_nh, type = "fractions") %>%
  dplyr::select(sample, XF, fraction_highTC) %>%
  dplyr::rename(nonhier_fraction = fraction_highTC)

compare_fxn <- inner_join(
  fractions_h,
  fractions_nh,
  by = c("sample", "XF")
)


ghvsnh_dens <- compare_fxn %>%
  dplyr::filter(n > 50) %>%
  dplyr::filter(sample == "K562_total_WT_120min_rep1") %>%
  dplyr::mutate(density = get_density(
    x = nonhier_fraction,
    y = hier_fraction,
    n = 200
  )) %>%
  ggplot(aes(x = nonhier_fraction,
             y = hier_fraction,
             color = density)) +
  geom_point(size = 0.25) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("non-hier fn est") +
  ylab("hier fn est") +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.7,
              linetype = "dotted") +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=10), #change font size of axis text
    axis.title=element_text(size=12), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=10), #change font size of legend text
    legend.title=element_text(size=12)) + #change font size of legend title
  theme(legend.position = "none")



ghvsnh_MT <- compare_fxn %>%
  dplyr::filter(n > 50) %>%
  dplyr::filter(sample == "K562_total_WT_120min_rep1") %>%
  dplyr::inner_join(mito_df, by = "XF") %>%
  dplyr::arrange(chrMT) %>%
  ggplot(aes(x = nonhier_fraction,
             y = hier_fraction,
             color = chrMT,
             size = chrMT)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("gray70", "darkred"))  +
  scale_size_manual(values = c(0.25, 0.75)) +
  xlab("non-hier fn est") +
  ylab("hier fn est") +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'darkred',
              linewidth = 0.7,
              linetype = "dotted") +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=10), #change font size of axis text
    axis.title=element_text(size=12), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=10), #change font size of legend text
    legend.title=element_text(size=12)) + #change font size of legend title
  theme(legend.position = "none")



setwd(savedir)
ggsave(filename = "nonhier_vs_hier_denscolor.pdf",
       plot = ghvsnh_dens,
       width = 3,
       height = 1.67*1.5)
ggsave(filename = "nonhier_vs_hier_MTcolor.pdf",
       plot = ghvsnh_MT,
       width = 3,
       height = 1.67*1.5)
