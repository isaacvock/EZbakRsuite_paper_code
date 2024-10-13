### PURPOSE OF THIS SCRIPT
## Reproduce Figure S3 of EZbakR-suite paper.
##
## Workflow:
## 1) Analyze real NR-seq data with pre-RNA processing model in Figure 5A.
##    Data from Ietswaart et al., 2024.
## 2) Make plots assessing premature RNA dynamics and trends
## Required data:
## 1) cB file for total RNA data from Ietswaart et al., 2024.
## 2) Annotation file used in making cB file, for length normalizing
##    read counts.
##
## Estimated total runtime: ~1.5-2 hours


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(data.table)
library(EZbakR)
library(rtracklayer)
library(tidyr)
library(readxl)


# Path to cB file
cB_path <- "G:/Shared drives/Matthew_Simon/IWV/EZbakR_paper/Data/Subcellular_TLseq/total_K562_lvl1and2/cB/cB.csv.gz"

# Path to GTF file
gtf_path <- "G:/Shared drives/Matthew_Simon/IWV/Annotations/Filtered_references/Hs_ensembl_lvl1_and_2.gtf"

# Path to save files and figures to
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

### Panel 1: Proportion intronic

# Load data
cB <- fread(cB_path)

gtf <- rtracklayer::import(gtf_path)


# Find intronless genes
intron_containing_genes <- dplyr::as_tibble(gtf) %>%
  dplyr::filter(type == "exon") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(intronless = ifelse(
    any(exon_number > 1),
    FALSE,
    TRUE
  )) %>%
  dplyr::filter(!intronless) %>%
  dplyr::rename(GF = gene_id)
# This is only flagging like 3 intronless genes


# Calculate intronic coverage per sample
cB[,.(intronic_reads = sum(n[XF == "__no_feature" & GF != "__no_feature"]),
      exonic_reads = sum(n[XF != "__no_feature"])),
   by = sample] %>%
  mutate(intronic_fraction = intronic_reads / (intronic_reads + exonic_reads))
# Around 38 percent

# Calculate intronic coverage per gene
intron_content <- cB[GF != "__no_feature"][sample == "K562_total_WT_0min_rep1"][,.(intronic_reads = sum(n[XF == "__no_feature"]),
                                                                                   exonic_reads = sum(n[XF != "__no_feature"])),
                                                                                by = .(sample, GF)] %>%
  mutate(intronic_fraction = intronic_reads / (intronic_reads + exonic_reads),
         reads = intronic_reads + exonic_reads) %>%
  dplyr::filter(reads > 50) %>%
  dplyr::inner_join(intron_containing_genes,
                    by = "GF")


intron_content %>%
  filter(intronic_fraction == 1) %>%
  arrange(-reads)

# Plot per-gene intronic coverage
gIC <- intron_content %>%
  ggplot(aes(x = intronic_fraction)) +
  geom_histogram(color = 'black',
                 fill = 'gray',
                 linewidth = 0.3) +
  theme_classic() +
  xlab("Fraction intronic") +
  ylab("Density") +
  coord_cartesian(xlim = c(0,1)) +
  geom_vline(
    xintercept = mean(intron_content$intronic_fraction),
    color = 'darkred',
    linetype = "dotted",
    linewidth = 0.75
  ) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")



setwd(savedir)
ggsave(filename = "Distribution_fraction_intronic.pdf",
       plot = gIC,
       width = 2,
       height = 1.67)



# Panel 2: Violin plots of kp and kdeg -----------------------------------------

if(!exists("cB")){
  cB <- fread(cB_path)
}

### Estimate fractions

metadf <- tibble(
  sample = unique(cB$sample)
) %>%
  dplyr::mutate(
    tl = case_when(
      grepl("_0min_", sample) ~ 0,
      grepl("_15min_", sample) ~ 15,
      grepl("_30min_", sample) ~ 30,
      grepl("_60min_", sample) ~ 60,
      grepl("_120min_", sample) ~ 120
    )
  ) %>% as_tibble()

ezbdo <- EZbakRData(cB,
                    metadf)

ezbdo <- EstimateFractions(ezbdo,
                           features = c("GF", "XF"),
                           filter_cols = "GF",
                           pold_from_nolabel = TRUE)


ezbdo <- CorrectDropout(ezbdo)


### Get feature lengths from gtf for normalization

# First, import the GTF-file
library(GenomicFeatures)
txdb <- makeTxDbFromGRanges(gtf)

# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")

# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))

# Then get gene lengths
genes <- genes(txdb)
gene_sizes <- width(genes)
names(gene_sizes) <- mcols(genes)$gene_id

# Combine to infer "intron" length
gene_df <- tibble(gene_width = gene_sizes,
                  gene_id = names(gene_sizes))
exon_df <- tibble(exon_width = exonic.gene.sizes,
                  gene_id = names(exonic.gene.sizes))

final_df <- gene_df %>%
  inner_join(exon_df,
             by = "gene_id") %>%
  mutate(intron_width = gene_width - exon_width)


lengths <- bind_rows(list(
  tibble(GF = final_df$gene_id,
         XF = "__no_feature",
         length = final_df$intron_width),
  tibble(GF = final_df$gene_id,
         XF = final_df$gene_id,
         length = final_df$exon_width)
))


ezbdo <- AverageAndRegularize(ezbdo,
                              formula_mean = ~ tl,
                              type = "fractions",
                              feature_lengths = lengths,
                              parameter = "logit_fraction_highTC")


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


ezbdo <- EZDynamics(ezbdo,
                    graph,
                    sub_features = c("GF", "XF"),
                    grouping_features = "GF",
                    modeled_to_measured = total_list)


ezbdo$cB <- NULL
setwd(savedir)
saveRDS(ezbdo,
        file = "DOcorrected_totalRNA_PtoM_ezfit_20240827.rds")


# Change names to reflect old naming convention used to make
# the plots originally
estimates <- ezbdo$dynamics$dynamics1 %>%
  dplyr::rename(
    k1 = logk1,
    k2 = logk2,
    k3 = logk3,
    k1_se = se_logk1,
    k2_se = se_logk2,
    k3_se = se_logk3
  )

gest_v <- estimates %>%
  filter(k2_se < 0.3 & k3_se < 0.3) %>%
  dplyr::select(GF, k2, k3) %>%
  tidyr::pivot_longer(
    cols = c(k2, k3),
    names_to = "parameter",
    values_to = "estimate"
  ) %>%
  dplyr::mutate(parameter = case_when(
    parameter == "k2" ~ "kp",
    parameter == "k3" ~ "kdeg"
  )) %>%
  dplyr::mutate(
    parameter = factor(parameter,
                       levels = c("kp", "kdeg"))
  ) %>%
  ggplot(aes(x = parameter,
             y = estimate)) +
  geom_violin(trim = F,
              color = 'black',
              fill = 'gray30',
              linewidth = 0.2) + # get rid of loud outlines
  geom_boxplot(width = 0.1,
               outlier.shape = NA, # Get rid of outlier points
               fill = 'white',
               color = NA, # Boxplot will show .25 and .75 quartiles.
               notch = T,
               linewidth = 0.2) +
  theme_classic() +
  xlab("Parameter") +
  ylab("log(estimate)") +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")

setwd(savedir)
ggsave(filename = "Parameter_est_violin_plot.pdf",
       plot = gest_v,
       width = 2,
       height = 1.67)


### Panel 3: Scatter plot comparing the two for each feature

# Change names to reflect old naming convention used to make
# the plots originally
estimates <- ezbdo$dynamics$dynamics1 %>%
  dplyr::rename(
    k1 = logk1,
    k2 = logk2,
    k3 = logk3,
    k1_se = se_logk1,
    k2_se = se_logk2,
    k3_se = se_logk3
  )

gest_s <- estimates %>%
  filter(k2_se < 0.3 & k3_se < 0.3) %>%
  dplyr::select(GF, k2, k3) %>%
  dplyr::rename(
    kp = k2,
    kdeg = k3
  ) %>%
  dplyr::mutate(
    density = get_density(
      x = kp,
      y = kdeg,
      n = 200
    )
  ) %>%
  ggplot(aes(x = kp,
             y = kdeg,
             color = density)) +
  geom_point(size = 0.3) +
  theme_classic() +
  xlab("log(kp)") +
  ylab("log(kdeg)") +
  scale_color_viridis_c() +
  geom_abline(slope = 1,
              color = 'darkred',
              linetype = 'dotted',
              linewidth = 0.75) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")

gest_s

setwd("C:/Users/isaac/Documents/Simon_Lab/EZbakR_paper/Figures/Supplemental_preRNA/")
ggsave(filename = "Parameter_est_scatter.pdf",
       plot = gest_s,
       width = 2,
       height = 1.67)


