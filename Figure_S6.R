### PURPOSE OF THIS SCRIPT
## Reproduce Figure S6 of EZbakR-suite paper.
##
## Workflow:
## 1) Analyze subcellular NR-seq data with EZbakR and model including nuclear
##    degradation
## 2) Make figures assessing prevalence of nuclear degradation
##
## Required data:
## 1) Arrow dataset made from total RNA, cytoplasmic RNA, and nuclear RNA
##    cB files from Ietswaart et al 2024 paper.
## 2) Annotation used to create cB files, for length normalization of read
##    counts.
## 3) Excel sheets from Ietswaart et al. from which previously identified
##    PUNDs can be identified.
##
## Estimated total runtime: ~1 hour

# Load dependencies ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(EZbakR)
library(MASS)
library(rtracklayer)
library(GenomicFeatures)
library(readxl)
library(arrow)

# Directory to save figures and data to
savedir <- getwd()

# Path to annotation file used for processing data
gtf_path <- "Hs_ensembl.gtf"

# Path to arrow dataset containing processed Subcellular TL-seq data from
# Ietswaart et al., 2024
arrow_dataset_path <- "subtlseq_dataset/"

# Table S1 from the Subcellular TimeLapse-seq paper
# Paper link: https://www.cell.com/molecular-cell/fulltext/S1097-2765(24)00511-2#:~:text=Thus%2C%20RNA%20flow%20impacts%20cell,processing%2C%20including%20splicing%20and%20polyadenylation.
previous_pund_sheets <- "mmc2.xlsx"

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


# Figure S6 --------------------------------------------------------------------

### STEP 1: ESTIMATE FRACTIONS

# Which feature to use?
feature_to_use <- "XF"
chr_list <- list(XF ~ CH)
nuc_list <- list(XF ~ NP + CH)
cyt_list <- list(XF ~ CY + PL)
pl_list <- list(XF ~ PL)
tot_list <- list(XF ~ CH + CY + NP + PL)


ds <- open_dataset(arrow_dataset_path)

comps <- c("cyto", "nuc", "total")
times <- c("0min", "15min", "30min", "60min", "120min")
metadf <- tibble(sample = paste0("K562_",
                                 rep(comps, each = 10),
                                 "_WT_",
                                 rep(rep(times, each = 2), times = 3),
                                 rep(c("_rep1", "_rep2"), times = 15))) %>%
  mutate(compartment = case_when(
    grepl("_nuc_", sample) ~ "nucleus",
    grepl("_cyto_", sample) ~ "cytoplasm",
    grepl("_total_", sample) ~ "total"
  ),
  tl = case_when(
    grepl("_0min_", sample) ~ 0,
    grepl("_15min_", sample) ~ 15,
    grepl("_30min_", sample) ~ 30,
    grepl("_60min_", sample) ~ 60,
    grepl("_120min_", sample) ~ 120
  ))


ezbado <- EZbakRArrowData(ds, metadf)

### STEP 3: ESTIMATE SCALE FACTORS

ezbado <- EstimateFractions(ezbado,
                            features = feature_to_use,
                            pold_from_nolabel = TRUE)

ezbado_c <- CorrectDropout(ezbado)


### STEP 4: DEFINE MODEL


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
total_list <- list(XF ~ C + N)
nuc_list <- list(XF ~ N)
cyt_list <- list(XF ~ C)


### STEP 5: GET FEATURE LENGTHS

gtf <- rtracklayer::import(gtf_path)
txdb <- makeTxDbFromGRanges(gtf)

exons.list.per.gene <- exonsBy(txdb,by="gene")

exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))

exon_df <- tibble(exon_width = exonic.gene.sizes,
                  gene_id = names(exonic.gene.sizes))


lengths <- exon_df %>%
  dplyr::rename(XF = gene_id,
                length = exon_width)



### STEP 6: AVERAGE


# With separate compartments as separate samples
ezbado_c <- AverageAndRegularize(ezbado_c,
                                 formula_mean = ~tl:compartment - 1,
                                 type = "fractions",
                                 feature_lengths = lengths,
                                 parameter = "logit_fraction_highTC")



### STEP 7: Fit model

ezbdo <- EZDynamics(ezbado_c,
                    graph = graph,
                    sub_features = "XF",
                    grouping_features = "XF",
                    sample_feature = "compartment",
                    modeled_to_measured = list(
                      total = total_list,
                      nucleus = nuc_list,
                      cytoplasm = cyt_list
                    ))

setwd(savedir)
saveRDS(ezbdo, "CytoNuc_XF_Nucdeg_RPK_normalized.rds")


### STEP 8: How common is nuclear deg, and how do the putative PUNDs look?


churchman_ests <- readxl::read_excel(previous_pund_sheets,
                                     sheet = 2)


ez_ests <- ezbdo$dynamics$dynamics1

ez_ests <- ez_ests %>%
  dplyr::rename(Gene = XF) %>%
  inner_join(churchman_ests %>%
               dplyr::select(Gene, Symbol, PUND),
             by = "Gene")


conf_cutoff <- 0.3

ez_ests <- ez_ests %>%
  dplyr::mutate(kdn_confidence = factor(case_when(
    se_logk2 < conf_cutoff ~ "high",
    se_logk2 >= conf_cutoff ~ "low",
  ), levels = c("low", "high")))

gkdn <- ez_ests %>%
  dplyr::mutate(PUND = factor(PUND, levels = c("TRUE", "FALSE"))) %>%
  dplyr::arrange(desc(PUND)) %>%
  na.omit() %>%
  ggplot(aes(x = logk2, y = logk3,
             color = PUND,
             alpha = kdn_confidence)) +
  geom_point(size = 0.5) +
  theme_classic() +
  scale_color_manual(
    values = c("darkred", "gray50")
  ) +
  scale_alpha_manual(
    values = c(0.1, 1)
  ) +
  xlab("log(kdn)") +
  ylab("log(kp)") +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'black',
              linetype = 'dotted',
              linewidth = 0.7) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=10), #change font size of axis text
    axis.title=element_text(size=12), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=10), #change font size of legend text
    legend.title=element_text(size=12))

# How many PUNDS are in low confience set
pund_perc_lc <- sum(ez_ests$PUND[ez_ests$kdn_confidence == "low"] == "TRUE")/sum(ez_ests$kdn_confidence == "low")
# 2% of these are PUNDs
notpund_perc_lc <- 1 - pund_perc_lc

pund_perc_hc <- sum(ez_ests$PUND[ez_ests$kdn_confidence == "high"] == "TRUE")/sum(ez_ests$kdn_confidence == "high")
# 43.2% of these are PUNDs
notpund_perc_hc <- 1 - pund_perc_hc

bar_df <- tibble(
  kdn_confidence = factor(rep(c("low", "high"), each = 2),
                          levels = c("low", "high")),
  perc = c(pund_perc_lc, notpund_perc_lc,
           pund_perc_hc, notpund_perc_hc),
  PUND = factor(c("TRUE", "FALSE",
                  "TRUE", "FALSE"),
                levels = c("TRUE", "FALSE"))
)

gpp <- bar_df %>%
  ggplot(aes(x = kdn_confidence,
             y = perc,
             fill = PUND)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("kdn confidence") +
  ylab("Percentage (%)") +
  scale_fill_manual(values = c("darkred", "gray50")) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")



## Save figures
setwd(savedir)
ggsave(filename = "kdn_vs_kp_scatter_RPKnorm.pdf",
       plot = gkdn,
       width = 4.5,
       height = 3)
ggsave(filename = "pund_percentage_barplot_RPKnorm.pdf",
       plot = gpp,
       width = 2,
       height = 1.67)

