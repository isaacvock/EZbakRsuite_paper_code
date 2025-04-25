### PURPOSE OF THIS SCRIPT
## Reproduce Figure S3 of EZbakR-suite paper.
##
## Workflow:
## 1) Analyze simulated data with INSPEcT
## 2) Analyze real NR-seq data with pre-RNA processing model in Figure 5A.
##    Data from Ietswaart et al., 2024.
## 3) Make plots assessing premature RNA dynamics and trends
## Required data:
## 1) cB file for total RNA data from Ietswaart et al., 2024.
## 2) Annotation file used in making cB file, for length normalizing
##    read counts.
##
## Estimated total runtime: ~1 hour


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(data.table)
library(EZbakR)
library(rtracklayer)
library(tidyr)
library(readxl)
library(dplyr)
library(INSPEcT)

# Path to save files and figures to
savedir <- getwd()

# Path to cB file
cB_path <- "cB_ensemblLvl1and2_totRNAsubtlseq.csv.gz"

# Path to GTF file
gtf_path <- "Hs_ensembl_lvl1_and_2.gtf"


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


# Figure S3A-B -----------------------------------------------------------------

##### Simulate data #####

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
                        tl = c(1, 1,
                               3, 3))

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



simdata <- SimulateDynamics(nfeatures = nfeatures,
                            graph = graph,
                            metadf = metadf,
                            formula_list = formula_list,
                            log_means = log_means,
                            log_sds = log_sds,
                            unassigned_name = unassigned_name,
                            seqdepth = seqdepth,
                            dispersion = dispersion,
                            lfn_sd = lfn_sd)


##### Run INSPEcT with mutation cutoff fraction new estimates #####


read_cts_nascent <- simdata$cB %>%
  dplyr::filter(
    sample %in% c("sampleC", "sampleD")
  ) %>%
  dplyr::group_by(sample, GF, XF) %>%
  dplyr::summarise(
    nascent_reads = sum(n[TC > 0])
  ) %>%
  dplyr::select(
    sample, GF, XF, nascent_reads
  ) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(
    names_from = "sample",
    values_from = nascent_reads,
    values_fill = 0
  )


nascent_exon <- read_cts_nascent %>%
  filter(!grepl("__", XF)) %>%
  dplyr::select(-GF)

nascent_intron <- read_cts_nascent %>%
  filter(!grepl("__", GF)) %>%
  dplyr::select(-XF)

nascent_exon_mat <- nascent_exon %>%
  dplyr::select(-XF) %>%
  as.matrix()
rownames(nascent_exon_mat) <- nascent_exon$XF

nascent_intron_mat <- nascent_intron %>%
  dplyr::select(-GF) %>%
  as.matrix()
rownames(nascent_intron_mat) <- nascent_intron$GF



read_cts_total <- simdata$cB %>%
  dplyr::filter(
    sample %in% c("sampleC", "sampleD")
  ) %>%
  dplyr::group_by(sample, GF, XF) %>%
  summarise(
    n = sum(n)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    sample, GF, XF, n
  ) %>%
  tidyr::pivot_wider(
    names_from = "sample",
    values_from = n,
    values_fill = 0
  )


total_exon <- read_cts_total %>%
  filter(!grepl("__", XF)) %>%
  dplyr::select(-GF)

total_intron <- read_cts_total %>%
  filter(!grepl("__", GF)) %>%
  dplyr::select(-XF)


total_exon_mat <- total_exon %>%
  dplyr::select(-XF) %>%
  as.matrix()
rownames(total_exon_mat) <- total_exon$XF

total_intron_mat <- total_intron %>%
  dplyr::select(-GF) %>%
  as.matrix()
rownames(total_intron_mat) <- total_intron$GF

exonwidths <- rep(1000, times = nrow(total_exon))
names(exonwidths) <- total_exon$XF

intronwidths <- rep(1000, times = nrow(total_intron))
names(intronwidths) <- total_intron$GF


matExp_DESeq2<-quantifyExpressionsFromTrCounts(
  allcounts=list(exonsCounts = total_exon_mat + total_intron_mat,
                 intronsCounts = total_intron_mat)
  ,exonsWidths=exonwidths
  ,intronsWidths=intronwidths
  ,experimentalDesign=c(1, 1))


nascExp_DESeq2<-quantifyExpressionsFromTrCounts(
  allcounts=list(exonsCounts = nascent_exon_mat + nascent_intron_mat,
                 intronsCounts = nascent_intron_mat)
  ,exonsWidths=exonwidths
  ,intronsWidths=intronwidths
  ,experimentalDesign=c(1, 1))

colSums(matExp_DESeq2$exonsExpressions)



NI <- newINSPEcT("A",
                 labeling_time = 3,
                 nascentExpressions = nascExp_DESeq2,
                 matureExpressions =  matExp_DESeq2,
                 degDuringPulse = TRUE,
                 simulatedData = FALSE)


ksyn_ests <- ratesFirstGuess(NI, "synthesis") %>%
  as_tibble(rownames = "feature")

kp_ests <- ratesFirstGuess(NI, "processing") %>%
  as_tibble(rownames = "feature")

kd_ests <- ratesFirstGuess(NI, "degradation") %>%
  as_tibble(rownames = "feature")



### Assess accuracy ###

gt <- simdata$ground_truth$parameter_truth

dynfit <- ksyn_ests %>%
  inner_join(
    kp_ests,
    by = "feature"
  ) %>%
  inner_join(
    kd_ests,
    by = "feature"
  )

compare <- dplyr::inner_join(dynfit, gt,
                             by = "feature")


scale_factor <- mean(compare$synthesis_A / compare$true_k1)


gPk1 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k1),
    y = log(synthesis_A/scale_factor),
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k1),
             y = log(synthesis_A/scale_factor),
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

gPk1


gPk2 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k2),
    y = log(processing_A),
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k2),
             y = log(processing_A),
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

gPk2

gPk3 <- compare %>%
  filter(degradation_A > 0) %>%
  dplyr::mutate(density = get_density(
    x = log(true_k3),
    y = log(degradation_A),
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k3),
             y = log(degradation_A),
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

gPk3


setwd(savedir)
ggsave(
  filename = "INSPEcT_ksyn_accuracy_mutation_cutoff.pdf",
  plot = gPk1,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "INSPEcT_kp_accuracy_mutation_cutoff.pdf",
  plot = gPk2,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "INSPEcT_kdeg_accuracy_mutation_cutoff.pdf",
  plot = gPk3,
  width = 2,
  height = 1.67
)



##### Run INSPEcT with EZbakR fraction news #####

cB <- simdata$cB %>%
  dplyr::mutate(feature = dplyr::case_when(
    GF == "__no_feature" ~ XF,
    .default = GF
  ))

metadf <- metadf

ezbdo <- EZbakRData(cB, metadf)

ezbdo <- EstimateFractions(ezbdo)


read_cts_nascent <- ezbdo$fractions$GF_XF_feature %>%
  dplyr::filter(
    sample %in% c("sampleC", "sampleD")
  ) %>%
  dplyr::mutate(
    nascent_reads = round(n*fraction_highTC)
  ) %>%
  dplyr::select(
    sample, GF, XF, nascent_reads
  ) %>%
  tidyr::pivot_wider(
    names_from = "sample",
    values_from = nascent_reads,
    values_fill = 0
  )


nascent_exon <- read_cts_nascent %>%
  filter(!grepl("__", XF)) %>%
  dplyr::select(-GF)

nascent_intron <- read_cts_nascent %>%
  filter(!grepl("__", GF)) %>%
  dplyr::select(-XF)

nascent_exon_mat <- nascent_exon %>%
  dplyr::select(-XF) %>%
  as.matrix()
rownames(nascent_exon_mat) <- nascent_exon$XF

nascent_intron_mat <- nascent_intron %>%
  dplyr::select(-GF) %>%
  as.matrix()
rownames(nascent_intron_mat) <- nascent_intron$GF



read_cts_total <- ezbdo$fractions$GF_XF_feature %>%
  dplyr::filter(
    sample %in% c("sampleC", "sampleD")
  ) %>%
  dplyr::select(
    sample, GF, XF, n
  ) %>%
  tidyr::pivot_wider(
    names_from = "sample",
    values_from = n,
    values_fill = 0
  )


total_exon <- read_cts_total %>%
  filter(!grepl("__", XF)) %>%
  dplyr::select(-GF)

total_intron <- read_cts_total %>%
  filter(!grepl("__", GF)) %>%
  dplyr::select(-XF)


total_exon_mat <- total_exon %>%
  dplyr::select(-XF) %>%
  as.matrix()
rownames(total_exon_mat) <- total_exon$XF

total_intron_mat <- total_intron %>%
  dplyr::select(-GF) %>%
  as.matrix()
rownames(total_intron_mat) <- total_intron$GF

exonwidths <- rep(1000, times = nrow(total_exon))
names(exonwidths) <- total_exon$XF

intronwidths <- rep(1000, times = nrow(total_intron))
names(intronwidths) <- total_intron$GF

matExp_DESeq2<-quantifyExpressionsFromTrCounts(
  allcounts=list(exonsCounts = total_exon_mat + total_intron_mat,
                 intronsCounts = total_intron_mat)
  ,exonsWidths=exonwidths
  ,intronsWidths=intronwidths
  ,experimentalDesign=c(1, 1))


nascExp_DESeq2<-quantifyExpressionsFromTrCounts(
  allcounts=list(exonsCounts = nascent_exon_mat + nascent_intron_mat,
                 intronsCounts = nascent_intron_mat)
  ,exonsWidths=exonwidths
  ,intronsWidths=intronwidths
  ,experimentalDesign=c(1, 1))

colSums(matExp_DESeq2$exonsExpressions)



NI <- newINSPEcT("A",
                 labeling_time = 3,
                 nascentExpressions = nascExp_DESeq2,
                 matureExpressions =  matExp_DESeq2,
                 degDuringPulse = TRUE,
                 simulatedData = FALSE)


ksyn_ests <- ratesFirstGuess(NI, "synthesis") %>%
  as_tibble(rownames = "feature")

kp_ests <- ratesFirstGuess(NI, "processing") %>%
  as_tibble(rownames = "feature")

kd_ests <- ratesFirstGuess(NI, "degradation") %>%
  as_tibble(rownames = "feature")



### Assess accuracy ###

gt <- simdata$ground_truth$parameter_truth

dynfit <- ksyn_ests %>%
  inner_join(
    kp_ests,
    by = "feature"
  ) %>%
  inner_join(
    kd_ests,
    by = "feature"
  )

compare <- dplyr::inner_join(dynfit, gt,
                             by = "feature")


scale_factor <- mean(compare$synthesis_A / compare$true_k1)


gPk1 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k1),
    y = log(synthesis_A/scale_factor),
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k1),
             y = log(synthesis_A/scale_factor),
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

gPk1


gPk2 <- compare %>%
  dplyr::mutate(density = get_density(
    x = log(true_k2),
    y = log(processing_A),
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k2),
             y = log(processing_A),
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

gPk2

gPk3 <- compare %>%
  filter(degradation_A > 0) %>%
  dplyr::mutate(density = get_density(
    x = log(true_k3),
    y = log(degradation_A),
    n = 200
  )) %>%
  ggplot(aes(x = log(true_k3),
             y = log(degradation_A),
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

gPk3


setwd(savedir)
ggsave(
  filename = "INSPEcT_ksyn_accuracy.pdf",
  plot = gPk1,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "INSPEcT_kp_accuracy.pdf",
  plot = gPk2,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "INSPEcT_kdeg_accuracy.pdf",
  plot = gPk3,
  width = 2,
  height = 1.67
)



# Figure S3C Panel 1: Proportion intronic --------------------------------------

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

setwd(savedir)
ggsave(filename = "Parameter_est_scatter.pdf",
       plot = gest_s,
       width = 2,
       height = 1.67)


