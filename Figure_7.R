### PURPOSE OF THIS SCRIPT
## Reproduce Figure 7 (Analysis of effect of DDX3X knockdown on
## subcellular kinetics)


# Load dependencies ------------------------------------------------------------

### Packages
library(devtools)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(data.table)
library(EZbakR)
library(arrow)
library(GenomicRanges)
library(rtracklayer)


### File paths for reproducing figure

# Path toi DDX3X eCLIP peak BED file from ENCODE
DDX3X_clip_peaks <- "C:/Users/isaac/Documents/Simon_Lab/EZbakR_paper/Data/Figure_4/ENCFF901BYH_DDX3X_eCLIP.bed"

# Path to hg38 RefSeq annotation GTF
hg38_gtf_path <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/Annotation_gamut_analyses/Annotations/hg38_refseq.gtf"

# Path to EZbakR analysis output object
ezbdo_path <- "C:/Users/isaac/Documents/Simon_Lab/EZbakR_paper/Data/subtlseq_analyses/EZbakR_Comparisons_NoChromatin_DONormalized.rds"

# Path to folder to which figures are saved
figure_savedir <- "C:/Users/isaac/Documents/Simon_Lab/EZbakR_paper/Figures/Figure_6/"


### File paths for also running full analysis

# Path to arrow partitioned dataset
arrow_dataset <- "C:/Users/isaac/Yale University/Simon Lab – RNA - Documents/IWV/EZbakR_paper/Data/Subcellular_TLseq_perturbations/arrow_dataset/"



# Reproduce figure -------------------------------------------------------------


##### Annotate peaks #####

DDX3X_clip_dt <- fread(DDX3X_clip_peaks)


peaks_dt <- DDX3X_clip_dt[V7 > 3 & V8 > 10][, .(chr = V1, start = V2, end = V3, strand = V6)]


peaks_gr <- makeGRangesFromDataFrame(
  peaks_dt,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  strand.field = "strand",
  keep.extra.columns = FALSE
)

refseq <- rtracklayer::import(hg38_gtf_path)

genes_gtf <- refseq[refseq$type ==  "exon"]


overlaps <- findOverlaps(peaks_gr, genes_gtf)

# Annotate the peaks with gene details
annotated_peaks <- tibble(
  peak_chr = as.character(seqnames(peaks_gr)[queryHits(overlaps)]),
  peak_start = start(peaks_gr)[queryHits(overlaps)],
  peak_end = end(peaks_gr)[queryHits(overlaps)],
  peak_strand = as.character(strand(peaks_gr)[queryHits(overlaps)]),
  gene_id = mcols(genes_gtf)$gene_id[subjectHits(overlaps)],
  gene_name = mcols(genes_gtf)$gene_name[subjectHits(overlaps)]
)


##### Load EZbakR analysis #####


setwd("C:/Users/isaac/Documents/Simon_Lab/EZbakR_paper/Data/subtlseq_analyses/")
ezbdo <- readRDS(ezbdo_path)



##### Does DDX3X binding correlate with kexp changes? #####

DDX3X_kexp <- EZget(
  ezbdo,
  type = "comparisons",
  parameter = "logk2",
  experimental = "DDX3X"
)


ke_volc <- DDX3X_kexp %>%
  mutate(
    L2FC = difference * log2(exp(1)),
    conclusion = factor(case_when(
      padj < 0.05 & L2FC < -1 ~ "Decreased",
      padj < 0.05 & L2FC > 1 ~ "Increased",
      .default = "Not sig."
    ),
    levels = c("Decreased", "Increased", "Not sig."))
  ) %>%
  ggplot(
    aes(x = difference * log2(exp(1)),
        y = -log10(padj),
        color = conclusion)
  ) +
  geom_point(size = 0.2) +
  theme_classic() +
  scale_color_manual(
    values = c("darkcyan", "darkorange", "darkgray")
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = -1,
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 1,
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  xlab("L2FC(kexp)") +
  ylab("-log10(padj)") +
  theme(legend.position="none") +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) +
  coord_cartesian(
    xlim = c(-2.25, 2.25)
  )




DDX3X_kexp_eC <- DDX3X_kexp %>%
  left_join(
    annotated_peaks %>%
      dplyr::select(gene_id, peak_chr) %>%
      dplyr::distinct() %>%
      dplyr::rename(XF = gene_id),
    by = "XF"
  ) %>%
  dplyr::mutate(
    DDX3X_target = factor(!is.na(peak_chr),
                          levels = c(FALSE, TRUE))
  ) %>%
  arrange(desc(DDX3X_target))


ke_box <- DDX3X_kexp_eC %>%
  ggplot(
    aes(x = DDX3X_target,
        y = difference * log2(exp(1)))
  ) +
  geom_jitter(width = 0.2,
              height = 0,
              size = 0.1,
              alpha = 0.1) +
  geom_boxplot(fill = "aquamarine3",
               color = "black",
               outliers = F,
               linewidth = 0.2,
               notch = T,
               width = 0.3) +
  theme_classic() +
  geom_hline(
    yintercept = median(DDX3X_kexp_eC$difference[DDX3X_kexp_eC$DDX3X_target == FALSE]),
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  xlab("DDX3X target (ENCODE eCLIP)") +
  ylab("L2FC(kexp)") +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) +
  coord_cartesian(
    ylim = c(-2.25, 2.25)
  )


setwd(figure_savedir)
ggsave(
  filename = "kexp_volcano_plot.pdf",
  plot = ke_volc,
  width = 2.5,
  height = 2
)

ggsave(
  filename = "kexp_box_plot.pdf",
  plot = ke_box,
  width = 2.5,
  height = 2
)


##### Does DDX3X binding correlate with kdeg changes? #####

DDX3X_kdeg <- EZget(
  ezbdo,
  type = "comparisons",
  parameter = "logk3",
  experimental = "DDX3X"
)


kd_volc <- DDX3X_kdeg %>%
  mutate(
    L2FC = difference * log2(exp(1)),
    conclusion = factor(case_when(
      padj < 0.05 & L2FC < -1 ~ "Decreased",
      padj < 0.05 & L2FC > 1 ~ "Increased",
      .default = "Not sig."
    ),
    levels = c("Decreased", "Increased", "Not sig."))
  ) %>%
  ggplot(
    aes(x = difference * log2(exp(1)),
        y = -log10(padj),
        color = conclusion)
  ) +
  geom_point(size = 0.2) +
  theme_classic() +
  scale_color_manual(
    values = c("darkcyan", "darkorange", "darkgray")
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = -1,
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 1,
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  xlab("L2FC(kdeg)") +
  ylab("-log10(padj)") +
  theme(legend.position="none") +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) +
  coord_cartesian(
    xlim = c(-2.25, 2.25)
  )




DDX3X_kdeg_eC <- DDX3X_kdeg %>%
  left_join(
    annotated_peaks %>%
      dplyr::select(gene_id, peak_chr) %>%
      dplyr::distinct() %>%
      dplyr::rename(XF = gene_id),
    by = "XF"
  ) %>%
  dplyr::mutate(
    DDX3X_target = factor(!is.na(peak_chr),
                          levels = c(FALSE, TRUE))
  ) %>%
  arrange(desc(DDX3X_target))


kd_box <- DDX3X_kdeg_eC %>%
  ggplot(
    aes(x = DDX3X_target,
        y = difference * log2(exp(1)))
  ) +
  geom_jitter(width = 0.2,
              height = 0,
              size = 0.1,
              alpha = 0.1) +
  geom_boxplot(fill = "aquamarine3",
               color = "black",
               outliers = F,
               linewidth = 0.2,
               notch = T,
               width = 0.3) +
  theme_classic() +
  geom_hline(
    yintercept = median(DDX3X_kdeg_eC$difference[DDX3X_kdeg_eC$DDX3X_target == FALSE]),
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  xlab("DDX3X target (ENCODE eCLIP)") +
  ylab("L2FC(kdeg)") +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) +
  coord_cartesian(
    ylim = c(-2.25, 2.25)
  )


setwd(figure_savedir)
ggsave(
  filename = "kdeg_volcano_plot.pdf",
  plot = kd_volc,
  width = 2.5,
  height = 2
)

ggsave(
  filename = "kdeg_box_plot.pdf",
  plot = kd_box,
  width = 2.5,
  height = 2
)



##### Does DDX3X binding correlate with ksyn changes? #####

DDX3X_ksyn <- EZget(
  ezbdo,
  type = "comparisons",
  parameter = "logk1",
  experimental = "DDX3X"
)


ks_volc <- DDX3X_ksyn %>%
  mutate(
    L2FC = difference * log2(exp(1)),
    conclusion = factor(case_when(
      padj < 0.05 & L2FC < -1 ~ "Decreased",
      padj < 0.05 & L2FC > 1 ~ "Increased",
      .default = "Not sig."
    ),
    levels = c("Decreased", "Increased", "Not sig."))
  ) %>%
  ggplot(
    aes(x = difference * log2(exp(1)),
        y = -log10(padj),
        color = conclusion)
  ) +
  geom_point(size = 0.2) +
  theme_classic() +
  scale_color_manual(
    values = c("darkcyan", "darkorange", "darkgray")
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = -1,
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 1,
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  xlab("L2FC(ksyn)") +
  ylab("-log10(padj)") +
  theme(legend.position="none") +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) +
  coord_cartesian(
    xlim = c(-2.25, 2.25)
  )




DDX3X_ksyn_eC <- DDX3X_ksyn %>%
  left_join(
    annotated_peaks %>%
      dplyr::select(gene_id, peak_chr) %>%
      dplyr::distinct() %>%
      dplyr::rename(XF = gene_id),
    by = "XF"
  ) %>%
  dplyr::mutate(
    DDX3X_target = factor(!is.na(peak_chr),
                          levels = c(FALSE, TRUE))
  ) %>%
  arrange(desc(DDX3X_target))


ks_box <- DDX3X_ksyn_eC %>%
  ggplot(
    aes(x = DDX3X_target,
        y = difference * log2(exp(1)))
  ) +
  geom_jitter(width = 0.2,
              height = 0,
              size = 0.1,
              alpha = 0.1) +
  geom_boxplot(fill = "aquamarine3",
               color = "black",
               outliers = F,
               linewidth = 0.2,
               notch = T,
               width = 0.3) +
  theme_classic() +
  geom_hline(
    yintercept = median(DDX3X_ksyn_eC$difference[DDX3X_ksyn_eC$DDX3X_target == FALSE]),
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.4
  ) +
  xlab("DDX3X target (ENCODE eCLIP)") +
  ylab("L2FC(ksyn)") +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) +
  coord_cartesian(
    ylim = c(-2.25, 2.25)
  )


setwd(figure_savedir)
ggsave(
  filename = "ksyn_volcano_plot.pdf",
  plot = ks_volc,
  width = 2.5,
  height = 2
)

ggsave(
  filename = "ksyn_box_plot.pdf",
  plot = ks_box,
  width = 2.5,
  height = 2
)


# Reproduce full analysis ------------------------------------------------------



### Create EZbakRData object
metadf <- dplyr::tibble(
  sample = c("SRR20080253", "SRR20080252", "SRR20080251", "SRR20080250",
             "SRR20080249", "SRR20080248", "SRR20080246", "SRR20080245",
             "SRR20080244", "SRR20080247", "SRR20080243", "SRR20080242",
             "SRR20080241", "SRR20080240", "SRR20080239", "SRR20080238",
             "SRR20080237", "SRR20080235", "SRR20080234", "SRR20080233",
             "SRR20080236", "SRR20080232", "SRR20080231", "SRR20080230",
             "SRR20080229", "SRR20080227", "SRR20080226", "SRR20080225",
             "SRR20080228", "SRR20080224"),
  tl = rep(c(0, 0, 60, 60, 60, 60, 60, 60, 60, 60),
           times = 3),
  siRNA = rep(c("scramble", "DDX3X", "PABPC4"),
              each = 10),
  compartment = rep(rep(c("total", "chromatin", "nuclear", "cytoplasm"),
                        times = c(4, 2, 2, 2)), times = 3)
) %>%
  dplyr::filter(compartment != "chromatin")


ds <- open_dataset(arrow_dataset)

ezbdo <- EZbakRArrowData(ds, metadf)

ezbdo <- EstimateFractions(ezbdo, pold_from_nolabel = TRUE,
                           features = "XF")



##### Normalize for dropout #####
# Now implemented in EZbakR's NormalizeForDropout() function

fractions <- ezbdo$fractions$XF


do_norm_ll <- function(param,
                       reffn,
                       fndo,
                       sig){

  pdo <- inv_logit(param[1])

  Efn <- fndo / ((1 - pdo) + (fndo*pdo))

  ll <- stats::dnorm(EZbakR:::logit(reffn),
                     EZbakR:::logit(Efn),
                     sig,
                     log = TRUE)

  return(-sum(ll))

}

ref_samples <- fractions %>%
  inner_join(metadf,
             by = "sample") %>%
  dplyr::filter(tl > 0) %>%
  dplyr::group_by(
    siRNA, compartment, sample
  ) %>%
  dplyr::filter(n > 50) %>%
  dplyr::summarise(
    avg_fn = mean(fraction_highTC)
  ) %>%
  dplyr::arrange(
    compartment,
    avg_fn
  ) %>%
  ungroup() %>%
  group_by(compartment) %>%
  filter(avg_fn == max(avg_fn)) %>%
  ungroup() %>%
  dplyr::select(sample, compartment)


fxn_wide <- fractions %>%
  dplyr::filter(n > 50) %>%
  dplyr::select(sample, XF, fraction_highTC) %>%
  pivot_wider(
    names_from = "sample",
    values_from = "fraction_highTC"
  )

pdos <- metadf %>%
  dplyr::filter(tl > 0) %>%
  dplyr::mutate(
    pdo = 0
  )

for(i in seq_along(pdos$sample)){

  samp <- pdos$sample[i]

  if(samp %in% ref_samples$sample) next

  reference <- ref_samples$sample[ref_samples$compartment == pdos$compartment[i]]

  pdos$pdo[i] <- fxn_wide %>%
    dplyr::select(
      XF,
      !!samp,
      !!reference
    ) %>%
    na.omit() %>%
    summarise(
      pdo = stats::optim(
        par = c(-2),
        fn = do_norm_ll,
        reffn = !!dplyr::sym(reference),
        fndo = !!dplyr::sym(samp),
        sig = 0.2,
        method = "L-BFGS-B",
        upper = 6,
        lower = -6
      )$par[1]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(pdo) %>%
    unlist() %>%
    unname() %>%
    EZbakR:::inv_logit()


}

fractions_normalized <- fractions %>%
  dplyr::left_join(
    pdos %>%
      dplyr::select(sample, pdo),
    by = "sample"
  ) %>%
  dplyr::mutate(
    pdo = ifelse(is.na(pdo), 0, pdo)
  ) %>%
  dplyr::mutate(
    fraction_highTC = (fraction_highTC)/((1 - pdo) + fraction_highTC * pdo),
  ) %>%
  dplyr::mutate(
    logit_fraction_highTC = EZbakR:::logit(fraction_highTC)
  ) %>%
  dplyr::group_by(
    sample
  ) %>%
  dplyr::mutate(
    num = mean(fraction_highTC)*(1-pdo) + (1-mean(fraction_highTC)),
    den = fraction_highTC*(1-pdo) + (1 - fraction_highTC),
    n = n*(num/den)
  ) %>%
  dplyr::select(
    -num, -den
  )



ezbdo$fractions$XF <- fractions_normalized

##### Analyze kinetics #####

ezbdo <- AverageAndRegularize(ezbdo,
                              parameter = "logit_fraction_highTC",
                              type = "fractions",
                              formula_mean = ~compartment:siRNA - 1)


graph <- matrix(c(0, 1, 0,
                  0, 0, 2,
                  3, 0, 0),
                nrow = 3,
                ncol = 3,
                byrow = TRUE)
rownames(graph) <- c("0", "N", "C")
colnames(graph) <- rownames(graph)


modeled_to_measured <- list(
  nuclear = list(XF ~ N),
  cytoplasm = list(XF ~ C),
  total = list(XF ~ C + N) # total RNA is a combination of C and N
)


ezbdo <- EZDynamics(ezbdo,
                    graph = graph,
                    sub_features = "XF",
                    grouping_features = "XF",
                    sample_feature = "compartment",
                    modeled_to_measured = modeled_to_measured)

saveRDS(ezbdo, ezbdo_path)

