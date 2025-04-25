### PURPOSE OF THIS SCRIPT
## Reproduce Figure S6 (comparison with kinetic parameters from
## Ietswaart et al. 2024)


# Load dependencies ------------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(devtools)
load_all("C:/Users/isaac/Documents/Simon_Lab/EZbakR/")

# Fit full model to subcellular TL-seq data ------------------------------------

### STEP 1: ESTIMATE FRACTIONS

# Which feature to use?
feature_to_use <- "XF"
chr_list <- list(XF ~ CH)
nuc_list <- list(XF ~ NP + CH)
cyt_list <- list(XF ~ CY + PL)
pl_list <- list(XF ~ PL)
tot_list <- list(XF ~ CH + CY + NP + PL)


ds <- open_dataset("C:/Users/")

comps <- c("chr", "cyto", "nuc", "poly", "total")
times <- c("0min", "15min", "30min", "60min", "120min")
metadf <- tibble(sample = paste0("K562_",
                                 rep(comps, each = 10),
                                 "_WT_",
                                 rep(rep(times, each = 2), times = 5),
                                 rep(c("_rep1", "_rep2"), times = 25))) %>%
  mutate(compartment = case_when(
    grepl("_chr_", sample) ~ "chromatin",
    grepl("_nuc_", sample) ~ "nucleus",
    grepl("_poly_", sample) ~ "polysome",
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


### STEP 4: ESTIMATE NORMALZIATION FACTORS


# graph
graph <- matrix(c(0, 1, 0, 0, 0,
                  0, 0, 2, 0, 0,
                  0, 0, 0, 3, 0,
                  4, 0, 0, 0, 5,
                  4, 0, 0, 0, 0),
                nrow = 5,
                ncol = 5,
                byrow = TRUE)

colnames(graph) <- c("0", "CH", "NP", "CY", "PL")
rownames(graph) <- colnames(graph)


### STEP 5: FIT MODEL

# With separate compartments as separate samples
ezbado_c <- AverageAndRegularize(ezbado_c,
                                 formula_mean = ~tl:compartment - 1,
                                 type = "fractions",
                                 parameter = "logit_fraction_highTC")


# Estimate parameters
ezbado_c <- EZDynamics(ezbado_c,
                       graph = graph,
                       sub_features = feature_to_use,
                       grouping_features = feature_to_use,
                       sample_feature = "compartment",
                       modeled_to_measured = list(chromatin = chr_list,
                                                  nucleus = nuc_list,
                                                  cytoplasm = cyt_list,
                                                  polysome = pl_list,
                                                  total = tot_list)
)


setwd("G:/Shared drives/Matthew_Simon/IWV/EZbakR_paper/Fits/Subcellular_TL/")
saveRDS(ezbado_c,
        file= "XF_fullDOC_fnonly_subtlseq_fit_20240910.rds")


# Figure S6 --------------------------------------------------------------------

churchman_ests <- readxl::read_excel("C:/Users/isaac/Documents/Simon_Lab/EZbakR_paper/Data/Figure_4/Churchman_estimates/mmc2.xlsx",
                                     sheet = 2)

setwd("C:/Users/isaac/Yale University/Simon Lab – RNA - Documents/IWV/EZbakR_paper/Fits/Subcellular_TL/")
ez_ests <- readRDS("XF_fullDOC_fnonly_subtlseq_fit_20240910.rds")
my_ests <- ez_ests$dynamics$dynamics1

if("GF" %in% colnames(my_ests)){

  my_ests <- my_ests %>%
    dplyr::rename(XF = GF)

}

if("k2_se" %in% colnames(my_ests)){

  my_ests <- my_ests %>%
    dplyr::rename(se_k2 = k2_se,
                  se_k3 = k3_se,
                  se_k4 = k4_se,
                  se_k5 = k5_se)

}


### Uncertainties

u_violin <- my_ests %>%
  dplyr::select(-k1, -se_k1) %>%
  dplyr::select(starts_with("se_")) %>%
  tidyr::pivot_longer(
    everything(),
    values_to = "uncertainty",
    names_to = "parameter"
  ) %>%
  dplyr::mutate(
    parameter = factor(case_when(
      parameter == "se_k2" ~ "kch",
      parameter == "se_k3" ~ "kexp",
      parameter == "se_k4" ~ "kdeg",
      parameter == "se_k5" ~ "kpl"
    ),
    levels = c("kch", "kexp", "kdeg", "kpl")
    )
  ) %>%
  dplyr::filter(
    uncertainty < 1.5
  ) %>%
  ggplot(
    aes(x = parameter,
        y = log(uncertainty),
        fill = parameter)
  ) +
  geom_jitter(
    width = 0.3,
    size = 0.1,
    alpha = 0.1,
    height = 0
  ) +
  geom_violin(color = "black",
              alpha = 1,
              linewidth = 0.5) +
  theme_classic() +
  scale_fill_manual(
    values = c("#6699CC",
               "gray60",
               "#CC6666",
               "#CC99FF")
  ) +
  xlab("Parameter") +
  ylab("log(Uncertainty)") +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10),
    legend.position = "none")



### Polysome loading rate

my_polys <- my_ests %>%
  dplyr::select(XF, k5, se_k5) %>%
  dplyr::filter(se_k5 < 1) %>%
  dplyr::rename(Gene = XF,
                my_kpoly = k5)

church_polys <- churchman_ests %>%
  dplyr::select(Gene, Symbol, PUND, half_life_poly_entry.Mean) %>%
  na.omit() %>%
  mutate(church_kpoly = log(2)/half_life_poly_entry.Mean)

compare_polys <- my_polys %>%
  inner_join(
    church_polys,
    by = c("Gene")
  )


kpl_plot <- compare_polys %>%
  dplyr::filter(PUND == "FALSE") %>%
  mutate(density = get_density(
    x = log(church_kpoly),
    y = my_kpoly,
    n = 200
  )) %>%
  ggplot(aes(x = log(church_kpoly),
             y = my_kpoly,
             color = density)) +
  geom_point(
    size = 0.2
  ) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(kpl churchman") +
  ylab("log(kpl EZbakR") +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = 'darkred',
    linetype = 'dotted',
    linewidth = 0.4
  ) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10),
    legend.position = "none")



### Cytoplasmic degradation rate

my_cyts <- my_ests %>%
  dplyr::select(XF, k4, se_k4) %>%
  dplyr::filter(se_k4 < 0.3) %>%
  dplyr::rename(Gene = XF,
                my_kcyto = k4)

church_cyts <- churchman_ests %>%
  dplyr::select(Gene, Symbol, PUND, half_life_cyto.Mean) %>%
  na.omit() %>%
  mutate(church_kcyto = log(2)/half_life_cyto.Mean)

compare_cyts <- my_cyts %>%
  inner_join(
    church_cyts,
    by = c("Gene")
  )


kdeg_plot <- compare_cyts %>%
  dplyr::filter(PUND == "FALSE") %>%
  mutate(density = get_density(
    x = log(church_kcyto),
    y = my_kcyto,
    n = 200
  )) %>%
  ggplot(aes(x = log(church_kcyto),
             y = my_kcyto,
             color = density)) +
  geom_point(
    size = 0.2
  ) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(kdeg churchman") +
  ylab("log(kdeg EZbakR") +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = 'darkred',
    linetype = 'dotted',
    linewidth = 0.4
  ) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10),
    legend.position = "none")


### Nuclear export rate

my_exps <- my_ests %>%
  dplyr::select(XF, k3, se_k3) %>%
  dplyr::filter(se_k3 < 0.3) %>%
  dplyr::rename(Gene = XF,
                my_kexp = k3)

church_exps <- churchman_ests %>%
  dplyr::select(Gene, Symbol, PUND, half_life_nucexp_from_nucres.Mean) %>%
  na.omit() %>%
  mutate(church_kexp = log(2)/half_life_nucexp_from_nucres.Mean)

compare_exps <- my_exps %>%
  inner_join(
    church_exps,
    by = c("Gene")
  )


kexp_plot <- compare_exps %>%
  dplyr::filter(PUND == "FALSE") %>%
  mutate(density = get_density(
    x = log(church_kexp),
    y = my_kexp,
    n = 200
  )) %>%
  ggplot(aes(x = log(church_kexp),
             y = my_kexp,
             color = density)) +
  geom_point(
    size = 0.2
  ) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(kexp churchman") +
  ylab("log(kexp EZbakR") +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = 'darkred',
    linetype = 'dotted',
    linewidth = 0.4
  ) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10),
    legend.position = "none")



### Chromatin off rate comparison
# Looks pretty darn similar
# Differences likely due to dropout correction
# and model differences

my_chrs <- my_ests %>%
  dplyr::select(XF, k2, se_k2) %>%
  dplyr::filter(se_k2 < 0.3) %>%
  dplyr::rename(Gene = XF,
                my_kchr = k2)

church_chrs <- churchman_ests %>%
  dplyr::select(Gene, Symbol, PUND, half_life_chr.Mean) %>%
  na.omit() %>%
  mutate(church_kchr = log(2)/half_life_chr.Mean)

compare_ests <- my_chrs %>%
  inner_join(
    church_chrs,
    by = c("Gene")
  )


kch_plot <- compare_ests %>%
  mutate(density = get_density(
    x = log(church_kchr),
    y = my_kchr,
    n = 200
  )) %>%
  ggplot(aes(x = log(church_kchr),
             y = my_kchr,
             color = density)) +
  geom_point(
    size = 0.2
  ) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(kch churchman") +
  ylab("log(kch EZbakR") +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = 'darkred',
    linetype = 'dotted',
    linewidth = 0.4
  ) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10),
    legend.position = "none")




### Save plots
setwd("C:/Users/isaac/Documents/Simon_Lab/EZbakR_paper/Figures/Supplemental_subtlseq/")
ggsave(
  filename = "kch_comparison.pdf",
  plot = kch_plot,
  width = 2.25,
  height = 2
)
ggsave(
  filename = "kexp_comparison.pdf",
  plot = kexp_plot,
  width = 2.25,
  height = 2
)
ggsave(
  filename = "kdeg_comparison.pdf",
  plot = kdeg_plot,
  width = 2.25,
  height = 2
)
ggsave(
  filename = "kpl_comparison.pdf",
  plot = kpl_plot,
  width = 2.25,
  height = 2
)
ggsave(
  filename = "uncertainty_comparison.pdf",
  plot = u_violin,
  width = 4.5,
  height = 2
)
