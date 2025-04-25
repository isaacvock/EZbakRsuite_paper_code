### PURPOSE OF THIS SCRIPT
## Reproduce Figure S9, analysis of nanodynamo processed data
## and comparison to nanodynamo

# Load dependencies ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(dplyr)
library(devtools)
library(readxl)
library(EZbakR)


### File paths ###

# Path to nanodynamo RDS object
nd_data_path <- "C:/Users/isaac/Documents/Simon_Lab/EZbakR_paper/Data/nanodynamo/inferedRatesUntreatedMerged_YesChpNpNoP_multi.rds"

# Folder to save figures to
figure_savedir <- "C:/Users/isaac/Documents/Simon_Lab/EZbakR_paper/Figures/Supplement_nanodynamo/"

### Function(s) used throughout ###

process_table <- function(df, treatment_value){

  cols <- colnames(df)
  cols_to_keep <- cols[grepl("_", cols) &
                         (grepl("^[^_]*\\d[^_]*_", cols)
                          |
                            cols == "Entrez_ID"
                         )]


  fractions <- df %>%
    dplyr::select(!!cols_to_keep) %>%
    tidyr::pivot_longer(
      cols = !Entrez_ID,
      names_to = c("compartment", "RNA", "tl", "Replicate"),
      values_to = "reads",
      names_pattern = c("^(ch|cy|n|py)(mn|mp|pn|pp|n|p)(0|0\\.33)_(Rep[12])$")
    ) %>%
    dplyr::mutate(
      species = ifelse(
        RNA %in% c("pn", "pp"),
        "premRNA",
        "mRNA"
      )
    ) %>%
    group_by(Entrez_ID, compartment, tl, Replicate, species) %>%
    summarise(fraction_new = sum(reads[RNA %in% c("mn", "pn", "n")])
              /
                (sum(reads[RNA %in% c("mn", "pn", "n")]) +
                   sum(reads[RNA %in% c("mp", "pp", "p")])),
              reads = sum(reads)) %>%
    mutate(
      GF = ifelse(
        species == "mRNA",
        "__no_feature",
        Entrez_ID
      ),
      XF = ifelse(
        species == "mRNA",
        Entrez_ID,
        "__no_feature"
      ),
      fraction_new = ifelse(
        tl > 0 & fraction_new == 0 ,
        1 / (reads + 2),
        ifelse(fraction_new == 1 & tl > 0,
               (reads + 1)/(reads + 2),
               fraction_new)
      )
    ) %>%
    dplyr::select(-species) %>%
    dplyr::mutate(
      sample = paste0("sample_",compartment,"_", tl, "_", Replicate, "_treatment_", treatment_value)
    ) %>%
    ungroup()


  metadf <- fractions %>%
    dplyr::select(sample, tl, compartment) %>%
    dplyr::distinct() %>%
    dplyr::mutate(treatment = treatment_value)

  return(
    list(fractions = fractions,
         metadf = metadf)
  )

}

# Fit nanodynamo processed data ------------------------------------------------

ndata_ut_rds <- readRDS(nd_data_path)

# Generate data from these inferred using a negative binomial
# variance model
df <- ndata_ut_rds$inferedData


# 1. Generate the new data for Rep1
ndata_ut_rep1 <- as.data.frame(lapply(as_tibble(df), function(x){
  rnbinom(n = length(x), size = 100, mu = x)
} ))


# 2. Rename columns to append "Rep1"
colnames(ndata_ut_rep1) <- paste0(colnames(df), "_Rep1")

# 3. Generate the new data for Rep2
ndata_ut_rep2 <- as.data.frame(lapply(as_tibble(df), function(x) rnbinom(n = length(x), size = 100, mu = x)))

# 4. Rename columns to append "Rep2"
colnames(ndata_ut_rep2) <- paste0(colnames(df), "_Rep2")
ndata_ut_rep1$Entrez_ID <- rownames(df)
ndata_ut_rep2$Entrez_ID <- rownames(df)


ut_list <- process_table(
  inner_join(
    ndata_ut_rep1,
    ndata_ut_rep2,
    by = "Entrez_ID"
  ),
  "Untreated"
)

# Generalized from alternative runs
combined_fractions <- bind_rows(
  list(
    ut_list$fractions
  )
) %>%
  dplyr::select(sample, Entrez_ID, GF, XF, fraction_new, reads) %>%
  dplyr::rename(
    fraction_highTC = fraction_new,
    n = reads
  ) %>%
  dplyr::mutate(
    n = ceiling(n),
    logit_fraction_highTC = EZbakR:::logit(fraction_highTC)
  ) %>%
  na.omit()

combined_metadf <- bind_rows(
  list(
    ut_list$metadf
  )
) %>%
  mutate(tl = as.numeric(tl))

ezbfo <- EZbakRFractions(
  combined_fractions,
  combined_metadf
)



# Bit tricker; there is a mix of compartments with intronic and exonic data
# and compartments with just exonic data. Thus, have to tell EZbakR how many
# samples for each feature type is expected
feature_nsamp <- ezbfo$fractions$EntrezID_GF_XF %>%
  dplyr::select(Entrez_ID) %>%
  dplyr::distinct() %>%
  mutate(
    nsamps = sum(combined_metadf$tl > 0) +
      sum(combined_metadf$tl[combined_metadf$compartment %in% c("ch", "n")] > 0)
  )


ezbfo <- AverageAndRegularize(
  ezbfo,
  type = "fractions",
  parameter = "logit_fraction_highTC",
  #formula_mean = ~compartment:treatment - 1,
  formula_mean = ~compartment - 1,
  feature_sample_counts = feature_nsamp,
  scale_factor_df = tibble(
    sample = combined_metadf$sample,
    scale = 1
  )
)


### Define and fit ODE model

# Graphical representation of model
graph <- matrix(
  c(0, 1, 0, 0, 0, 0,
    0, 0, 2, 4, 0, 0,
    0, 0, 0, 0, 3, 0,
    0, 0, 0, 0, 5, 0,
    0, 0, 0, 0, 0, 6,
    7, 0, 0, 0, 0, 0),
  nrow = 6,
  ncol = 6,
  byrow = TRUE
)
rownames(graph) <- c("0", "CHP", "CHM", "NP", "NM", "C")
colnames(graph) <- rownames(graph)

modeled_to_measured <- list(
  ch = list(
    GF ~ CHP,
    XF ~ CHM
  ),
  n = list(
    GF ~ NP,
    XF ~ NM
  ),
  cy = list(
    XF ~ C
  )
)

# Fit model; nanodynamo developers claimed to have
# normalized data already. I don't think their normalization
# strategy is a good one, but for comparing EZbakR to nanodynamo,
# it's important that the same assumptions are made about the data.
# In addition, there is no better way to normalize the data given
# the model being fit and the data available.
ezbfo <- EZDynamics(ezbfo,
                    graph = graph,
                    sub_features = c("GF", "XF"),
                    grouping_features = c("Entrez_ID"),
                    sample_feature = "compartment",
                    modeled_to_measured = modeled_to_measured,
                    scale_factors = tibble(
                      compartment = c("ch", "n", "cy"),
                      scale = 1
                    ),
                    prior_means = c(3, rep(-3, times = max(graph) - 1)))



##### Compare to nanodynamo estimates


nd_rds <- readRDS(nd_data_path)


nd_ests <- nd_rds$inferedRates %>%
  as_tibble(
    rownames = "Entrez_ID"
  )


cols <- colnames(ezbfo$dynamics$dynamics1)
se_cols <- cols[grepl("^se_", cols)]
comp_ests <- ezbfo$dynamics$dynamics1 %>%
  filter(if_any(all_of(se_cols), ~ .x < 0.5)) %>%
  dplyr::inner_join(
    nd_ests,
    #by = c("Entrez_ID", "treatment")
    by = "Entrez_ID"
  )



k1plot <- comp_ests %>%
  na.omit() %>%
  dplyr::mutate(
    density = EZbakR:::get_density(
      x = log(k1),
      y = logk1,
      n = 200
    )
  ) %>%
  ggplot(
    aes(x = log(k1),
        y = logk1,
        color = density)
  ) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(nanodynamo k1)") +
  ylab("log(EZbkaR k1)") +
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






k2plot <- comp_ests %>%
  na.omit() %>%
  dplyr::mutate(
    density = EZbakR:::get_density(
      x = log(k2),
      y = logk2,
      n = 200
    )
  ) %>%
  ggplot(
    aes(x = log(k2),
        y = logk2,
        color = density)
  ) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(nanodynamo k2)") +
  ylab("log(EZbkaR k2)") +
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


k3plot <- comp_ests %>%
  na.omit() %>%
  dplyr::mutate(
    density = EZbakR:::get_density(
      x = log(k3),
      y = logk3,
      n = 200
    )
  ) %>%
  ggplot(
    aes(x = log(k3),
        y = logk3,
        color = density)
  ) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(nanodynamo k3)") +
  ylab("log(EZbkaR k3)") +
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



k4plot <- comp_ests %>%
  na.omit() %>%
  dplyr::mutate(
    density = EZbakR:::get_density(
      x = log(k4),
      y = logk4,
      n = 200
    )
  ) %>%
  ggplot(
    aes(x = log(k4),
        y = logk4,
        color = density)
  ) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(nanodynamo k4)") +
  ylab("log(EZbkaR k4)") +
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


k5plot <- comp_ests %>%
  na.omit() %>%
  dplyr::mutate(
    density = EZbakR:::get_density(
      x = log(k5),
      y = logk5,
      n = 200
    )
  ) %>%
  ggplot(
    aes(x = log(k5),
        y = logk5,
        color = density)
  ) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(nanodynamo k5)") +
  ylab("log(EZbkaR k5)") +
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



k6plot <- comp_ests %>%
  na.omit() %>%
  dplyr::mutate(
    density = EZbakR:::get_density(
      x = log(k6),
      y = logk6,
      n = 200
    )
  ) %>%
  ggplot(
    aes(x = log(k6),
        y = logk6,
        color = density)
  ) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(nanodynamo k6)") +
  ylab("log(EZbkaR k6)") +
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



k7plot <- comp_ests %>%
  na.omit() %>%
  dplyr::mutate(
    density = EZbakR:::get_density(
      x = log(k8),
      y = logk7,
      n = 200
    )
  ) %>%
  ggplot(
    aes(x = log(k8),
        y = logk7,
        color = density)
  ) +
  geom_point(size=0.4) +
  theme_classic() +
  scale_color_viridis_c() +
  xlab("log(nanodynamo k7)") +
  ylab("log(EZbkaR k7)") +
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


setwd(figure_savedir)
ggsave(
  filename = "k1plot_noP.pdf",
  plot = k1plot,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "k2plot_noP.pdf",
  plot = k2plot,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "k3plot_noP.pdf",
  plot = k3plot,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "k4plot_noP.pdf",
  plot = k4plot,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "k5plot_noP.pdf",
  plot = k5plot,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "k6plot_noP.pdf",
  plot = k6plot,
  width = 2,
  height = 1.67
)
ggsave(
  filename = "k7plot_noP.pdf",
  plot = k7plot,
  width = 2,
  height = 1.67
)
