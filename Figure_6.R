### PURPOSE OF THIS SCRIPT
## Reproduce Figure 6 of EZbakR-suite paper.
##
## Workflow:
## 1) Analyze simulated data from bakR paper with EZbakR
## 2) Compare performance of EZbakR and bakR implementations
##
## Required data:
## 1) Simulated datasets from bakR paper (5 replicates)
##
## Estimated total runtime: ~30 minutes



# Load dependencies ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(data.table)
library(readr)
library(EZbakR)
library(bakR)


# Directory to save figures to
savedir <- getwd()

# Directory containing all replicates of simulated data
simdata_dir <- "G:/Shared drives/Matthew_Simon/IWV/bakR_Paper/Data/Fits/Simulation_Replicates/2_reps/"


# Use to get Matthew's Correlation Coefficient
calc_MCC <- function(N, P, FDR, TPR){

  TP <- TPR*P
  FP <- (FDR*(TP))/(1-FDR)
  PP <- FP + TP
  FN <- P - TP
  TN <- N - FP
  PN <- FN + TN

  TNR <- TN/N
  NPV <- TN/PN
  FNR <- FN/P
  FPR <- FP/N

  #MCC <-sqrt(TPR*TNR*(1-FDR)*NPV) - sqrt(FNR*FPR*(1-NPV)*FDR)

  MCC <- (TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))


}


# Analyze bakR simulated data (moderate variance trend) ------------------------

rep_names <- c("1st", "2nd", "3rd", "4th", "5th")
final_df <- tibble()
for(i in 1:5){


  ##### ASSESS EZBAKR PERFORMANCE

  simdata <- readRDS(paste0(simdata_dir, "/Rep_", i,"/MLE_fit_2_reps_",  rep_names[i],"_rep_2022_04_01.rds"))
  mcmc_sim <- readRDS(paste0(simdata_dir, "/Rep_", i,"/Fit_MCMC_2reps_",  rep_names[i],"_rep_2022_04_01.rds"))

  cB <- as_tibble(simdata$Data_lists$Fast_df) %>%
    group_by(sample, XF, TC, nT) %>%
    summarise(n = sum(n))
  metadf <- simdata$Data_lists$Fast_df %>%
    dplyr::select(sample, mut, tl) %>%
    dplyr::distinct() %>%
    dplyr::rename(Exp_ID = mut)
  metadf$tl <- c(60, 60, 60, 60, 0, 0)

  ezbdo <- EZbakRData(cB,
                      metadf)

  ezbdo <- EstimateFractions(ezbdo)
  ezbdo <- EstimateKinetics(ezbdo)
  ezbdo <- AverageAndRegularize(ezbdo)
  ezbdo <- CompareParameters(ezbdo,
                             design_factor = "Exp_ID",
                             reference = "1",
                             experimental = "2")


  comparisons <- ezbdo$comparisons$comparison1

  hyperparam <- simdata$Fast_Fit$Hyper_Parameters[1]
  mle_res <- simdata$Fast_Fit$Effects_df %>%
    dplyr::mutate(pval = 2*pt(-abs(effect/se), df = 2 + 2*hyperparam)) %>%
    ungroup() %>%
    dplyr::mutate(padj = p.adjust(pval, method = "BH"))


  mcmc_res <- mcmc_sim$Stan_Fit$Effects_df

  cutoffs <- seq(from = 0, to = 1, length.out = 10000)
  Powers <- rep(0, times = length(cutoffs))
  FDRs <- rep(0, times = length(cutoffs))
  Powers_mle <- Powers
  FDRs_mle <- FDRs
  MCCs <- MCCs_mle <- MCCs_mcmc <- Powers_mcmc <- FDRs_mcmc <- FDRs


  nsig <- 1000
  nnotsig <- 9000

  for(c in seq_along(cutoffs)){

    sig_genes <- comparisons %>%
      filter(padj < cutoffs[c]) %>%
      dplyr::select(XF) %>%
      unlist() %>%
      unname() %>%
      as.numeric()

    sig_genes_mle <- mle_res %>%
      filter(padj < cutoffs[c]) %>%
      dplyr::select(XF) %>%
      unlist() %>%
      unname() %>%
      as.numeric()


    sig_genes_mcmc <- mcmc_res %>%
      filter(padj < cutoffs[c]) %>%
      dplyr::select(XF) %>%
      unlist() %>%
      unname() %>%
      as.numeric()


    if(length(sig_genes) == 0){

      FDRs[c] <- 0

    }else{

      FDRs[c] <- sum(sig_genes <= nnotsig)/length(sig_genes)


    }

    if(length(sig_genes_mle) == 0){

      FDRs_mle[c] <- 0

    }else{

      FDRs_mle[c] <- sum(sig_genes_mle <= nnotsig)/length(sig_genes_mle)

    }

    if(length(sig_genes_mcmc) == 0){

      FDRs_mcmc[c] <- 0

    }else{

      FDRs_mcmc[c] <- sum(sig_genes_mcmc <= nnotsig)/length(sig_genes_mcmc)

    }


    Powers[c] <- sum(sig_genes > nnotsig)/nsig
    Powers_mle[c] <- sum(sig_genes_mle > nnotsig)/nsig
    Powers_mcmc[c] <- sum(sig_genes_mcmc > nnotsig)/nsig

    MCCs[c] <- calc_MCC(N = nnotsig,
                        P = nsig,
                        FDR = FDRs[c],
                        TPR = Powers[c])

    MCCs_mle[c] <- calc_MCC(N = nnotsig,
                            P = nsig,
                            FDR = FDRs_mle[c],
                            TPR = Powers_mle[c])

    MCCs_mcmc[c] <- calc_MCC(N = nnotsig,
                             P = nsig,
                             FDR = FDRs_mcmc[c],
                             TPR = Powers_mcmc[c])

  }

  results_df <- tibble(FDR_ezbakR = FDRs,
                       Power_ezbakR = Powers,
                       MCC_ezbakR = MCCs,
                       FDR_mle = FDRs_mle,
                       Power_mle = Powers_mle,
                       MCC_mle = MCCs_mle,
                       FDR_mcmc = FDRs_mcmc,
                       Power_mcmc = Powers_mcmc,
                       MCC_mcmc = MCCs_mcmc,
                       cutoff = cutoffs,
                       Replicate = i)
  final_df <- bind_rows(final_df, results_df)

}


setwd(savedir)
write_csv(final_df,
          "Two_reps_final_df.csv")



##### START HERE TO AVOID RERUNNING ANALYSIS


### MCC plot

cutoffs <- final_df %>%
  filter(cutoff > 0.049 & cutoff < 0.051) %>%
  dplyr::select(cutoff) %>%
  unlist() %>%
  unname()

if(length(cutoffs) %% 2 == 1){
  cutoff_touse <- median(cutoffs)
}else{
  index <- round(length(cutoffs)/2)
  cutoff_touse <- cutoffs[index]
}


MCCs <- final_df %>%
  filter(cutoff == cutoff_touse) %>%
  dplyr::select(Replicate, starts_with("MCC")) %>%
  tidyr::pivot_longer(cols = starts_with("MCC"),
                      names_to = "model",
                      values_to = "MCC",
                      names_pattern = "MCC_(.*)") %>%
  mutate(model = case_when(
    model == "ezbakR" ~ "EZbakR",
    model == "mle" ~ "bakR (MLE)",
    model == "mcmc" ~ "bakR (MCMC)"
  ))

avg_MCCs <- final_df %>%
  filter(cutoff == cutoff_touse) %>%
  dplyr::select(Replicate, starts_with("MCC")) %>%
  tidyr::pivot_longer(cols = starts_with("MCC"),
                      names_to = "model",
                      values_to = "MCC",
                      names_pattern = "MCC_(.*)") %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(avg_MCC = mean(MCC),
                   sd_MCC = sd(MCC)) %>%
  mutate(model = case_when(
    model == "ezbakR" ~ "EZbakR",
    model == "mle" ~ "bakR (MLE)",
    model == "mcmc" ~ "bakR (MCMC)"
  ))



set.seed(44)
gMCC <- ggplot(MCCs,
               aes(x = model, y = MCC, color = model)) +
  geom_errorbar(data = avg_MCCs, mapping = aes(x = model, ymin = avg_MCC - 2*sd_MCC/sqrt(5), ymax = avg_MCC + 2*sd_MCC/sqrt(5)),
                width = 0.4, inherit.aes = FALSE,
                linewidth = 0.2) +
  geom_point(data = avg_MCCs, mapping = aes(x = model, y = avg_MCC), color = "black", inherit.aes = FALSE,
             size = 1)  +
  geom_jitter( height = 0, width = 0.4, size = 0.8,
               color = "gray50") +  theme_classic() +
  xlab("Model") +
  ylab("MCC") +
  ylim(c(0.35, 0.5))  +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10), #change font size of legend title
    axis.text.x = element_text(angle = 45, vjust = 0.6, hjust=0.5)) +
  theme(legend.position = "none")

gMCC

### Make Power plot as an MCC-like scatter plot

cutoffs <- final_df %>%
  filter(cutoff > 0.049 & cutoff < 0.051) %>%
  dplyr::select(cutoff) %>%
  unlist() %>%
  unname()

if(length(cutoffs) %% 2 == 1){
  cutoff_touse <- median(cutoffs)
}else{
  index <- round(length(cutoffs)/2)
  cutoff_touse <- cutoffs[index]
}


Powers <- final_df %>%
  filter(cutoff == cutoff_touse) %>%
  dplyr::select(Replicate, starts_with("Power")) %>%
  tidyr::pivot_longer(cols = starts_with("Power"),
                      names_to = "model",
                      values_to = "Power",
                      names_pattern = "Power_(.*)") %>%
  mutate(model = case_when(
    model == "ezbakR" ~ "EZbakR",
    model == "mle" ~ "bakR (MLE)",
    model == "mcmc" ~ "bakR (MCMC)"
  ))

avg_Powers <- final_df %>%
  filter(cutoff == cutoff_touse) %>%
  dplyr::select(Replicate, starts_with("Power")) %>%
  tidyr::pivot_longer(cols = starts_with("Power"),
                      names_to = "model",
                      values_to = "Power",
                      names_pattern = "Power_(.*)") %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(avg_Power = mean(Power),
                   sd_Power = sd(Power)) %>%
  mutate(model = case_when(
    model == "ezbakR" ~ "EZbakR",
    model == "mle" ~ "bakR (MLE)",
    model == "mcmc" ~ "bakR (MCMC)"
  ))



set.seed(44)
gPower_s <- ggplot(Powers,
                   aes(x = model, y = Power)) +
  geom_errorbar(data = avg_Powers, mapping = aes(x = model, ymin = avg_Power - 2*sd_Power/sqrt(5),
                                                 ymax = avg_Power + 2*sd_Power/sqrt(5)),
                width = 0.4, inherit.aes = FALSE,
                linewidth = 0.2) +
  geom_point(data = avg_Powers, mapping = aes(x = model, y = avg_Power), color = "black", inherit.aes = FALSE,
             size = 1)  +
  geom_jitter( height = 0, width = 0.4, size = 0.8,
               color = "gray50") +
  theme_classic() +
  xlab("Model") +
  ylab("Power") +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10), #change font size of legend title
    axis.text.x = element_text(angle = 45, vjust = 0.6, hjust=0.5)) +
  theme(legend.position = "none") +
  ylim(c(0, 0.3))

gPower_s


### Make FDR plot as an MCC-like scatter plot

cutoffs <- final_df %>%
  filter(cutoff > 0.049 & cutoff < 0.051) %>%
  dplyr::select(cutoff) %>%
  unlist() %>%
  unname()

if(length(cutoffs) %% 2 == 1){
  cutoff_touse <- median(cutoffs)
}else{
  index <- round(length(cutoffs)/2)
  cutoff_touse <- cutoffs[index]
}


FDRs <- final_df %>%
  filter(cutoff == cutoff_touse) %>%
  dplyr::select(Replicate, starts_with("FDR")) %>%
  tidyr::pivot_longer(cols = starts_with("FDR"),
                      names_to = "model",
                      values_to = "FDR",
                      names_pattern = "FDR_(.*)") %>%
  mutate(model = case_when(
    model == "ezbakR" ~ "EZbakR",
    model == "mle" ~ "bakR (MLE)",
    model == "mcmc" ~ "bakR (MCMC)"
  ))

avg_FDRs <- final_df %>%
  filter(cutoff == cutoff_touse) %>%
  dplyr::select(Replicate, starts_with("FDR")) %>%
  tidyr::pivot_longer(cols = starts_with("FDR"),
                      names_to = "model",
                      values_to = "FDR",
                      names_pattern = "FDR_(.*)") %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(avg_FDR = mean(FDR),
                   sd_FDR = sd(FDR)) %>%
  mutate(model = case_when(
    model == "ezbakR" ~ "EZbakR",
    model == "mle" ~ "bakR (MLE)",
    model == "mcmc" ~ "bakR (MCMC)"
  ))



set.seed(44)
gFDR_s <- ggplot(FDRs,
                 aes(x = model, y = FDR)) +
  geom_errorbar(data = avg_FDRs, mapping = aes(x = model, ymin = avg_FDR - 2*sd_FDR/sqrt(5),
                                               ymax = avg_FDR + 2*sd_FDR/sqrt(5)),
                width = 0.4, inherit.aes = FALSE,
                linewidth = 0.2) +
  geom_point(data = avg_FDRs, mapping = aes(x = model, y = avg_FDR), color = "black", inherit.aes = FALSE,
             size = 1)  +
  geom_jitter( height = 0, width = 0.4, size = 0.8,
               color = "gray50") +
  xlab("Model") +
  ylab("FDR") +
  theme_classic() +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10), #change font size of legend title
    axis.text.x = element_text(angle = 45, vjust = 0.6, hjust=0.5)) +
  theme(legend.position = "none") +
  ylim(c(0, 0.06)) +
  geom_hline(yintercept = 0.05,
             color = 'darkred',
             linetype = 'dotted',
             linewidth = 0.5)

### Save plots

setwd(savedir)
ggsave(plot = gMCC,
       filename = "MCC_2reps_rerun.pdf",
       width = 2.15,
       height = 2.35)
ggsave(plot = gPower_s,
       filename = "Power_scatter_2reps_rerun.pdf",
       width = 2.15,
       height = 2.35)
ggsave(plot = gFDR_s,
       filename = "FDR_scatter_2reps_rerun.pdf",
       width = 2.15,
       height = 2.35)



# Collect runtime data ---------------------------------------------------------

trials <- 5

final_time <- tibble()
for(i in 1:trials){

  simdata <- EZSimulate(10000, nreps = 2)


  start_ez <- Sys.time()

  ezbdo <- EZbakRData(simdata$cB, simdata$metadf)
  ezbdo <- EstimateFractions(ezbdo)
  ezbdo <- EstimateKinetics(ezbdo)
  ezbdo <- AverageAndRegularize(ezbdo)
  ezbdo <- CompareParameters(ezbdo,
                             condition = 'treatment',
                             reference = 'treatment1',
                             experimental = 'treatment2')

  end_ez <- Sys.time()



  start_b <- Sys.time()

  bmeta <- data.frame(
    tl = c(2, 2, 2, 2, 0, 0),
    Exp_ID = c(1, 1, 2, 2, 1, 2)
  )
  rownames(bmeta) <- simdata$metadf$sample
  bdo <- bakRData(simdata$cB %>%
                    rename(XF = feature), bmeta)
  bfo <- bakRFit(bdo)

  end_b <- Sys.time()


  time_df <- tibble(bakR_time = end_b - start_b,
                    EZbakR_time = end_ez - start_ez,
                    replicate = i)

  final_time <- bind_rows(final_time, time_df)


}


setwd(savedir)
write_csv(final_time, "EZbakR_and_bakR_runtimes.csv")



##### START HERE TO AVOID RERUNING

# MCMC time comes from runs on McCleary
final_time <- final_time %>%
  mutate(bakR_MCMC_time = c(105720,
                            105420,
                            89880,
                            86580,
                            83460),
         bakR_time = as.numeric(bakR_time),
         EZbakR_time = as.numeric(EZbakR_time)) %>%
  rename(Replicate = replicate)



Times <- final_time %>%
  dplyr::select(Replicate, ends_with("time")) %>%
  tidyr::pivot_longer(cols = ends_with("time"),
                      names_to = "model",
                      values_to = "Runtime",
                      names_pattern = "(.*)_time") %>%
  mutate(model = case_when(
    model == "EZbakR" ~ "EZbakR",
    model == "bakR" ~ "bakR (MLE)",
    model == "bakR_MCMC" ~ "bakR (MCMC)"
  ))

avg_Times <- final_time %>%
  dplyr::select(Replicate, ends_with("time")) %>%
  tidyr::pivot_longer(cols = ends_with("time"),
                      names_to = "model",
                      values_to = "Runtime",
                      names_pattern = "(.*)_time") %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(avg_runtime = mean(Runtime),
                   sd_runtime = sd(Runtime)) %>%
  mutate(model = case_when(
    model == "EZbakR" ~ "EZbakR",
    model == "bakR" ~ "bakR (MLE)",
    model == "bakR_MCMC" ~ "bakR (MCMC)"
  ))



set.seed(4)
gTime <- ggplot(Times,
                aes(x = model, y = Runtime, color = model)) +
  geom_errorbar(data = avg_Times, mapping = aes(x = model, ymin = avg_runtime - 2*sd_runtime/sqrt(5),
                                                ymax = avg_runtime + 2*sd_runtime/sqrt(5)),
                width = 0.4, inherit.aes = FALSE,
                linewidth = 0.2) +
  geom_point(data = avg_Times, mapping = aes(x = model, y = avg_runtime), color = "black", inherit.aes = FALSE,
             size = 1)  +
  geom_jitter( height = 0, width = 0.4, size = 0.8) +
  theme_classic() +
  xlab("Model") +
  ylab("Runtime (s)") +
  scale_color_manual(values = c("magenta", "darkgray", "darkgreen")) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10), #change font size of legend title
    axis.text.x = element_text(angle = 45, vjust = 0.6, hjust=0.5)) +
  theme(legend.position = "none") +
  scale_y_log10(breaks = scales::log_breaks(n = 5))



set.seed(421)
main_timeplot <- ggplot(Times,
                        aes(x = model, y = Runtime)) +
  geom_errorbar(data = avg_Times, mapping = aes(x = model, ymin = avg_runtime - 2*sd_runtime/sqrt(5),
                                                ymax = avg_runtime + 2*sd_runtime/sqrt(5)),
                width = 0.4, inherit.aes = FALSE,
                linewidth = 0.2) +
  geom_point(data = avg_Times, mapping = aes(x = model, y = avg_runtime), color = "black", inherit.aes = FALSE,
             size = 1)  +
  geom_jitter( height = 0, width = 0.4, size = 0.8,
               color = "gray50") +
  theme_classic() +
  xlab("Model") +
  ylab("Runtime (s)") +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10), #change font size of legend title
    axis.text.x = element_text(angle = 45, vjust = 0.6, hjust=0.5)) +
  theme(legend.position = "none")


set.seed(421)
inset_timeplot <- ggplot(Times %>%
                           filter(model != "bakR (MCMC)") %>%
                           mutate(model = factor(model,
                                                 levels = c("bakR (MLE)", "EZbakR"))),
                         aes(x = model, y = Runtime, color = model)) +
  geom_errorbar(data = avg_Times%>%
                  filter(model != "bakR (MCMC)") %>%
                  mutate(model = factor(model,
                                        levels = c("bakR (MLE)", "EZbakR"))),
                mapping = aes(x = model,
                              ymin = avg_runtime - 2*sd_runtime/sqrt(5),
                              ymax = avg_runtime + 2*sd_runtime/sqrt(5)),
                width = 0.4, inherit.aes = FALSE,
                linewidth = 0.2) +
  geom_point(data = avg_Times%>%
               filter(model != "bakR (MCMC)") %>%
               mutate(model = factor(model,
                                     levels = c("bakR (MLE)", "EZbakR"))),
             mapping = aes(x = model,
                           y = avg_runtime),
             color = "black",
             inherit.aes = FALSE,
             size = 1)  +
  geom_jitter( height = 0, width = 0.4, size = 0.8,
               color = "gray50") +
  theme_classic() +
  xlab("Model") +
  ylab("Runtime (s)") +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 75)) +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10), #change font size of legend title
    axis.text.x = element_text(angle = 45, vjust = 0.6, hjust=0.5)) +
  theme(legend.position = "none")



setwd(savedir)
ggsave(filename = "Comparing_runtimes_full.pdf",
       plot = main_timeplot,
       width = 2.25,
       height = 2.65)
ggsave(filename = "Comparing_runtimes_zoom.pdf",
       plot = inset_timeplot,
       width = 1.5,
       height = 2.08)

