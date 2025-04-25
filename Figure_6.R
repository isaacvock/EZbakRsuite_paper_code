### PURPOSE OF THIS SCRIPT
## Reproduce Figure 6 of EZbakR-suite paper.
##
## Workflow:
## 1) Analyze simulated data from bakR paper with EZbakR
## 2) Compare performance of EZbakR, bakR, and grandR
##
## Required data:
## 1) Simulated datasets from bakR paper (5 replicates)
## 2) grandR output for all simulated data
##
## Estimated total runtime: ~30 minutes



# Load dependencies ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(data.table)
library(readr)
library(EZbakR)
library(bakR)
library(grandR)


### File paths for panels B-E

# Path to save files and figures to
EZbakR_and_bakR_results <- "C:/Users/isaac/Documents/Simon_Lab/EZbakR_paper/Data/Figure_5/Two_reps_final_df_with_grandR.csv"


# Directory containing all replicates of simulated data
simdata_dir <- "C:/Users/isaac/Documents/Simon_Lab/EZbakR_paper/Data/bakR_Simulations/"


# Path to grandR + bakR + EZbakR output processed as detailed below
all_processed_output <- "C:/Users/isaac/Documents/Simon_Lab/EZbakR_paper/Data/Figure_5/Two_reps_final_df_rerun_2024_10_02.csv"


# Folder to save figures to
savedir <- "C:/Users/isaac/Documents/Simon_Lab/EZbakR_paper/Figures/Figure_5/"


### File paths for panel F

# Path to Courvan et al. 2024 reanalysis
courvan_ezbdo_path <- "C:/Users/isaac/Documents/Simon_Lab/EZbakR_paper/Data/courvan_data/EZbakRFit_courvan_nocB.rds"


### Functions

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


# Analyze bakR simulated data with bakR and EZbakR -----------------------------

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


write_csv(final_df,
          EZbakR_and_bakR_results)



# grandR results exploration ---------------------------------------------------

results_df <- tibble()
ROPEs <- rep(0, times = 5)
for(i in 1:5){

  gres <- readRDS(paste0(simdata_dir, "grandR_Fit_Rep", i, ".rds"))

  gres_df <- GetAnalysisTable(gres, "Regulation.Condition1 vs Condition2")

  gres_df %>% as_tibble() %>%
    dplyr::select(`Regulation.Condition1 vs Condition2.s.ROPE`)

  colnames(gres_df %>% as_tibble())

  sig_genes <- paste0("Gene", 9001:10000)
  nsig_genes <- paste0("Gene", 1:9000)


  ### ROPE oracle strategy
  cutoffs <- seq(from = 0.5, to = 0.99, length.out = 1000)

  for(c in cutoffs){


    gres_df <- gres_df %>% as_tibble() %>%
      mutate(
        ROPE_conclusion = ifelse(abs(`Regulation.Condition1 vs Condition2.HL.ROPE`) > c,
                                 "Significant",
                                 "Not Sig.")
      )


    FDR <- gres_df %>%
      summarise(
        FDR_ROPE = sum(ROPE_conclusion == "Significant" & Gene %in% nsig_genes)/sum(ROPE_conclusion == "Significant"),
        Power_ROPE = sum(ROPE_conclusion == "Significant" & Gene %in% sig_genes)/length(sig_genes)
      ) %>%
      dplyr::select(FDR_ROPE) %>% unlist() %>% unname()

    if(FDR <= 0.041){
      break
    }


  }


  ROPE_c <- c
  ROPEs[i] <- c


  gres_df <- gres_df %>% as_tibble() %>%
    mutate(
      ROPE_conclusion = ifelse(abs(`Regulation.Condition1 vs Condition2.HL.ROPE`) > 0.95,
                               "Significant",
                               "Not Sig."),
      ROPE_oracle = ifelse(abs(`Regulation.Condition1 vs Condition2.HL.ROPE`) > ROPE_c,
                           "Significant",
                           "Not Sig."),
      CI_conclusion = ifelse(sign(`Regulation.Condition1 vs Condition2.HL.cred.lower`) == sign(`Regulation.Condition1 vs Condition2.HL.cred.upper`),
                             "Significant",
                             "Not Sig."),
    )

  results <- gres_df %>%
    summarise(
      FDR_grandR = sum(ROPE_conclusion == "Significant" & Gene %in% nsig_genes)/sum(ROPE_conclusion == "Significant"),
      Power_grandR = sum(ROPE_conclusion == "Significant" & Gene %in% sig_genes)/length(sig_genes),
      FDR_grandR_o = sum(ROPE_oracle == "Significant" & Gene %in% nsig_genes)/sum(ROPE_oracle == "Significant"),
      Power_grandR_o = sum(ROPE_oracle == "Significant" & Gene %in% sig_genes)/length(sig_genes),
      FDR_grandR_CI = sum(CI_conclusion == "Significant" & Gene %in% nsig_genes)/sum(CI_conclusion == "Significant"),
      Power_grandR_CI = sum(CI_conclusion == "Significant" & Gene %in% sig_genes)/length(sig_genes)
    ) %>%
    mutate(
      MCC_grandR = calc_MCC(length(nsig_genes), length(sig_genes), FDR_grandR, Power_grandR),
      MCC_grandR_o = calc_MCC(length(nsig_genes), length(sig_genes), FDR_grandR_o, Power_grandR_o),
      MCC_grandR_CI = calc_MCC(length(nsig_genes), length(sig_genes), FDR_grandR_CI, Power_grandR_CI),
      Replicate = i
    )

  results_df <- bind_rows(results,
                          results_df)



}


final_df <- read_csv(EZbakR_and_bakR_results) %>%
  filter(cutoff > 0.0499 & cutoff < 0.04999) %>%
  dplyr::select(-cutoff) %>%
  inner_join(results_df,
             by = "Replicate")
write_csv(final_df, all_processed_output)




# Reproduce figures ------------------------------------------------------------

final_df <- read_csv(all_processed_output)


MCCs <- final_df %>%
  dplyr::select(Replicate, starts_with("MCC")) %>%
  tidyr::pivot_longer(cols = starts_with("MCC"),
                      names_to = "model",
                      values_to = "MCC",
                      names_pattern = "MCC_(.*)") %>%
  filter(model != "grandR") %>%
  mutate(model = case_when(
    model == "ezbakR" ~ "EZbakR",
    model == "mle" ~ "bakR (MLE)",
    model == "mcmc" ~ "bakR (MCMC)",
    model == "grandR_CI" ~ "grandR (CI)",
    model == "grandR_o" ~ "grandR (oracle)"
  )) %>%
  mutate(
    model = factor(model,
                   levels = c("bakR (MCMC)",
                              "bakR (MLE)",
                              "grandR (CI)",
                              "grandR (oracle)",
                              "EZbakR"))
  )

avg_MCCs <- final_df %>%
  dplyr::select(Replicate, starts_with("MCC")) %>%
  tidyr::pivot_longer(cols = starts_with("MCC"),
                      names_to = "model",
                      values_to = "MCC",
                      names_pattern = "MCC_(.*)") %>%
  filter(model != "grandR") %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(avg_MCC = mean(MCC),
                   sd_MCC = sd(MCC)) %>%
  mutate(model = case_when(
    model == "ezbakR" ~ "EZbakR",
    model == "mle" ~ "bakR (MLE)",
    model == "mcmc" ~ "bakR (MCMC)",
    model == "grandR_CI" ~ "grandR (CI)",
    model == "grandR_o" ~ "grandR (oracle)"
  )) %>%
  mutate(
    model = factor(model,
                   levels = c("bakR (MCMC)",
                              "bakR (MLE)",
                              "grandR (CI)",
                              "grandR (oracle)",
                              "EZbakR"))
  )



set.seed(44)
gMCC <- ggplot(MCCs,
               aes(x = model, y = MCC, color = model)) +
  geom_errorbar(data = avg_MCCs, mapping = aes(x = model, ymin = avg_MCC - 2*sd_MCC/sqrt(5), ymax = avg_MCC + 2*sd_MCC/sqrt(5)),
                width = 0.4, inherit.aes = FALSE,
                linewidth = 0.2) +
  geom_point(data = avg_MCCs, mapping = aes(x = model, y = avg_MCC), color = "black", inherit.aes = FALSE,
             size = 1)  +
  geom_jitter( height = 0, width = 0.4, size = 0.5,
               color = "gray50") +  theme_classic() +
  xlab("Model") +
  ylab("MCC") +
  ylim(c(0.2, 0.5))  +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10), #change font size of legend title
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position = "none")

gMCC

### Make Power plot as an MCC-like scatter plot

Powers <- final_df %>%
  dplyr::select(Replicate, starts_with("Power")) %>%
  tidyr::pivot_longer(cols = starts_with("Power"),
                      names_to = "model",
                      values_to = "Power",
                      names_pattern = "Power_(.*)") %>%
  filter(model != "grandR") %>%
  mutate(model = case_when(
    model == "ezbakR" ~ "EZbakR",
    model == "mle" ~ "bakR (MLE)",
    model == "mcmc" ~ "bakR (MCMC)",
    model == "grandR_CI" ~ "grandR (CI)",
    model == "grandR_o" ~ "grandR (oracle)"
  )) %>%
  mutate(
    model = factor(model,
                   levels = c("bakR (MCMC)",
                              "bakR (MLE)",
                              "grandR (CI)",
                              "grandR (oracle)",
                              "EZbakR"))
  )

avg_Powers <- final_df %>%
  dplyr::select(Replicate, starts_with("Power")) %>%
  tidyr::pivot_longer(cols = starts_with("Power"),
                      names_to = "model",
                      values_to = "Power",
                      names_pattern = "Power_(.*)") %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(avg_Power = mean(Power),
                   sd_Power = sd(Power)) %>%
  filter(model != "grandR") %>%
  mutate(model = case_when(
    model == "ezbakR" ~ "EZbakR",
    model == "mle" ~ "bakR (MLE)",
    model == "mcmc" ~ "bakR (MCMC)",
    model == "grandR_CI" ~ "grandR (CI)",
    model == "grandR_o" ~ "grandR (oracle)"
  )) %>%
  mutate(
    model = factor(model,
                   levels = c("bakR (MCMC)",
                              "bakR (MLE)",
                              "grandR (CI)",
                              "grandR (oracle)",
                              "EZbakR"))
  )


set.seed(44)
gPower_s <- ggplot(Powers,
                   aes(x = model, y = Power)) +
  geom_errorbar(data = avg_Powers, mapping = aes(x = model, ymin = avg_Power - 2*sd_Power/sqrt(5),
                                                 ymax = avg_Power + 2*sd_Power/sqrt(5)),
                width = 0.4, inherit.aes = FALSE,
                linewidth = 0.2) +
  geom_point(data = avg_Powers, mapping = aes(x = model, y = avg_Power), color = "black", inherit.aes = FALSE,
             size = 1)  +
  geom_jitter( height = 0, width = 0.4, size = 0.5,
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
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  ylim(c(0, 0.3))

gPower_s


### Make FDR plot as an MCC-like scatter plot


FDRs <- final_df %>%
  dplyr::select(Replicate, starts_with("FDR")) %>%
  tidyr::pivot_longer(cols = starts_with("FDR"),
                      names_to = "model",
                      values_to = "FDR",
                      names_pattern = "FDR_(.*)") %>%
  filter(model != "grandR") %>%
  mutate(model = case_when(
    model == "ezbakR" ~ "EZbakR",
    model == "mle" ~ "bakR (MLE)",
    model == "mcmc" ~ "bakR (MCMC)",
    model == "grandR_CI" ~ "grandR (CI)",
    model == "grandR_o" ~ "grandR (oracle)"
  )) %>%
  mutate(
    model = factor(model,
                   levels = c("bakR (MCMC)",
                              "bakR (MLE)",
                              "grandR (CI)",
                              "grandR (oracle)",
                              "EZbakR"))
  )


avg_FDRs <- final_df %>%
  dplyr::select(Replicate, starts_with("FDR")) %>%
  tidyr::pivot_longer(cols = starts_with("FDR"),
                      names_to = "model",
                      values_to = "FDR",
                      names_pattern = "FDR_(.*)") %>%
  filter(model != "grandR") %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(avg_FDR = mean(FDR),
                   sd_FDR = sd(FDR)) %>%
  mutate(model = case_when(
    model == "ezbakR" ~ "EZbakR",
    model == "mle" ~ "bakR (MLE)",
    model == "mcmc" ~ "bakR (MCMC)",
    model == "grandR_CI" ~ "grandR (CI)",
    model == "grandR_o" ~ "grandR (oracle)"
  )) %>%
  mutate(
    model = factor(model,
                   levels = c("bakR (MCMC)",
                              "bakR (MLE)",
                              "grandR (CI)",
                              "grandR (oracle)",
                              "EZbakR"))
  )


set.seed(44)
gFDR_s <- ggplot(FDRs,
                 aes(x = model, y = FDR)) +
  geom_errorbar(data = avg_FDRs, mapping = aes(x = model, ymin = avg_FDR - 2*sd_FDR/sqrt(5),
                                               ymax = avg_FDR + 2*sd_FDR/sqrt(5)),
                width = 0.4, inherit.aes = FALSE,
                linewidth = 0.2) +
  geom_point(data = avg_FDRs, mapping = aes(x = model, y = avg_FDR), color = "black", inherit.aes = FALSE,
             size = 1)  +
  geom_jitter( height = 0, width = 0.4, size = 0.5,
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
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  ylim(c(0, 0.06)) +
  geom_hline(yintercept = 0.05,
             color = 'darkred',
             linetype = 'dotted',
             linewidth = 0.5)


### Save plots

setwd(savedir)
ggsave(plot = gMCC,
       filename = "MCC_2reps_grandR.pdf",
       width = 2.15,
       height = 2.35)
ggsave(plot = gPower_s,
       filename = "Power_scatter_2reps_grandR.pdf",
       width = 2.15,
       height = 2.35)
ggsave(plot = gFDR_s,
       filename = "FDR_scatter_2reps_grandR.pdf",
       width = 2.15,
       height = 2.35)

ggsave(plot = gMCC,
       filename = "MCC_2reps_grandR.png",
       width = 2.15,
       height = 2.35)
ggsave(plot = gPower_s,
       filename = "Power_scatter_2reps_grandR.png",
       width = 2.15,
       height = 2.35)
ggsave(plot = gFDR_s,
       filename = "FDR_scatter_2reps_grandR.png",
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

# MCMC times come from runs on McCleary
# So do grandR times
final_time <- final_time %>%
  mutate(bakR_MCMC_time = c(105720,
                            105420,
                            89880,
                            86580,
                            83460),
         grandR_time = c(18360,
                         18300,
                         18360, # Same within rounding to rep 1
                         17880,
                         17940
         ),
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
    model == "bakR_MCMC" ~ "bakR (MCMC)",
    model == "grandR" ~ "grandR"
  )) %>%
  mutate(
    model = factor(model,
                   levels = c(
                     "bakR (MCMC)",
                     "grandR",
                     "bakR (MLE)",
                     "EZbakR"
                   ))
  )


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
    model == "bakR_MCMC" ~ "bakR (MCMC)",
    model == "grandR" ~ "grandR"
  )) %>%
  mutate(
    model = factor(model,
                   levels = c(
                     "bakR (MCMC)",
                     "grandR",
                     "bakR (MLE)",
                     "EZbakR"
                   ))
  )



set.seed(4)
gTime <- ggplot(Times,
                aes(x = model, y = Runtime)) +
  geom_errorbar(data = avg_Times, mapping = aes(x = model, ymin = avg_runtime - 2*sd_runtime/sqrt(5),
                                                ymax = avg_runtime + 2*sd_runtime/sqrt(5)),
                width = 0.4, inherit.aes = FALSE,
                linewidth = 0.2) +
  geom_point(data = avg_Times, mapping = aes(x = model, y = avg_runtime), color = "black", inherit.aes = FALSE,
             size = 1)  +
  geom_jitter( height = 0, width = 0.4, size = 0.5) +
  theme_classic() +
  xlab("Model") +
  ylab("Runtime (s)") +
  theme(#text=element_text(size=20), #change font size of all text
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    #plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10), #change font size of legend title
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  scale_y_log10(breaks = scales::log_breaks(n = 5))


gTime


set.seed(421)
main_timeplot <- ggplot(Times,
                        aes(x = model, y = Runtime)) +
  geom_errorbar(data = avg_Times, mapping = aes(x = model, ymin = avg_runtime - 2*sd_runtime/sqrt(5),
                                                ymax = avg_runtime + 2*sd_runtime/sqrt(5)),
                width = 0.4, inherit.aes = FALSE,
                linewidth = 0.2) +
  geom_point(data = avg_Times, mapping = aes(x = model, y = avg_runtime), color = "black", inherit.aes = FALSE,
             size = 1)  +
  geom_jitter( height = 0, width = 0.4, size = 0.5,
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
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position = "none")

set.seed(421)
inset_timeplot <- ggplot(Times %>%
                           filter(!(model %in% c("bakR (MCMC)", "grandR"))) %>%
                           mutate(model = factor(model,
                                                 levels = c("bakR (MLE)", "EZbakR"))),
                         aes(x = model, y = Runtime, color = model)) +
  geom_errorbar(data = avg_Times%>%
                  filter(!(model %in% c("bakR (MCMC)", "grandR"))) %>%
                  mutate(model = factor(model,
                                        levels = c("bakR (MLE)", "EZbakR"))),
                mapping = aes(x = model,
                              ymin = avg_runtime - 2*sd_runtime/sqrt(5),
                              ymax = avg_runtime + 2*sd_runtime/sqrt(5)),
                width = 0.4, inherit.aes = FALSE,
                linewidth = 0.2) +
  geom_point(data = avg_Times%>%
               filter(!(model %in% c("bakR (MCMC)", "grandR"))) %>%
               mutate(model = factor(model,
                                     levels = c("bakR (MLE)", "EZbakR"))),
             mapping = aes(x = model,
                           y = avg_runtime),
             color = "black",
             inherit.aes = FALSE,
             size = 1)  +
  geom_jitter( height = 0, width = 0.4, size = 0.5,
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
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position = "none")





setwd(savedir)
ggsave(filename = "Comparing_runtimes_full_grandR.pdf",
       plot = main_timeplot,
       width = 2.25,
       height = 2.65)
ggsave(filename = "Comparing_runtimes_zoom_grandR.pdf",
       plot = inset_timeplot,
       width = 1.5,
       height = 2.08)



# Reanalyze Courvan et al. 2024 data -------------------------------------------


##### Reanalyze data #####

# Download data from GEO
cB1 <- fread("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE245750&format=file&file=GSE245750%5FcB%5FTimeLapse%5FS01andS19%2Ecsv%2Egz")
cB2 <- fread("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE245750&format=file&file=GSE245750%5FcB%5FTimeLapse%5FS10toS18%2Ecsv%2Egz")
cB3 <- fread("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE245750&format=file&file=GSE245750%5FcB%5FTimeLapse%5FS20toS27%2Ecsv%2Egz")

cB <- bind_rows(list(cB1, cB2, cB3))


cB %>% group_by(sample) %>% summarise(mutrate = sum(TC*n)/sum(nT*n), reads = sum(n)) %>%
  arrange()


metadf <- tibble(
  sample = c("S1", paste0("S", 10:27)),
  tl = c(2, 2, 2, 2,
         2, 2, 2,
         2, 2, 2,
         2, 2, 2, 2,
         2, 2, 2,
         0, 0),
  state = factor(c("resting", "resting", "resting", "resting",
                   "activated", "activated", "activated",
                   "activated", "activated", "activated",
                   "resting", "resting", "resting",
                   "resting", "resting", "resting",
                   "resting", "resting", "resting"),
                 levels = c("activated", "resting")),
  oxygen = factor(c("normal", "normal", "normal", "normal",
                    "normal", "normal", "normal",
                    "low", "low", "low",
                    "normal", "normal", "normal", "normal",
                    "low", "low", "low",
                    "normal", "low"),
                  levels = c("normal", "low"))
) %>%
  mutate(
    tl = ifelse(
      sample == "S19",
      0, tl
    )
  )


ezbdo <- EZbakRData(cB %>% filter(TC/nT < 0.25 & nT > 5),
                    metadf)

ezbdo <- EstimateFractions(ezbdo,
                           features = "XF",
                           pold_from_nolabel = TRUE)

ezbdo <- EstimateKinetics(ezbdo)

ezbdo <- AverageAndRegularize(ezbdo,
                              parameter = "log_kdeg",
                              formula_mean = ~state*oxygen)

ezbdo <- AverageAndRegularize(ezbdo,
                              parameter = "log_ksyn",
                              formula_mean = ~state*oxygen)


### Resting + normoxia vs. activated + normoxia
ezbdo$metadata$averages$logkdeg_XF$fit_params

ezbdo <- CompareParameters(
  ezbdo,
  param_name = "stateresting:oxygenlow"
)

ezbdo <- CompareParameters(
  ezbdo,
  parameter = "log_ksyn",
  param_name = "stateresting:oxygenlow"
)



### Resting + normoxia vs. resting + hypoxia


ezbdo <- CompareParameters(
  ezbdo,
  design_factor = "",
  reference = "stateresting:oxygennormal",
  experimental = "stateresting:oxygenlow"
)

ezbdo <- CompareParameters(
  ezbdo,
  parameter = "log_ksyn",
  design_factor = "",
  reference = "stateresting:oxygennormal",
  experimental = "stateresting:oxygenlow"
)



### Activated + normoxia vs. activated + hypoxia

ezbdo <- CompareParameters(
  ezbdo,
  design_factor = "",
  reference = "stateactivated:oxygennormal",
  experimental = "stateactivated:oxygenlow"
)

ezbdo <- CompareParameters(
  ezbdo,
  parameter = "log_ksyn",
  design_factor = "",
  reference = "stateactivated:oxygennormal",
  experimental = "stateactivated:oxygenlow"
)



ezbdo$cB <- NULL
saveRDS(
  ezbdo,
  courvan_ezbdo_path
)


##### Reproduce figures ######

ezbdo <- readRDS(courvan_ezbdo_path)



kdeg1 <- EZget(
  ezbdo,
  type = "comparisons",
  parameter = "log_kdeg",
  design_factor = "",
  reference = "stateresting:oxygennormal",
  experimental = "stateactivated:oxygennormal"
)


kd1_volc <- kdeg1 %>%
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
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_manual(
    values = c("darkcyan", "darkorange", "darkgray")
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = -1,
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = 1,
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.3
  ) +
  xlab("L2FC(kdeg)") +
  ylab("-log10(padj)") +
  theme(legend.position="none") +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10))




kdeg2 <- EZget(
  ezbdo,
  type = "comparisons",
  parameter = "log_kdeg",
  design_factor = "",
  reference = "stateresting:oxygennormal",
  experimental = "stateresting:oxygenlow"
)


kd2_volc <- kdeg2 %>%
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
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_manual(
    values = c("darkcyan", "darkorange", "darkgray")
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = -1,
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = 1,
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.3
  ) +
  xlab("L2FC(kdeg)") +
  ylab("-log10(padj)") +
  theme(legend.position="none") +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10))



kdeg3 <- EZget(
  ezbdo,
  type = "comparisons",
  parameter = "log_kdeg",
  design_factor = "",
  reference = "stateactivated:oxygennormal",
  experimental = "stateactivated:oxygenlow"
)


kd3_volc <- kdeg3 %>%
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
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_manual(
    values = c("darkorange", "darkgray")
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = -1,
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = 1,
    color = "darkred",
    linetype = "dotted",
    linewidth = 0.3
  ) +
  xlab("L2FC(kdeg)") +
  ylab("-log10(padj)") +
  theme(legend.position="none") +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10))


setwd(savedir)
ggsave(
  filename = "RN_vs_AN.pdf",
  plot = kd1_volc,
  width = 2,
  height = 1.6
)
ggsave(
  filename = "RN_vs_RH.pdf",
  plot = kd2_volc,
  width = 2,
  height = 1.6
)
ggsave(
  filename = "AN_vs_AH.pdf",
  plot = kd3_volc,
  width = 2,
  height = 1.6
)


