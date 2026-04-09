rm(list=ls())
library(tidyverse)
library(magrittr)
library(scoringRules)
library(spiox)

g <- glue::glue


load("jura/spiox_jura_fit.RData")
load("jura/spioxmcmc_jura_fit.RData")
load("jura/GpGpm_jura_preds.RData")
load("jura/lmc_jura_fit.RData")

library(gstat)
data(jura)
jura.all <- rbind(jura.pred, jura.val)
n_all <- nrow(jura.all)
coords <- jura.all %>% dplyr::select(Xloc, Yloc) %>%
  apply(2, \(x) (x-min(x))/(max(x)-min(x))) %>% as.matrix()
Y_all <- jura.all %>% dplyr::select(Cd, Co, Cr, Cu, Ni, Pb, Zn) %>% as.matrix() %>% log()
Y_all <- Y_all %>% apply(2, \(y) y - mean(y))
q <- ncol(Y_all)

set.seed(2026)
n_na <- 200
all_cells <- expand.grid(row=1:n_all, col=1:q)
na_idx <- all_cells[sample(nrow(all_cells), n_na), ] %>% as.matrix()
Y_m <- Y_all; Y_m[na_idx] <- NA

compute_metrics <- function(samps, truth, q, idx_col){
  purrr::map_dfr(1:q, function(j){
    wj <- idx_col == j
    if(sum(wj)==0) return(NULL)
    samps_j <- samps[wj, , drop=FALSE]
    truth_j <- truth[wj]
    low <- apply(samps_j, 1, quantile, 0.025)
    upp <- apply(samps_j, 1, quantile, 0.975)
    tibble(outcome=j,
           rmspe=sqrt(mean((rowMeans(samps_j) - truth_j)^2)),
           coverage=mean(truth_j >= low & truth_j <= upp),
           crps=mean(crps_sample(truth_j, samps_j)))
  })
}

compute_metrics <- function(samps, truth, q, idx_col){
  per_outcome <- purrr::map_dfr(1:q, function(j){
    wj <- idx_col == j
    if(sum(wj)==0) return(NULL)
    samps_j <- samps[wj, , drop=FALSE]
    truth_j <- truth[wj]
    low <- apply(samps_j, 1, quantile, 0.025)
    upp <- apply(samps_j, 1, quantile, 0.975)
    tibble(outcome=j,
           rmspe=sqrt(mean((rowMeans(samps_j) - truth_j)^2)),
           coverage=mean(truth_j >= low & truth_j <= upp),
           crps=mean(crps_sample(truth_j, samps_j)))
  })
  pmean <- rowMeans(samps)
  low_all <- apply(samps, 1, quantile, 0.025)
  upp_all <- apply(samps, 1, quantile, 0.975)
  avg <- tibble(outcome=0L,
                rmspe=sqrt(mean((pmean - truth)^2)),
                coverage=mean(truth >= low_all & truth <= upp_all),
                crps=mean(crps_sample(truth, samps)))
  bind_rows(per_outcome, avg)
}

################################################################################
# IOX fix theta
################################################################################
iox_fix_m <- compute_metrics(
  iox_fix_mcmc$Y_missing_samples, Y_all[is.na(Y_m)], q, iox_fix_mcmc$Y_missing_col
) %>% mutate(model="IOXfix")

################################################################################
# IOX full MCMC
################################################################################
iox_full_m <- compute_metrics(
  iox_all_mcmc$Y_missing_samples, Y_all[is.na(Y_m)], q, iox_all_mcmc$Y_missing_col
) %>% mutate(model="IOXfull")

################################################################################
# Parsimonious Matern
################################################################################
pmatern_m <- compute_metrics(
  pmatern_samps, Y_all[na_idx], q, na_idx[,"col"]
) %>% mutate(model="pMatern")

################################################################################
# LMC
################################################################################
meshed_order <- lmc_meshed$savedata$coords_blocking$ix
yhat_lmc <- lmc_meshed$yhat_mcmc %>% abind::abind(along=3) %>% `[`(order(meshed_order),,)

lmc_samps <- do.call(rbind, lapply(1:q, function(j) yhat_lmc[is.na(Y_m[,j]), j, ]))
lmc_truth <- unlist(lapply(1:q, function(j) Y_all[is.na(Y_m[,j]), j]))
lmc_col   <- rep(1:q, sapply(1:q, function(j) sum(is.na(Y_m[,j]))))

lmc_m <- compute_metrics(lmc_samps, lmc_truth, q, lmc_col) %>% mutate(model="LMC")

################################################################################
# Combine and summarize
################################################################################
all_m <- bind_rows(iox_full_m, pmatern_m, lmc_m)

outcome_names <- c("Overall", "Cd","Co","Cr","Cu","Ni","Pb","Zn")
all_m$outcome_name <- outcome_names[all_m$outcome+1]

times <- tibble(
  model = c("IOXfix","IOXfull","pMatern","LMC"),
  time  = c(iox_fix_mcmc$time_elapsed, iox_all_mcmc$time_elapsed,
            pmatern_time_elapsed, lmc_meshed$time_elapsed)
)

all_m <- all_m %>% left_join(times, by="model")

print(all_m %>% arrange(outcome, model), n=Inf)

save(all_m, file="jura/jura_pred_results.RData")

################################################################################
# Table for LaTeX
################################################################################
tab_wide <- all_m %>% dplyr::filter(model %in% c("IOXfull", "pMatern", "LMC")) %>% dplyr::select(-coverage) %>%
  dplyr::select(outcome_name, model, rmspe, crps) %>%
  pivot_longer(c(rmspe, crps), names_to="metric") %>%
  pivot_wider(names_from=model, values_from=value) %>%
  pivot_wider(names_from = metric, values_from = c(IOXfull, pMatern, LMC),
              names_glue = "{metric}_{.value}") %>%
  select(outcome_name, starts_with("rmspe"), starts_with("crps")) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  kbl(format = "latex", booktabs = TRUE, escape = FALSE,
      col.names = c("Outcome", rep(c("IOX", "pMatern", "LMC"), 2))) %>%
  add_header_above(c(" " = 1, "RMSPE" = 3, "CRPS" = 3))

print(tab_wide, n=Inf)
