rm(list=ls())
library(tidyverse)
library(magrittr)
library(scoringRules)

g <- glue::glue

source("functions.R")

for(s in 1:200){
  
  load(g("sim_pmatern/reps/data_{s}.RData"))
  load(g("sim_pmatern/reps/spiox_pmatern_fit5_{s}.RData"))
  load(g("sim_pmatern/reps/spioxmcmc_pmatern_fit5_{s}.RData"))
  load(g("sim_pmatern/reps/GpGpm_pmatern_preds_{s}.RData"))
  load(g("sim_pmatern/reps/lmc_pmatern_fit5_{s}.RData"))
  
  ################################################################################
  # Helper
  ################################################################################
  compute_metrics <- function(samps, truth, q, idx_col){
    rmspe <- coverage <- crps_val <- numeric(q)
    for(j in 1:q){
      wj <- idx_col == j
      samps_j <- samps[wj, , drop=FALSE]
      truth_j <- truth[wj]
      rmspe[j] <- sqrt(mean((rowMeans(samps_j) - truth_j)^2))
      low_j <- apply(samps_j, 1, quantile, 0.025)
      upp_j <- apply(samps_j, 1, quantile, 0.975)
      coverage[j] <- mean(truth_j >= low_j & truth_j <= upp_j)
      crps_val[j] <- mean(crps_sample(truth_j, samps_j))
    }
    data.frame(outcome=1:q, rmspe, coverage, crps=crps_val)
  }
  
  ################################################################################
  # IOX fix theta
  ################################################################################
  iox_fix_metrics <- compute_metrics(
    iox_fix_mcmc$Y_missing_samples, Y[na_idx], q, na_idx[,"col"]
  )
  
  ################################################################################
  # IOX fix theta
  ################################################################################
  iox_full_metrics <- compute_metrics(
    iox_all_mcmc$Y_missing_samples, Y[na_idx], q, na_idx[,"col"]
  )
  
  ################################################################################
  # Parsimonious Matern
  ################################################################################
  pmatern_metrics <- compute_metrics(
    pmatern_samps, Y[na_idx], q, na_idx[,"col"]
  )
  
  ################################################################################
  # LMC
  ################################################################################
  meshed_order <- lmc_meshed$savedata$coords_blocking$ix
  yhat_lmc <- lmc_meshed$yhat_mcmc %>% abind::abind(along=3) %>% `[`(order(meshed_order),,)
  
  lmc_samps_list <- lapply(1:q, function(j){
    na_j <- is.na(Y_m[, j])
    yhat_lmc[na_j, j, ]
  })
  lmc_samps <- do.call(rbind, lmc_samps_list)
  lmc_truth <- unlist(lapply(1:q, function(j) Y[is.na(Y_m[,j]), j]))
  lmc_col <- rep(1:q, sapply(1:q, function(j) sum(is.na(Y_m[,j]))))
  lmc_metrics <- compute_metrics(lmc_samps, lmc_truth, q, lmc_col)
  
  ################################################################################
  # Summary table
  ################################################################################
  results <- data.frame(
    outcome          = 1:q,
    rmspe_IOX_full        = iox_full_metrics$rmspe,
    rmspe_IOX_fix        = iox_fix_metrics$rmspe,
    rmspe_pMatern    = pmatern_metrics$rmspe,
    rmspe_LMC        = lmc_metrics$rmspe,
    coverage_IOX_full     = iox_full_metrics$coverage,
    coverage_IOX_fix     = iox_fix_metrics$coverage,
    coverage_pMatern = pmatern_metrics$coverage,
    coverage_LMC     = lmc_metrics$coverage,
    crps_IOX_full         = iox_full_metrics$crps,
    crps_IOX_fix         = iox_fix_metrics$crps,
    crps_pMatern     = pmatern_metrics$crps,
    crps_LMC         = lmc_metrics$crps,
    time_IOX_full         = iox_all_mcmc$time_elapsed,
    time_IOX_fix         = iox_fix_mcmc$time_elapsed,
    time_pMatern     = pmatern_time_elapsed,
    time_LMC         = lmc_meshed$time_elapsed
  )
  print(results)
  
  save(results, file=g("sim_pmatern/reps/results_{s}.RData"))
}





################################################################################
# Aggregate across replications
################################################################################
all_results <- list()
for(ss in 1:200){
  f <- g("sim_pmatern/reps/results_{ss}.RData")
  if(file.exists(f)){
    load(f)
    results$rep <- ss
    all_results[[ss]] <- results
  }
}
all_results <- bind_rows(all_results) 

colnames(all_results) <- c('outcome','rmspe_IOXfull','rmspe_IOXfix','rmspe_pMatern','rmspe_LMC','coverage_IOXfull','coverage_IOXfix','coverage_pMatern','coverage_LMC','crps_IOXfull','crps_IOXfix','crps_pMatern','crps_LMC','time_IOXfull','time_IOXfix','time_pMatern','time_LMC','rep')

summary_table <- all_results %>%
  pivot_longer(cols = -c(outcome, rep), names_to = c("metric", "model"),
               names_pattern = "(.+)_(.+)") %>%
  group_by(metric, model, outcome) %>%
  summarise(mean = mean(value), sd = sd(value), .groups="drop") %>%
  pivot_wider(names_from = outcome, values_from = c(mean, sd),
              names_glue = "Y{outcome}_{.value}") %>%
  arrange(metric, model)

print(summary_table, n=Inf, width=Inf)

save(all_results, summary_table,
     file="sim_pmatern/prediction_results.RData")

all_pred_results <- all_results %>% dplyr::select(contains("pMatern"), contains("IOXfull"), contains("LMC"), rep, outcome) %>% head()

## estimating coregionalization matrix and partial correlation matrix at distance zero

make_part_corr <- function(Q){
  temp <- -cov2cor(Q)
  diag(temp) <- 1
  return(temp)
}

outfile <- g("sim_pmatern/estimation_results.RData")

if(file.exists(outfile)){
  load(outfile)
  done <- unique(df$rep)
} else {
  df <- tibble()
  done <- integer(0)
}

for(s in 1:200){
  if(s %in% done){ cat(s, "skip\n"); next }
  cat(s, "\n")
  
  load(g("sim_pmatern/reps/data_{s}.RData"))
  load(g("sim_pmatern/reps/spioxmcmc_pmatern_fit5_{s}.RData"))
  load(g("sim_pmatern/reps/GpGpm_pmatern_fit5_{s}.RData"))
  load(g("sim_pmatern/reps/lmc_pmatern_fit5_{s}.RData"))
  
  S0 <- Sigma * Gamma_nu
  Qpc <- make_part_corr(solve(Sigma))
  Qpc0 <- Qpc * Gamma_nu
  
  mcmc <- dim(iox_all_mcmc$Sigma)[3]
  
  tail_mcmc <- 200
  iox_fij0 <- spiox:::iox_fij0(iox_all_mcmc, n_threads=16, matern=1, tail_mcmc=tail_mcmc)
  
  iox_S0_mcmc <- tail(iox_all_mcmc$Sigma, c(NA,NA,tail_mcmc)) * iox_fij0
  iox_Qpc0_mcmc <- iox_fij0 * (tail(iox_all_mcmc$Sigma, c(NA,NA,tail_mcmc)) %>%
                                 apply(3, \(x) make_part_corr(solve(x))) %>%
                                 array(dim=c(q,q,tail_mcmc)))
  
  pmatern_S0 <- gpgpm_fit$gamma_nu * gpgpm_fit$Sigma
  pmatern_Qpc0 <- gpgpm_fit$gamma_nu * make_part_corr(solve(gpgpm_fit$Sigma))
  
  lmc_S0_mcmc <- lmc_meshed$lambda_mcmc %>% meshed:::cube_tcrossprod()
  
  results_list <- list()
  
  for(i in 1:q){
    for(j in i:q){
      type <- if(i == j) "S0_diag" else "S0_offdiag"
      truth <- S0[i,j]
      
      samp_iox <- iox_S0_mcmc[i,j,]
      ci_iox <- quantile(samp_iox, c(0.025, 0.975))
      results_list[[length(results_list)+1]] <- tibble(
        rep=s, i=i, j=j, type=type, model="IOX",
        truth=truth, estimate=mean(samp_iox),
        sq_err=(mean(samp_iox)-truth)^2,
        cover=as.numeric(truth >= ci_iox[1] & truth <= ci_iox[2]),
        crps=crps_sample(y=truth, dat=samp_iox))
      
      samp_lmc <- lmc_S0_mcmc[i,j,]
      ci_lmc <- quantile(samp_lmc, c(0.025, 0.975))
      results_list[[length(results_list)+1]] <- tibble(
        rep=s, i=i, j=j, type=type, model="LMC",
        truth=truth, estimate=mean(samp_lmc),
        sq_err=(mean(samp_lmc)-truth)^2,
        cover=as.numeric(truth >= ci_lmc[1] & truth <= ci_lmc[2]),
        crps=crps_sample(y=truth, dat=samp_lmc))
      
      results_list[[length(results_list)+1]] <- tibble(
        rep=s, i=i, j=j, type=type, model="pMatern",
        truth=truth, estimate=pmatern_S0[i,j],
        sq_err=(pmatern_S0[i,j]-truth)^2,
        cover=NA_real_, crps=NA_real_)
    }
  }
  
  for(i in 1:q){
    for(j in (i+1):q){
      if(j > q) next
      truth <- Qpc0[i,j]
      
      samp <- iox_Qpc0_mcmc[i,j,]
      ci <- quantile(samp, c(0.025, 0.975))
      results_list[[length(results_list)+1]] <- tibble(
        rep=s, i=i, j=j, type="Qpc0_offdiag", model="IOX",
        truth=truth, estimate=mean(samp),
        sq_err=(mean(samp)-truth)^2,
        cover=as.numeric(truth >= ci[1] & truth <= ci[2]),
        crps=crps_sample(y=truth, dat=samp))
      
      results_list[[length(results_list)+1]] <- tibble(
        rep=s, i=i, j=j, type="Qpc0_offdiag", model="pMatern",
        truth=truth, estimate=pmatern_Qpc0[i,j],
        sq_err=(pmatern_Qpc0[i,j]-truth)^2,
        cover=NA_real_, crps=NA_real_)
    }
  }
  
  df <- bind_rows(df, bind_rows(results_list))
  save(df, file=outfile)
}

# --- Summary table ---
summary_tab <- df %>%
  group_by(type, model) %>%
  summarise(
    RMSE = sqrt(mean(sq_err)),
    Coverage = mean(cover, na.rm=TRUE),
    CRPS = mean(crps, na.rm=TRUE),
    .groups="drop") %>%
  arrange(type, model)
print(summary_tab, n=40)

save(df, summary_tab, file=outfile)




## summary table

estimation_summary <- df %>% group_by(type, model) %>% summarise(rmse = sqrt(mean(sq_err)), crps = mean(crps))

prediction_summary <- all_pred_results %>% dplyr::select(-rep, -outcome) %>%
  mutate(row = row_number()) %>%
  pivot_longer(-row, names_to = c("measure", "model"), names_pattern = "(.+)_(.+)") %>%
  select(-row) %>% 
  group_by(measure, model) %>% 
  summarise_all(mean) %>% 
  dplyr::filter(measure %in% c("rmspe", "crps")) %>%
  pivot_wider(id_cols = model, names_from = measure, values_from = value) %>% 
  mutate(type="Predict") %>%
  rename(rmse = rmspe)
prediction_summary[prediction_summary$model=="IOXfull", "model"] <- "IOX"

library(kableExtra)

df_table_latex <- bind_rows(estimation_summary, prediction_summary) %>%
  pivot_wider(names_from = model, values_from = c(rmse, crps), names_glue = "{.value}_{model}") 

df_table_latex$type <- c("Partial correlations", "Marginal variance",  "Cross-covariance", "Predictions")
df_table_latex %>%
  kbl(format = "latex", booktabs = TRUE, escape = FALSE,
      col.names = c("Type", rep(c("IOX", "Pars. Matérn", "LMC"), 2)),
      digits=3) %>%
  add_header_above(c(" " = 1, "RMSE" = 3, "CRPS" = 3))



