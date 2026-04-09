library(spiox)
library(tidyverse)
library(magrittr)
library(meshed)
g <- glue::glue

args <- commandArgs(trailingOnly = TRUE)
s_start <- as.integer(args[1])
s_end   <- as.integer(args[2])

mcmc <- 5000
for(s in s_start:s_end){
  load(g("sim_pmatern/reps/data_{s}.RData"))
  na_idx <- which(is.na(Y_m), arr.ind = TRUE)
  
  failures <- character(0)
  if(T){
  ################################################################################
  # IOX via spiox
  ################################################################################
    tryCatch({
      set.seed(2026)
      marginal_fit_time <- system.time({
        marginal_fix_Theta <- autostart(Y_m, X[,1,drop=F], coords, m=30, method = "response", nu = 0)
        })["elapsed"]
      set.seed(2026)
      nrow(coords) -> n
      spiox_time_mcmc <- system.time({
        iox_fix_mcmc <- spiox(Y = Y_m, X = X[,1,drop=F], coords = coords, m = 30,
                               method = "response", fit = "mcmc", iter = mcmc,
                               starting = marginal_fix_Theta,
                               opts=list(update_Theta=c(0,0,0), nu=0, num_threads=2),
                               print_every = 100)
      })["elapsed"]
      iox_fix_mcmc$time_elapsed <- spiox_time_mcmc + marginal_fit_time
      save(iox_fix_mcmc, file=g("sim_pmatern/reps/spiox_pmatern_fit5_{s}.RData"))
    }, error = function(e){
      print(conditionMessage(e))
      failures <<- c(failures, paste0("IOX: ", conditionMessage(e)))
    })
  }
  
  if(T){
    ################################################################################
    # IOX full MCMC via spiox
    ################################################################################
    tryCatch({
      set.seed(2026)
      nrow(coords) -> n
      spiox_time_mcmc <- system.time({
        iox_all_mcmc <- spiox(Y = Y_m, X = X[,1,drop=F], coords = coords, m = 30,
                               method = "response", fit = "mcmc", iter = mcmc,
                               starting = list(Theta=rbind(rep(40,q), 1, 1e-2)),
                               opts=list(update_Theta=c(1,1,1), nu=0, num_threads=16),
                               print_every = 100)
      })["elapsed"]
      iox_all_mcmc$time_elapsed <- spiox_time_mcmc
      save(iox_all_mcmc, file=g("sim_pmatern/reps/spioxmcmc_pmatern_fit5_{s}.RData"))
    }, error = function(e){
      print(conditionMessage(e))
      failures <<- c(failures, paste0("IOX full mcmc: ", conditionMessage(e)))
    })
  }
  
  if(F){
  ################################################################################
  # Parsimonious Matern via GpGpm
  ################################################################################
    tryCatch({
      source("fit_parsimonious_matern.R")
      set.seed(2026)
      gpgpm_fit <- parsimonious_matern(Y_m, X, coords, m=30,
                                       filename=g("sim_pmatern/reps/GpGpm_pmatern_fit5_{s}.RData"))
      
      # build prediction inputs for NA locations
      n_test <- nrow(na_idx)
      Xo_list <- list()
      locs_list <- list()
      for(k in 1:n_test){
        i <- na_idx[k, "row"]
        j <- na_idx[k, "col"]
        Xo_list[[k]] <- c(1, j)
        locs_list[[k]] <- c(as.numeric(coords[i, ]), j)
      }
      Xo_pmatern <- do.call(rbind, Xo_list) %>% as.data.frame()
      colnames(Xo_pmatern) <- c("Intercept", "j")
      Xo_pmatern$j <- factor(Xo_pmatern$j)
      Xo_pmatern <- model.matrix(data=Xo_pmatern, ~.-1)[,-(q+1)]
      colnames(Xo_pmatern)[1] <- "Intercept"
      locs_test <- do.call(rbind, locs_list)
      
      set.seed(2026)
      pmatern_sim_time <- system.time({
        pmatern_samps <- GpGpm::cond_sim(
          fit = gpgpm_fit$fit,
          locs_pred = as.matrix(locs_test),
          X_pred = Xo_pmatern,
          nsims = mcmc/4
        )
      })["elapsed"]
      pmatern_time_elapsed <- gpgpm_fit$fit$time_elapsed + pmatern_sim_time
      
      save(pmatern_samps, pmatern_time_elapsed, na_idx,
           file=g("sim_pmatern/reps/GpGpm_pmatern_preds_{s}.RData"))
    }, error = function(e){
      failures <<- c(failures, paste0("pMatern: ", conditionMessage(e)))
    })
  }
  
  if(F){
  ################################################################################
  # LMC via meshed
  ################################################################################
    tryCatch({
      set.seed(2026)
      lmc_time_mcmc <- system.time({
        lmc_meshed <- spmeshed(y=Y_m, x=X,
                               coords=coords, k=q, block_size=40, n_samples=mcmc, n_threads=16, verbose=10,
                               prior=list(phi=c(0.3, 200)))
      })
      lmc_meshed$v_mcmc <- NULL
      lmc_meshed$w_mcmc <- NULL
      lmc_meshed$lp_mcmc <- NULL
      lmc_meshed$time_elapsed <- lmc_time_mcmc["elapsed"]
      save(lmc_meshed, file=g("sim_pmatern/reps/lmc_pmatern_fit5_{s}.RData"))
    }, error = function(e){
      failures <<- c(failures, paste0("LMC: ", conditionMessage(e)))
    })
  }
  ################################################################################
  # Log failures
  ################################################################################
  if(length(failures) > 0){
    writeLines(
      c(paste0("rep: ", s), failures),
      g("sim_pmatern/reps/failures_{s}.txt")
    )
    cat("Failures for rep", s, ":\n", paste(failures, collapse="\n"), "\n")
  } else {
    cat("All methods completed successfully for rep", s, "\n")
  }

}
