rm(list=ls())
library(tidyverse)
library(magrittr)
library(xtable)
g <- glue::glue

load("sim_pmatern/pred_results.RData")
load("sim_pmatern/estim_results.RData")

q <- 5

fmt <- function(mn, sd) sprintf("%.3f {\\scriptsize(%.3f)}", mn, sd)

##############################################################################
# Table 1: Predictive performance
##############################################################################
pred_tab <- pred_df %>%
  pivot_longer(c(rmspe, coverage, crps), names_to="metric") %>%
  group_by(model, outcome, metric) %>%
  summarise(mn=mean(value), sd=sd(value), .groups="drop") %>%
  mutate(cell = fmt(mn, sd)) %>%
  select(model, outcome, metric, cell) %>%
  pivot_wider(names_from=outcome, values_from=cell, names_prefix="$Y_") %>%
  rename_with(~paste0(.x, "$"), starts_with("$Y_")) %>%
  arrange(metric, model)

cat("\n% Predictive performance\n")
print(xtable(pred_tab, caption="Predictive performance: mean (sd) across 200 replications.",
             label="tab:pred"), 
      sanitize.text.function=identity, include.rownames=FALSE,
      file="sim_pmatern/tab_pred.tex")

##############################################################################
# Table 2: Estimation of S0
##############################################################################
estim_S0 <- estim_df %>%
  filter(type %in% c("S0_diag","S0_offdiag")) %>%
  group_by(rep, type, model) %>%
  summarise(rmse=sqrt(mean(sq_err)),
            cover=mean(cover, na.rm=TRUE),
            crps=mean(crps, na.rm=TRUE), .groups="drop") %>%
  pivot_longer(c(rmse, cover, crps), names_to="metric") %>%
  group_by(type, model, metric) %>%
  summarise(mn=mean(value, na.rm=TRUE), sd=sd(value, na.rm=TRUE), .groups="drop") %>%
  mutate(cell = ifelse(is.nan(mn), "---", fmt(mn, sd))) %>%
  select(type, model, metric, cell) %>%
  pivot_wider(names_from=metric, values_from=cell) %>%
  arrange(type, model)

cat("\n% Estimation of S(0)\n")
print(xtable(estim_S0, caption="Estimation of $\\mathbf{S}(\\mathbf{0})$: mean (sd) across 200 replications.",
             label="tab:S0"),
      sanitize.text.function=identity, include.rownames=FALSE,
      file="sim_pmatern/tab_S0.tex")

##############################################################################
# Table 3: Estimation of Qpc0 (off-diagonal)
##############################################################################
estim_Qpc0 <- estim_df %>%
  filter(type == "Qpc0_offdiag") %>%
  group_by(rep, model) %>%
  summarise(rmse=sqrt(mean(sq_err)),
            cover=mean(cover, na.rm=TRUE),
            crps=mean(crps, na.rm=TRUE), .groups="drop") %>%
  pivot_longer(c(rmse, cover, crps), names_to="metric") %>%
  group_by(model, metric) %>%
  summarise(mn=mean(value, na.rm=TRUE), sd=sd(value, na.rm=TRUE), .groups="drop") %>%
  mutate(cell = ifelse(is.nan(mn), "---", fmt(mn, sd))) %>%
  select(model, metric, cell) %>%
  pivot_wider(names_from=metric, values_from=cell) %>%
  arrange(model)

cat("\n% Estimation of Qpc(0)\n")
print(xtable(estim_Qpc0, caption="Estimation of off-diagonal $\\mathbf{R}(\\mathbf{0})$: mean (sd) across 200 replications.",
             label="tab:Qpc0"),
      sanitize.text.function=identity, include.rownames=FALSE,
      file="sim_pmatern/tab_Qpc0.tex")

##############################################################################
# Boxplots
##############################################################################
model_cols <- c("IOX"="#2166AC", "pMatern"="#B2182B", "LMC"="#7CAE00")
model_labs <- c("IOX"="IOX", "pMatern"="Pars. Mat\u00e9rn", "LMC"="LMC")

theme_box <- theme_minimal(base_size=12) +
  theme(legend.position="bottom",
        strip.text=element_text(face="bold"),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank())

# --- Predictive ---
pred_long <- pred_df %>%
  pivot_longer(c(rmspe, crps), names_to="metric") %>%
  mutate(outcome=factor(outcome, labels=paste0("Y",1:q)),
         model=factor(model, levels=names(model_labs), labels=model_labs),
         metric=factor(metric, levels=c("rmspe","crps"), labels=c("RMSPE","CRPS")))

gg_pred <- ggplot(pred_long, aes(x=outcome, y=value, fill=model)) +
  geom_boxplot(outlier.size=.5, linewidth=.3, alpha=.85,
               position=position_dodge(.75), width=.65) +
  facet_wrap(~metric, scales="free_y", nrow=1) +
  scale_fill_manual(values=setNames(model_cols, model_labs)) +
  scale_y_log10() +
  labs(x="Outcome", y=NULL, fill=NULL) + theme_box

ggsave("sim_pmatern/fig_pred_boxplot.pdf", gg_pred, width=8, height=3.5)

# --- Estimation S0 ---
estim_S0_long <- estim_df %>%
  filter(type %in% c("S0_diag","S0_offdiag")) %>%
  group_by(rep, type, model) %>%
  summarise(rmse=sqrt(mean(sq_err)),
            crps=mean(crps, na.rm=TRUE), .groups="drop") %>%
  pivot_longer(c(rmse, crps), names_to="metric") %>%
  filter(!(metric=="crps" & model=="pMatern")) %>%
  mutate(model=factor(model, levels=names(model_labs), labels=model_labs),
         metric=factor(metric, levels=c("rmse","crps"), labels=c("RMSE","CRPS")))

gg_S0 <- ggplot(estim_S0_long, aes(x=type, y=value, fill=model)) +
  geom_boxplot(outlier.size=.5, linewidth=.3, alpha=.85,
               position=position_dodge(.75), width=.65) +
  facet_wrap(~metric, scales="free_y", nrow=1) +
  scale_fill_manual(values=setNames(model_cols, model_labs)) +
  scale_x_discrete(labels=c("S0_diag"="Diagonal","S0_offdiag"="Off-diagonal")) +
  scale_y_log10() +
  labs(x=NULL, y=NULL, fill=NULL) + theme_box

ggsave("sim_pmatern/fig_S0_boxplot.pdf", gg_S0, width=8, height=3.5)

# --- Estimation Qpc0 ---
estim_Qpc0_long <- estim_df %>%
  filter(type=="Qpc0_offdiag") %>%
  group_by(rep, model) %>%
  summarise(rmse=sqrt(mean(sq_err)),
            crps=mean(crps, na.rm=TRUE), .groups="drop") %>%
  pivot_longer(c(rmse, crps), names_to="metric") %>%
  filter(!(metric=="crps" & model=="pMatern")) %>%
  mutate(model=factor(model, levels=names(model_labs), labels=model_labs),
         metric=factor(metric, levels=c("rmse","crps"), labels=c("RMSE","CRPS")))

gg_Qpc0 <- ggplot(estim_Qpc0_long, aes(x=model, y=value, fill=model)) +
  geom_boxplot(outlier.size=.5, linewidth=.3, alpha=.85, width=.5) +
  facet_wrap(~metric, scales="free_y", nrow=1) +
  scale_fill_manual(values=setNames(model_cols, model_labs)) +
  labs(x=NULL, y=NULL, fill=NULL) + theme_box +
  theme(legend.position="none")

ggsave("sim_pmatern/fig_Qpc0_boxplot.pdf", gg_Qpc0, width=6, height=3.5)