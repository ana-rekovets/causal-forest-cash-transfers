# ==============================================================================
# Causal Forest Analysis: Mothers' Pension Program and Child Longevity
# ==============================================================================
# Description:
#   This script estimates heterogeneous treatment effects of the Mothers'
#   Pension (MP) program on child longevity using Generalized Random Forests
#   (GRF) as developed by Athey, Tibshirani, and Wager (2019) and Wager and 
#   Athey (2018). The analysis replicates and extends the baseline estimates 
#   from Table 4 (Panel A) of the original study by Aizer et al. (2016). 
#   It introduces causal machine learning to explore treatment heterogeneity, 
#   incorporating robustness checks via alternative covariate sets, propensity 
#   score trimming, and orthogonalization.
#
# Data:
#   - MP_1940.dta      : 1940 Census linkage
#   - MP_controls.dta  : Pre-treatment controls
#   - MP_data.dta      : Main analysis dataset
#   - MP_Ohio.dta      : Ohio sub-sample
#   - MP_WII.dta       : WWII linkage
#
# Output:
#   - ATE estimates across covariate specifications (log points + years)
#   - Propensity score overlap plots
#   - Heterogeneity tests by county and individual characteristics
#   - Comparison: orthogonalized vs. non-orthogonalized forest
#
# Requirements:
#   R >= 4.0, grf >= 0.10.2
#
# Author:  Anastasiia Rekovets
# Date:    July 2025
# License: MIT
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------

# Install packages if not already installed
packages <- c("haven", "readr", "dplyr", "grf", "Hmisc", "xfun",
               "labelled", "sandwich", "lmtest", "ggplot2", "xtable")

installed <- packages %in% rownames(installed.packages())
#if (any(!installed)) install.packages(packages[!installed])

library(haven)
library(dplyr)
library(grf)
library(Hmisc)
library(ggplot2)
library(labelled)
library(xtable)

if (packageVersion("grf") < "0.10.2") {
  stop("This script requires grf >= 0.10.2. Please update the package.")
}

rm(list = ls())
set.seed(1)


# ------------------------------------------------------------------------------
# 1. Load Data
# ------------------------------------------------------------------------------

base_path <- "./data/"

mp_data <- read_dta(paste0(base_path, "MP_data.dta"))
# mp_1940     <- read_dta(paste0(base_path, "MP_1940.dta"))     # Unused
# mp_controls <- read_dta(paste0(base_path, "MP_controls.dta")) # Unused
# mp_ohio     <- read_dta(paste0(base_path, "MP_Ohio.dta"))     # Unused
# mp_wii      <- read_dta(paste0(base_path, "MP_WII.dta"))      # Unused

# Print variable labels for reference if needed
#label_table <- data.frame(
#  Variable    = names(var_label(mp_data)),
#  Description = as.character(var_label(mp_data))
#)
#print(label_table)

# Filter to unique matches (nmatches == 1), as in paper Table 4 Panel A
mp_data_unique <- subset(mp_data, nmatches == 1)


# ------------------------------------------------------------------------------
# 2. Define Covariate Sets
# ------------------------------------------------------------------------------

# Mother characteristics
mom_vars <- c("divorced", "husbandaway", "marst_miss")

# Child characteristics
kid_vars <- c("childageyears", "length_name",
              "sib2", "sib3", "sib4", "sib5", "sib6", "sib7", "sib8",
              "maxage", "minage")

# County-level continuous controls (used in spec. without county FE)
county10_vars <- c("sei_mean", "sei_sd", "p_urban", "p_widows",
                   "children_singlemom", "poverty", "fem_lfp",
                   "child_labor", "val_farm")

# County fixed effects (dummy indicators CID2–CID75)
countyd_vars <- grep("^CID([2-9]|[1-6][0-9]|7[0-5])$",
                     names(mp_data), value = TRUE)

# Matching variable
match_vars <- c("datemiss")

# State fixed effects (dummy indicators S2–S11)
state_vars <- grep("^S([2-9]|10|11)$", names(mp_data), value = TRUE)

# State-year policy variables
state_year_vars <- c("manwrat", "ageent", "labage", "contschl",
                     "gen_total", "char_total", "tot_edu_schools",
                     "work_required", "limited_duration",
                     "monthlyallowfirstchild", "monthlyallowaddchild")

# Birth cohort dummies (BY2–BY26)
cohortd_vars <- grep("^BY([2-9]|1[0-9]|2[0-6])$",
                     names(mp_data), value = TRUE)

# Continuous birth-year controls
cohort_vars <- c("yob1", "yob2")

# Full covariate set used for propensity score model
X_vars_full <- c(kid_vars, mom_vars, match_vars, countyd_vars,
                 state_year_vars, cohortd_vars)


# ------------------------------------------------------------------------------
# 3. Propensity Score Estimation
# ------------------------------------------------------------------------------

# Helper: fit logistic propensity score on a given dataset
fit_pscore <- function(data, formula_vars, outcome = "accepted") {
  f <- as.formula(paste(outcome, "~", paste(formula_vars, collapse = " + ")))
  df <- model.frame(f, data = data)
  model <- glm(f, data = df, family = binomial)
  df$pscore <- predict(model, type = "response")
  df
}

## 3.1 Unique-matches sample
pscore_df <- fit_pscore(mp_data_unique, X_vars_full)
pscore_df$W <- pscore_df$accepted

cat("\n--- Propensity score summary (unique matches) ---\n")
print(summary(pscore_df$pscore))
cat("N =", nrow(pscore_df), "\n")

## 3.2 Full sample
pscore_df2 <- fit_pscore(mp_data, X_vars_full)

cat("\n--- Propensity score summary (full sample) ---\n")
print(summary(pscore_df2$pscore))
cat("N =", nrow(pscore_df2), "\n")

## 3.3 Overlap plot
plot_overlap <- ggplot(pscore_df, aes(x = pscore, fill = factor(W))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("red", "skyblue"),
                    labels = c("Control", "Treated"),
                    name = "Group") +
  labs(title = "Propensity Score Distribution by Treatment Status",
       x = "Propensity Score", y = "Density") +
  theme_minimal(base_size = 14)

print(plot_overlap)


# ------------------------------------------------------------------------------
# 4. Causal Forest Estimation
# ------------------------------------------------------------------------------

# Helper: run a two-step causal forest with variable selection
#   1. Fit raw forest on all X
#   2. Select variables above mean importance
#   3. Re-fit tuned forest on selected variables
run_causal_forest <- function(X, Y, W, clusters,
                               equalize_weights = FALSE) {
  # Step 1: Nuisance forests
  Y.forest <- regression_forest(X, Y,
                                clusters = clusters,
                                equalize.cluster.weights = equalize_weights)
  Y.hat <- predict(Y.forest)$predictions

  W.forest <- regression_forest(X, W,
                                clusters = clusters,
                                equalize.cluster.weights = equalize_weights)
  W.hat <- predict(W.forest)$predictions

  # Step 2: Raw causal forest for variable importance
  cf_raw <- causal_forest(X, Y, W,
                           Y.hat = Y.hat,
                           W.hat = W.hat,
                           clusters = clusters,
                           equalize.cluster.weights = equalize_weights)

  varimp <- variable_importance(cf_raw)
  selected_idx <- which(varimp > mean(varimp))

  # Fallback: use top 10 if nothing is selected
  if (length(selected_idx) == 0) {
    warning("No variables above mean importance; using top 10.")
    selected_idx <- order(varimp, decreasing = TRUE)[1:10]
  }

  # Step 3: Tuned causal forest on selected variables
  cf <- causal_forest(X[, selected_idx], Y, W,
                       Y.hat = Y.hat,
                       W.hat = W.hat,
                       clusters = clusters,
                       equalize.cluster.weights = equalize_weights,
                       tune.parameters = "all")

  list(forest      = cf,
       tau.hat     = predict(cf)$predictions,
       Y.hat       = Y.hat,
       W.hat       = W.hat,
       selected_idx = selected_idx,
       varimp      = varimp)
}

# Helper: print ATE with 95% CI
print_ate <- function(ate, label = "") {
  ci <- round(qnorm(0.975) * ate[2], 3)
  cat(label, "ATE:", round(ate[1], 3), "+/-", ci,
      "  [", round(ate[1] - ci, 3), ";", round(ate[1] + ci, 3), "]\n")
}

# Helper: convert log-scale ATE to additional life-years
ate_to_years <- function(ate_est, ate_se, mean_log_age) {
  mean_age   <- exp(mean_log_age)
  point      <- mean_age * (exp(ate_est) - 1)
  lower      <- mean_age * (exp(ate_est - 1.96 * ate_se) - 1)
  upper      <- mean_age * (exp(ate_est + 1.96 * ate_se) - 1)
  cat(sprintf("  %.2f years (95%% CI: %.2f to %.2f)\n", point, lower, upper))
  invisible(c(point = point, lower = lower, upper = upper))
}


## 4.1 Specification 1: State + cohort FEs only
# (minimal controls; equalize.cluster.weights = FALSE)

cat("\n======== Spec 1: State + Cohort FEs ========\n")

X1 <- model.matrix(~ . - 1,
                   data = mp_data_unique[, c(state_vars, cohortd_vars)])
Y1 <- as.numeric(mp_data_unique$logageatdeath)
W1 <- as.numeric(mp_data_unique$accepted)
cl1 <- mp_data_unique$fips

res1 <- run_causal_forest(X1, Y1, W1, cl1)

ATE1 <- average_treatment_effect(res1$forest)
print_ate(ATE1, "Spec 1 ATE:")

ATE1_ctrl <- average_treatment_effect(res1$forest, target.sample = "control")
print_ate(ATE1_ctrl, "Spec 1 ATE (control):")


## 4.2 Specification 2: All controls, NO county FEs
# (county-level continuous controls instead of dummies)

cat("\n======== Spec 2: All Controls, No County FEs ========\n")

X_vars2 <- c(kid_vars, mom_vars, match_vars, county10_vars,
             state_year_vars, cohortd_vars)

df2 <- na.omit(mp_data_unique[, c(X_vars2, "logageatdeath", "accepted", "fips")])
X2  <- model.matrix(~ . - 1, data = df2[, X_vars2])
Y2  <- as.numeric(df2$logageatdeath)
W2  <- as.numeric(df2$accepted)
cl2 <- df2$fips

res2 <- run_causal_forest(X2, Y2, W2, cl2)

ATE2 <- average_treatment_effect(res2$forest)
print_ate(ATE2, "Spec 2 ATE:")

ATE2_ctrl <- average_treatment_effect(res2$forest, target.sample = "control")
print_ate(ATE2_ctrl, "Spec 2 ATE (control):")


## 4.3 Specification 3: All controls WITH county FEs  [main specification]

cat("\n======== Spec 3: All Controls + County FEs (Main) ========\n")

X_vars3 <- c(kid_vars, mom_vars, match_vars, countyd_vars,
             state_year_vars, cohortd_vars)

df3 <- na.omit(mp_data_unique[, c(X_vars3, "logageatdeath", "accepted", "fips")])
X3  <- model.matrix(~ . - 1, data = df3[, X_vars3])
Y3  <- as.numeric(df3$logageatdeath)
W3  <- as.numeric(df3$accepted)
cl3 <- df3$fips

res3 <- run_causal_forest(X3, Y3, W3, cl3)

ATE3 <- average_treatment_effect(res3$forest)
print_ate(ATE3, "Spec 3 ATE:")

ATE3_ctrl <- average_treatment_effect(res3$forest, target.sample = "control")
print_ate(ATE3_ctrl, "Spec 3 ATE (control):")

cat("Spec 3 effect in years:\n")
ate_to_years(ATE3[1], ATE3[2], mean(Y3))


## 4.4 Specification 4: SSA date-of-birth outcome (logageatdeath_sa)

cat("\n======== Spec 4: SSA Outcome ========\n")

df4 <- na.omit(mp_data_unique[, c(X_vars3, "logageatdeath_sa", "accepted", "fips")])
X4  <- model.matrix(~ . - 1, data = df4[, X_vars3])
Y4  <- as.numeric(df4$logageatdeath_sa)
W4  <- as.numeric(df4$accepted)
cl4 <- df4$fips

res4 <- run_causal_forest(X4, Y4, W4, cl4)

ATE4 <- average_treatment_effect(res4$forest)
print_ate(ATE4, "Spec 4 ATE:")

ATE4_ctrl <- average_treatment_effect(res4$forest, target.sample = "control")
print_ate(ATE4_ctrl, "Spec 4 ATE (control):")


# ------------------------------------------------------------------------------
# 5. Propensity Score Trimming
# ------------------------------------------------------------------------------

# Attach outcome and cluster to pscore data frame
pscore_df$Y    <- mp_data_unique$logageatdeath[as.numeric(rownames(pscore_df))]
pscore_df$fips <- mp_data_unique$fips[as.numeric(rownames(pscore_df))]

cat("\n--- Propensity score quantiles by treatment group ---\n")
cat("Treated:\n");  print(quantile(pscore_df$pscore[pscore_df$W == 1], probs = seq(0, 1, 0.05)))
cat("Control:\n");  print(quantile(pscore_df$pscore[pscore_df$W == 0], probs = seq(0, 1, 0.05)))

# Trim to region of common support (chosen from overlap plot inspection)
lower_trim <- 0.65
upper_trim <- 0.97

pscore_df$keep <- with(pscore_df,
                        pscore >= lower_trim & pscore <= upper_trim)

cat("\nObservations before / after trimming:\n")
print(table(pscore_df$keep))

# Overlap plot: full vs. trimmed
plot_trimmed <- ggplot(pscore_df, aes(x = pscore, fill = factor(W))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ keep, labeller = labeller(keep = c("FALSE" = "Dropped", "TRUE" = "Kept"))) +
  scale_fill_manual(values = c("red", "skyblue"),
                    labels = c("Control", "Treated"), name = "Group") +
  labs(title = "Propensity Score: Trimmed vs. Full Sample",
       x = "Propensity Score", y = "Density") +
  theme_minimal(base_size = 14)

print(plot_trimmed)


# ------------------------------------------------------------------------------
# 6. Causal Forest on Trimmed Sample  [primary heterogeneity analysis]
# ------------------------------------------------------------------------------

df_trim <- subset(pscore_df, keep == TRUE)

X_trim   <- model.matrix(~ . - 1, data = df_trim[, X_vars_full])
Y_trim   <- as.numeric(df_trim$Y)
W_trim   <- as.numeric(df_trim$W)
cl_trim  <- df_trim$fips

res_trim <- run_causal_forest(X_trim, Y_trim, W_trim, cl_trim)

cf_trim      <- res_trim$forest
tau_hat_trim <- res_trim$tau.hat
Y_hat_trim   <- res_trim$Y.hat
W_hat_trim   <- res_trim$W.hat
sel_idx      <- res_trim$selected_idx

ATE_trim <- average_treatment_effect(cf_trim)
print_ate(ATE_trim, "Trimmed sample ATE:")

cat("Trimmed sample effect in years:\n")
ate_to_years(ATE_trim[1], ATE_trim[2], mean(log(exp(pscore_df$Y))))

cat("\nSelected variables:\n")
print(colnames(X_trim)[sel_idx])


# ------------------------------------------------------------------------------
# 7. Heterogeneity Tests
# ------------------------------------------------------------------------------

## 7.1 Calibration test
cat("\n--- Calibration test (trimmed forest) ---\n")
print(test_calibration(cf_trim))

## 7.2 High vs. Low CATE split
high_effect <- tau_hat_trim > median(tau_hat_trim)

ate_high <- average_treatment_effect(cf_trim, subset =  high_effect)
ate_low  <- average_treatment_effect(cf_trim, subset = !high_effect)

diff_ate <- ate_high[1] - ate_low[1]
diff_se  <- sqrt(ate_high[2]^2 + ate_low[2]^2)

cat("\n--- High vs. Low CATE group ATE difference ---\n")
cat(sprintf("Difference: %.3f +/- %.3f\n",
            diff_ate, round(qnorm(0.975) * diff_se, 3)))

## 7.3 Doubly-robust (DR) score per county

# Compute DR score (Nie & Wager, 2021)
dr_score <- tau_hat_trim +
  W_trim / cf_trim$W.hat *
  (Y_trim - cf_trim$Y.hat - (1 - cf_trim$W.hat) * tau_hat_trim) -
  (1 - W_trim) / (1 - cf_trim$W.hat) *
  (Y_trim - cf_trim$Y.hat + cf_trim$W.hat * tau_hat_trim)

# Aggregate DR score to county level
cluster_mat   <- model.matrix(~ factor(cl_trim) - 1)
cluster_size  <- as.numeric(table(cl_trim))
cluster_score <- as.vector(t(cluster_mat) %*% dr_score / cluster_size)

# County-level histogram
hist_county <- ggplot(data.frame(cluster_score = cluster_score),
                      aes(x = cluster_score)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.8) +
  labs(title = "Distribution of Estimated Treatment Effects by County",
       x = "County-Level DR Score", y = "Frequency") +
  theme_minimal(base_size = 14)

print(hist_county)

## 7.4 Formal t-tests for heterogeneity by selected covariate

results_list <- lapply(colnames(X_trim)[sel_idx], function(var_name) {
  cluster_var_mean <- as.vector(t(cluster_mat) %*% X_trim[, var_name]) / cluster_size
  high_group <- cluster_var_mean > median(cluster_var_mean)
  t.test(cluster_score[high_group], cluster_score[!high_group])
})
names(results_list) <- colnames(X_trim)[sel_idx]

# Assemble results table
result_table <- do.call(rbind, lapply(names(results_list), function(var_name) {
  res <- results_list[[var_name]]
  p   <- res$p.value
  data.frame(
    Variable  = var_name,
    Mean_High = round(res$estimate[1], 4),
    Mean_Low  = round(res$estimate[2], 4),
    Diff      = round(res$estimate[1] - res$estimate[2], 4),
    CI_Lower  = round(res$conf.int[1], 4),
    CI_Upper  = round(res$conf.int[2], 4),
    P_Value   = round(p, 4),
    Sig       = ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                ifelse(p < 0.05,  "*",
                ifelse(p < 0.1,   ".", " "))))
  )
}))

result_table <- result_table[order(result_table$P_Value), ]

cat("\n--- Heterogeneity tests by covariate (sorted by p-value) ---\n")
print(result_table, row.names = FALSE)


# ------------------------------------------------------------------------------
# 8. Propensity Score Variation Plots
# ------------------------------------------------------------------------------

## 8.1 By child age
df_what <- data.frame(
  W.hat          = cf_trim$W.hat,
  childageyears  = df_trim$childageyears
)

boxplot(W.hat ~ childageyears, data = df_what,
        ylab = "Propensity Score",
        xlab = "Child's Age (years)",
        main = "Propensity Score by Child Age")

lines(smooth.spline(df_what$childageyears, df_what$W.hat),
      lwd = 2, col = "blue")

## 8.2 By expenditures on schools (quintiles)
pscore_df$tot_edu_schools_q <- cut(
  pscore_df$tot_edu_schools,
  breaks = quantile(pscore_df$tot_edu_schools,
                    probs = seq(0, 1, 0.2), na.rm = TRUE),
  include.lowest = TRUE, dig.lab = 4
)

plot_schools <- ggplot(pscore_df,
                       aes(x = tot_edu_schools_q, y = pscore)) +
  geom_boxplot(fill = "lightgray", color = "black",
               outlier.shape = 1, width = 0.6) +
  stat_summary(fun = median, aes(group = 1),
               geom = "line", color = "blue", linewidth = 1.2) +
  labs(title = "Propensity Score by Expenditures on Schools (quintiles)",
       x = "Expenditures on Schools (quintiles)", y = "Propensity Score") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(plot_schools)

## 8.3 By state money to charities (quintiles)
pscore_df$char_total_q <- cut(
  pscore_df$char_total,
  breaks = quantile(pscore_df$char_total,
                    probs = seq(0, 1, 0.2), na.rm = TRUE),
  include.lowest = TRUE, dig.lab = 4
)

plot_charities <- ggplot(pscore_df,
                         aes(x = char_total_q, y = pscore)) +
  geom_boxplot(fill = "lightgray", color = "black",
               outlier.shape = 1, width = 0.6) +
  stat_summary(fun = median, aes(group = 1),
               geom = "line", color = "blue", linewidth = 1.2) +
  labs(title = "Propensity Score by State Money to Charities (quintiles)",
       x = "State Money to Charities (quintiles)", y = "Propensity Score") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(plot_charities)


# ------------------------------------------------------------------------------
# 9. Sanity Check: County-Level ATE vs. Forest ATE
# ------------------------------------------------------------------------------

ate_county    <- mean(cluster_score)
se_county     <- sqrt(var(cluster_score) / (length(cluster_score) - 1))
ci_low_county <- ate_county - 1.96 * se_county
ci_hi_county  <- ate_county + 1.96 * se_county

mean_log_age <- mean(Y_trim)

to_years <- function(est, mean_log_age) (exp(est) - 1) * exp(mean_log_age)

ATE_forest <- average_treatment_effect(cf_trim)

comparison_table <- data.frame(
  Method  = c("County-level mean (sanity check)", "Causal Forest ATE"),
  ATE_log = c(round(ate_county, 4),    round(ATE_forest[1], 4)),
  CI_log  = c(
    sprintf("[%.3f; %.3f]", ci_low_county, ci_hi_county),
    sprintf("[%.3f; %.3f]",
            ATE_forest[1] - 1.96 * ATE_forest[2],
            ATE_forest[1] + 1.96 * ATE_forest[2])
  ),
  ATE_years = c(
    round(to_years(ate_county,    mean_log_age), 4),
    round(to_years(ATE_forest[1], mean_log_age), 4)
  ),
  CI_years = c(
    sprintf("[%.3f; %.3f]",
            to_years(ci_low_county, mean_log_age),
            to_years(ci_hi_county,  mean_log_age)),
    sprintf("[%.3f; %.3f]",
            to_years(ATE_forest[1] - 1.96 * ATE_forest[2], mean_log_age),
            to_years(ATE_forest[1] + 1.96 * ATE_forest[2], mean_log_age))
  )
)

cat("\n--- ATE comparison table ---\n")
print(comparison_table)

# LaTeX export
print(xtable(comparison_table,
             caption = "ATE Comparison: County-Level Sanity Check vs. Causal Forest",
             label   = "tab:ate_comparison"),
      include.rownames = FALSE)


# ------------------------------------------------------------------------------
# 10. Robustness: No-Cluster and Non-Orthogonalized Forests
# ------------------------------------------------------------------------------

## 10.1 Forest without cluster adjustment
cf_noclust <- causal_forest(
  X_trim[, sel_idx], Y_trim, W_trim,
  Y.hat = Y_hat_trim, W.hat = W_hat_trim,
  tune.parameters = "all"
)

ATE_noclust <- average_treatment_effect(cf_noclust)
print_ate(ATE_noclust, "No-cluster ATE:")

cat("\n--- Calibration test (no clusters) ---\n")
print(test_calibration(cf_noclust))

tau_hat_noclust <- predict(cf_noclust)$predictions

# Treatment effects by county — ordered by median
df_county <- data.frame(county = factor(cl_trim),
                        tau    = tau_hat_noclust)

county_med <- aggregate(tau ~ county, data = df_county, FUN = median)
df_county$county <- factor(df_county$county,
                           levels = county_med$county[order(county_med$tau)])

plot_county <- ggplot(df_county, aes(x = county, y = tau)) +
  geom_boxplot(fill = "lightgray", outlier.shape = 1, width = 0.6) +
  labs(title = "Treatment Effects by County (No Clusters, Ordered)",
       x = "County (ordered by median effect)",
       y = "Estimated Treatment Effect") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))

print(plot_county)

## 10.2 R-loss comparison (clustered vs. no-cluster)
Rloss_clustered <- mean(((Y_trim - Y_hat_trim) -
                           tau_hat_trim * (W_trim - W_hat_trim))^2)
Rloss_noclust   <- mean(((Y_trim - Y_hat_trim) -
                           tau_hat_noclust * (W_trim - W_hat_trim))^2)

cat(sprintf("\nR-loss difference (no-cluster minus clustered): %.6f\n",
            Rloss_noclust - Rloss_clustered))

## 10.3 ANOVA: are county-level effects jointly significant?
cat("\n--- ANOVA: DR score ~ county ---\n")
print(summary(aov(dr_score ~ factor(cl_trim))))

## 10.4 Forest without orthogonalization (W.hat = marginal mean)
cf_noprop <- causal_forest(
  X_trim[, sel_idx], Y_trim, W_trim,
  Y.hat    = Y_hat_trim,
  W.hat    = mean(W_trim),           # no propensity adjustment
  tune.parameters = "all",
  equalize.cluster.weights = TRUE,
  clusters = cl_trim
)

ATE_noprop      <- average_treatment_effect(cf_noprop)
tau_hat_noprop  <- predict(cf_noprop)$predictions

print_ate(ATE_noprop, "Non-orthogonalized ATE:")

plot_compare <- ggplot(
  data.frame(tau_orth   = tau_hat_trim,
             tau_noprop = tau_hat_noprop),
  aes(x = tau_orth, y = tau_noprop)
) +
  geom_point(alpha = 0.4, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1,
              color = "red", linetype = "dashed", linewidth = 1.2) +
  labs(title = "Orthogonalized vs. Non-Orthogonalized CATE Estimates",
       x = "Orthogonalized CATE",
       y = "Non-Orthogonalized CATE") +
  theme_minimal(base_size = 14)

print(plot_compare)

# ==============================================================================
# 11. Export Figures
# ==============================================================================
cat("\n--- Saving figures to ./figures/ ---\n")

# Create the figures directory if it doesn't exist yet
if (!dir.exists("./figures")) {
  dir.create("./figures")
}

# 1. Save Propensity Score Overlap Plot
ggsave("./figures/1_propensity_overlap.png", plot = plot_overlap, 
       width = 8, height = 6, dpi = 300, bg = "white")

# 2. Save Trimmed Sample Overlap Plot
ggsave("./figures/2_trimmed_overlap.png", plot = plot_trimmed, 
       width = 10, height = 6, dpi = 300, bg = "white")

# 3. Save County Treatment Effects Histogram
ggsave("./figures/3_county_effects_hist.png", plot = hist_county, 
       width = 8, height = 6, dpi = 300, bg = "white")

# 4. Save Ordered County Effects Boxplot (No Clusters)
ggsave("./figures/4_county_effects_boxplot.png", plot = plot_county, 
       width = 12, height = 6, dpi = 300, bg = "white")

# 5. Save Orthogonalized vs Non-Orthogonalized Comparison
ggsave("./figures/5_orthogonalization_comparison.png", plot = plot_compare, 
       width = 8, height = 8, dpi = 300, bg = "white")

# 6. Save Base R Boxplot (Child Age vs P-score)
png("./figures/6_pscore_by_age.png", width = 800, height = 600, res = 100)
boxplot(W.hat ~ childageyears, data = df_what,
        ylab = "Propensity Score",
        xlab = "Child's Age (years)",
        main = "Propensity Score by Child Age")
lines(smooth.spline(df_what$childageyears, df_what$W.hat), lwd = 2, col = "blue")
dev.off()

# 7. Save Expenditures on Schools Plot
ggsave("./figures/7_pscore_by_schools.png", plot = plot_schools, 
       width = 10, height = 6, dpi = 300, bg = "white")

# 8. Save State Money to Charities Plot
ggsave("./figures/8_pscore_by_charities.png", plot = plot_charities, 
       width = 10, height = 6, dpi = 300, bg = "white")

cat("All figures successfully saved!\n")
# ==============================================================================
# End of script
# ==============================================================================
