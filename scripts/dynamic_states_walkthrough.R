############################################################
# Dynamic States Walkthrough (Step 1–3)
# Reuses existing repo patterns: imputeTS + nlme residuals +
# doremi::calculate.gold for derivatives (GOLD method).
#
# Inputs: data/ARI PII.3 S1_T0123_v7_noqual.sav (default)
#         columns expected: ID, week, sday, GrL, GrF, e2, e3, e6
#
# Outputs:
#  - A tidy feature file with derivatives per ID x context x window:
#      data/ds_features.csv
#  - Printed summaries and example multilevel models
#
# Notes:
#  - Embedding dimension is the GOLD tau window length (>= 3; default 4)
#  - Imputation can be 'linear' or 'spline' (imputeTS::na_interpolation)
#  - Residuals are level-1 (within-person) from nlme::lme per context
#  - Derivatives are computed within ID x context groups
#  - Operationalization clarification:
#      We operationalize “frequency of fluctuation” via the magnitude of the
#      first derivative (velocity) and “speed of return” via the second derivative
#      (acceleration) using the GOLD method. These are descriptive, within-person
#      features; formal damping ratio and natural frequency parameters require a
#      second-order ODE fit (e.g., x'' = -2*zeta*omega*x' - omega^2*x).
#
# Method selection guide (when to use what):
#  - Use derivative proxies (this script) when:
#      * You have short, sparse, or daily series (e.g., 7 points/week) and want
#        micro-level “ebb/flow” and “return speed” features.
#      * You want a transparent, low-assumption, reproducible pipeline.
#      * Your research question is short-term, within-person dynamics and event links.
#  - Prefer ODE/ctsem/state-space when:
#      * You need interpretable system parameters (omega, zeta) for macro-level or
#        cross-study comparability; you have longer/denser series or multiple windows.
#      * You want to model inputs explicitly and handle serial correlation/noise.
#      * You aim for population-level inference on dynamics (hierarchical ω, ζ).
#  - Population-level inference with this script:
#      * Aggregate local features per ID×context (e.g., mean |velocity|, mean accel).
#      * Fit multilevel models on those summaries (random intercepts/slopes, moderators).
#  - Timescale matters: Derivatives depend on sampling and embedding; report units and
#    sensitivity (e.g., tau = 3–5; linear vs spline interpolation; optional normalization
#    like mean|velocity| / sd(residual)).
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(nlme)
  library(imputeTS)
  library(doremi)     # for calculate.gold (GOLD derivatives)
  library(lmerTest)   # multilevel models for examples
})

## ----------------------
## User parameters
## ----------------------
embedding_dim <- 4         # GOLD window length (tau); use 3–5 per theory/design
impute_option <- "linear"   # 'linear' or 'spline' (see Step 1 rationale)
start_time_at_zero <- TRUE  # start each context at t = 0..6 (keeps step size consistent)
run_sensitivity <- FALSE    # set TRUE to compare embeddings and imputation methods

# Data source selection:
#  - Set data_mode <- "simulate" to generate example data
#  - Or set data_mode <- "csv" and specify csv_path below
data_mode <- "simulate"  # "simulate" or "csv"
csv_path  <- here("data", "identity_timeseries.csv")

## ----------------------
## Data load (defaults to SPSS in data/)
## ----------------------
if (identical(data_mode, "csv") && file.exists(csv_path)) {
  message("Loading CSV data: ", csv_path)
  raw <- readr::read_csv(csv_path, show_col_types = FALSE)
} else {
  message("Using simulated data (set data_mode='csv' to use your own file)")
  sim_src <- here("scripts", "simulate_dynamic_identity.R")
  if (!file.exists(sim_src)) stop("Missing simulation helper: ", sim_src)
  source(sim_src)
  raw <- simulate_identity_data(n_id = 30)
}

## Keep core variables and build event composite
# Validate and prepare columns
req_cols_min <- c("ID", "week", "sday", "GrL", "GrF")
if (!all(req_cols_min %in% names(raw))) {
  stop("Input must contain at least columns: ", paste(req_cols_min, collapse = ", "))
}

# Accept either e2,e3,e6 or a pre-computed event_agg
if (!("event_agg" %in% names(raw))) {
  if (all(c("e2","e3","e6") %in% names(raw))) {
    raw <- raw %>% mutate(event_agg = (e2 + e3 + e6)/3)
  } else {
    stop("Provide either columns e2,e3,e6 to compute event_agg or a column named event_agg.")
  }
}

dat <- raw %>%
  select(ID, week, sday, GrL, GrF, event_agg) %>%
  arrange(ID, week, sday)

## ----------------------
## Step 1: Handle missingness via time-series interpolation within ID x week
##         Caveat: Interpolation assumes MAR and can smooth/induce artifacts.
##         We recommend reporting robustness across 'linear' vs 'spline'.
## ----------------------
fill_by_group <- function(df, method = impute_option) {
  # imputeTS::na_interpolation supports option = 'linear' or 'spline'
  na_interpolation(df, option = method)
}

dat_i <- dat %>%
  group_by(ID, week) %>%
  group_modify(~ fill_by_group(.x)) %>%
  ungroup()

## Time variable per context
dat_i <- dat_i %>% mutate(TIME = sday - 1)

if (start_time_at_zero) {
  # Keep per-week time starting at 0..6 (already true with TIME = sday - 1)
  # Optionally, make continuous across weeks by subtracting 7 and 14, but
  # GOLD is computed within weeks, so not required.
}

## ----------------------
## Step 2: Linear de-trending via mixed models (level-1 residuals)
##         Residuals reflect within-person deviation from equilibrium.
##         Note: Linear trend removal centers around a time-invariant equilibrium;
##         a local level/state-space variant could be used as a sensitivity check.
## ----------------------
detres <- function(df, y) {
  # df: data frame for a single context (week)
  # y:  column name (string)
  f <- as.formula(paste0(y, " ~ TIME"))
  m <- lme(f, random = ~ TIME | ID, data = df, control = list(opt = "optim"),
           na.action = na.omit)
  residuals(m, type = "response")
}

# Split by context/week to keep time scale consistent within each block
ctx_a <- dat_i %>% filter(week == 1)
ctx_b <- dat_i %>% filter(week == 2)
ctx_c <- dat_i %>% filter(week == 3)

# Compute residuals per context for GrL, GrF, and event_agg
res_a_l <- detres(ctx_a, "GrL")
res_b_l <- detres(ctx_b, "GrL")
res_c_l <- detres(ctx_c, "GrL")

res_a_f <- detres(ctx_a, "GrF")
res_b_f <- detres(ctx_b, "GrF")
res_c_f <- detres(ctx_c, "GrF")

res_a_e <- detres(ctx_a, "event_agg")
res_b_e <- detres(ctx_b, "event_agg")
res_c_e <- detres(ctx_c, "event_agg")

# Merge residuals back
dat_i2 <- bind_rows(ctx_a, ctx_b, ctx_c) %>%
  arrange(ID, week, TIME) %>%
  mutate(
    l_res = c(res_a_l, res_b_l, res_c_l),
    f_res = c(res_a_f, res_b_f, res_c_f),
    e_res = c(res_a_e, res_b_e, res_c_e)
  )

## ----------------------
## Step 3: GOLD derivatives
##         - We operationalize frequency of fluctuation via |first derivative|
##           (velocity) and speed of return via second derivative (acceleration).
##         - Units are per time step (here, per day).
##         - For formal damping ratio (zeta) and natural frequency (omega), fit
##           a second-order ODE (e.g., regress x'' on x' and x to extract params).
##         Calculate per ID x week windows using doremi::calculate.gold
## ----------------------

gold_velocity <- function(x, t, emb = embedding_dim) {
  # Returns first derivatives for a local window (length = emb)
  if (length(t) < emb) return(NA_real_)
  g <- calculate.gold(x, t, embedding = emb)
  as.numeric(g$dsignal[seq_len(emb), 2])
}

gold_accel <- function(x, t, emb = embedding_dim) {
  # Returns second derivatives for a local window (length = emb)
  if (length(t) < emb) return(NA_real_)
  g <- calculate.gold(x, t, embedding = emb)
  as.numeric(g$dsignal[seq_len(emb), 3])
}

rep_id <- function(id, t, emb = embedding_dim) {
  if (length(t) < emb) return(NA)
  rep(unique(id), each = emb)
}

rep_week <- function(w, t, emb = embedding_dim) {
  if (length(t) < emb) return(NA)
  rep(unique(w), each = emb)
}

# Compute derivatives within ID x week for leader/follower/event residuals
feat <- dat_i2 %>%
  arrange(ID, week, TIME) %>%
  group_by(ID, week) %>%
  group_modify(~ {
    tibble(
      leader_velo = gold_velocity(.x$l_res, .x$TIME),
      follower_velo = gold_velocity(.x$f_res, .x$TIME),
      event_velo = gold_velocity(.x$e_res, .x$TIME),
      leader_accel = gold_accel(.x$l_res, .x$TIME),
      follower_accel = gold_accel(.x$f_res, .x$TIME),
      id_rep = rep_id(.x$ID, .x$TIME),
      week_rep = rep_week(.x$week, .x$TIME)
    )
  }) %>%
  ungroup() %>%
  rename(ID_out = id_rep, context = week_rep) %>%
  mutate(ID_out = factor(ID_out), context = factor(context))

# Clean NA rows created when groups lack enough observations
feat <- feat %>% drop_na(leader_velo, event_velo)

## Optional: Aggregate to per ID x context summary parameters
## (Mean absolute velocity as a rate-of-change proxy; mean acceleration as dampening proxy)
feat_summary <- feat %>%
  group_by(ID_out, context) %>%
  summarise(
    leader_freq_param = mean(abs(leader_velo), na.rm = TRUE),
    event_freq_param  = mean(abs(event_velo), na.rm = TRUE),
    leader_damp_param = mean(leader_accel, na.rm = TRUE),
    follower_damp_param = mean(follower_accel, na.rm = TRUE),
    .groups = "drop"
  )

## ----------------------
## Example models (illustrative; align with chapter narrative)
## ----------------------

message("Model: leader velocity ~ event velocity + follower velocity (random intercept for ID)")
mod_velo <- lmer(leader_velo ~ event_velo + follower_velo + (1 | ID_out), data = feat)
print(summary(mod_velo))

message("Model: leader acceleration ~ event velocity + follower acceleration (random intercept for ID)")
mod_accel <- lmer(leader_accel ~ event_velo + follower_accel + (1 | ID_out), data = feat)
print(summary(mod_accel))

## ----------------------
## Optional: Sensitivity checks over embedding and imputation choices
## ----------------------
if (isTRUE(run_sensitivity)) {
  emb_grid <- c(3, 4, 5)
  imp_grid <- c("linear", "spline")
  sens <- list()
  k <- 1
  for (emb in emb_grid) {
    for (imp in imp_grid) {
      message("Running sensitivity: embedding=", emb, ", impute=", imp)
      dat_i_s <- dat %>%
        group_by(ID, week) %>%
        group_modify(~ na_interpolation(.x, option = imp)) %>%
        ungroup() %>% mutate(TIME = sday - 1)
      a <- dat_i_s %>% filter(week == 1)
      b <- dat_i_s %>% filter(week == 2)
      c <- dat_i_s %>% filter(week == 3)
      r_al <- detres(a, "GrL"); r_bl <- detres(b, "GrL"); r_cl <- detres(c, "GrL")
      r_ae <- detres(a, "event_agg"); r_be <- detres(b, "event_agg"); r_ce <- detres(c, "event_agg")
      di2 <- bind_rows(a, b, c) %>% arrange(ID, week, TIME) %>%
        mutate(l_res = c(r_al, r_bl, r_cl), e_res = c(r_ae, r_be, r_ce))
      ft <- di2 %>% arrange(ID, week, TIME) %>% group_by(ID, week) %>%
        group_modify(~ tibble(
          l_v = gold_velocity(.x$l_res, .x$TIME, emb),
          e_v = gold_velocity(.x$e_res, .x$TIME, emb)
        )) %>% ungroup() %>% drop_na()
      sens[[k]] <- ft %>% mutate(embedding = emb, impute = imp)
      k <- k + 1
    }
  }
  sens_df <- bind_rows(sens)
  message("Sensitivity summary (mean |l_v| by setting):")
  print(sens_df %>% group_by(embedding, impute) %>% summarise(mean_abs_lv = mean(abs(l_v)), .groups = "drop"))
}

## ----------------------
## Save outputs
## ----------------------
out_feat <- here("data", "ds_features.csv")
out_sum  <- here("data", "ds_features_summary.csv")
readr::write_csv(feat, out_feat)
readr::write_csv(feat_summary, out_sum)

message("Wrote feature rows to: ", out_feat)
message("Wrote per ID x context summaries to: ", out_sum)

## ----------------------
## Quick sanity summaries
## ----------------------
print(feat %>% group_by(context) %>% summarise(across(c(leader_velo, leader_accel, event_velo), ~ mean(., na.rm = TRUE))))

if (requireNamespace("sessioninfo", quietly = TRUE)) {
  sessioninfo::session_info()
} else {
  message("Package 'sessioninfo' not installed; skipping session info.")
}
