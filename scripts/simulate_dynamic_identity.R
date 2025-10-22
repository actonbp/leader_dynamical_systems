############################################################
# Simulate dynamic identity data with damped oscillatory dynamics
# for leader (GrL) and follower (GrF) identities plus event components
# (e2, e3, e6). Produces 3 contexts (weeks) of 7 days each per ID.
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

simulate_identity_data <- function(
  n_id = 40,
  n_context = 3,
  days_per_ctx = 7,
  seed = 123,
  # oscillator parameters (heterogeneous across IDs)
  omega_mean = 0.9, omega_sd = 0.2,      # natural frequency (radians/day)
  zeta_mean = 0.3, zeta_sd = 0.1,        # damping ratio
  event_influence = 0.5,                 # effect of event on leader identity
  follow_coupling = -0.2,                # coupling of follower onto leader change
  noise_sd = 0.3,                        # process noise
  miss_rate = 0.08                       # probability a day is missing per variable
) {
  set.seed(seed)

  total_days <- n_context * days_per_ctx
  id_df <- tibble(ID = 1:n_id)

  # Per-ID parameters
  pars <- id_df %>% mutate(
    omega = pmax(0.2, rnorm(n_id, omega_mean, omega_sd)),
    zeta  = pmin(0.9, pmax(0.05, rnorm(n_id, zeta_mean, zeta_sd))),
    baseL = rnorm(n_id, 4.5, 0.6),  # base level for GrL
    baseF = rnorm(n_id, 4.5, 0.6)   # base level for GrF
  )

  sim_id <- function(id_row) {
    ID <- id_row$ID; omega <- id_row$omega; zeta <- id_row$zeta
    baseL <- id_row$baseL; baseF <- id_row$baseF

    # Discrete-time parameters from continuous damped oscillator
    dt <- 1
    phi <- exp(-zeta * omega * dt)
    omega_d <- omega * sqrt(pmax(1e-6, 1 - zeta^2))
    a1 <- 2 * phi * cos(omega_d * dt)
    a2 <- -phi^2

    # Events per context: trending/volatile patterns
    days <- 1:total_days
    week <- rep(1:n_context, each = days_per_ctx)
    sday <- rep(1:days_per_ctx, times = n_context)

    # Generate an event driver with context-varying structure
    ev_base <- rnorm(total_days, 0, 1)
    ev_ctx <- ev_base + rep(c(0.5, -0.3, 0.0), each = days_per_ctx) +
      cumsum(rnorm(total_days, 0, 0.1))
    # Split into components to emulate e2/e3/e6
    e2 <- scale(ev_ctx + rnorm(total_days, 0, 0.2))[,1]
    e3 <- scale(0.8 * ev_ctx + rnorm(total_days, 0, 0.3))[,1]
    e6 <- scale(0.6 * ev_ctx + rnorm(total_days, 0, 0.4))[,1]

    event_agg <- (e2 + e3 + e6) / 3

    # Initialize identities
    L <- numeric(total_days)
    F <- numeric(total_days)
    L[1:2] <- baseL + rnorm(2, 0, 0.2)
    F[1:2] <- baseF + rnorm(2, 0, 0.2)

    # Recurrence with event input and follower coupling
    for (t in 3:total_days) {
      input <- event_influence * event_agg[t - 1] + follow_coupling * (F[t - 1] - baseF)
      L[t] <- baseL + a1 * (L[t - 1] - baseL) + a2 * (L[t - 2] - baseL) + input + rnorm(1, 0, noise_sd)
      # Follower identity as slower, less oscillatory AR(1) with event influence
      F[t] <- baseF + 0.6 * (F[t - 1] - baseF) + 0.2 * event_agg[t - 1] + rnorm(1, 0, 0.25)
    }

    tibble(ID = ID, week = week, sday = sday, GrL = L, GrF = F, e2 = as.numeric(e2), e3 = as.numeric(e3), e6 = as.numeric(e6))
  }

  sim <- purrr::map_dfr(seq_len(n_id), ~ sim_id(pars[., ]))

  # Inject missingness (MCAR for demo)
  mask <- function(x) ifelse(runif(length(x)) < miss_rate, NA, x)
  sim <- sim %>% mutate(
    GrL = mask(GrL), GrF = mask(GrF),
    e2 = mask(e2), e3 = mask(e3), e6 = mask(e6)
  )

  sim
}

