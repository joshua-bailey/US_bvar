###############################################################################
## functions_forecasting_v6.R – Simulation, filtering & output utilities    ##
###############################################################################
## Roles of this file
##   • Convert BVAR draws to a state-space form
##   • Run Kalman filter / Carter-Kohn smoother draw-by-draw
##   • Produce unconditional and conditional forecast arrays
##   • Re-attach forecasts to historical data in raw units
###############################################################################

# ---------------------------------------------------------------------------
# 1.  State-space representation
# ---------------------------------------------------------------------------
#  Given BVAR coefficients (β) and residual covariance (Σ), build the
#  companion-form matrices
#        y_t = G x_t + v_t,   v_t ~ N(0,R)
#        x_t = F x_{t-1} + u_t, u_t ~ N(0,Q)
#
state_space <- function(B, Sig, p) {
  
  n   <- dim(B)[2]          # number of observed variables
  np1 <- n * (p + 1)        # dimension of companion state
  np  <- n * p
  
  R <- matrix(0, n, n)      # no additional measurement noise
  G <- matrix(0, n, np1); G[, 1:n] <- diag(n)
  
  Q <- matrix(0, np1, np1)
  Q[1:n, 1:n] <- Sig        # process noise only enters first n eqns
  
  F <- matrix(0, np1, np1)
  F[1:n, 1:np] <- t(B[2:(np + 1), ])         # lag coefficients
  F[1:n, (np + 1):np1] <- diag(n)            # intercept as state
  if (p > 1)
    F[(n + 1):np, 1:(n * (p - 1))] <- diag(n * (p - 1))  # shift lags
  F[(np + 1):np1, (np + 1):np1] <- diag(n)               # keep intercept
  
  list(F = F, Q = Q, G = G, R = R)
}

# ---------------------------------------------------------------------------
# 2.  Kalman filter – forward pass over history + forecast horizon
# ---------------------------------------------------------------------------
kalman <- function(F, Q, G, Y, c_init, p, h) {
  
  T <- nrow(Y)                       # total obs incl. NA forecast rows
  n <- nrow(G);  m <- nrow(F)
  
  # Build initial state: last p observations plus intercept
  x00 <- matrix(0, m, 1)
  for (i in 1:p)
    x00[(1 + n * (i - 1)):(n * i)] <- Y[T - h - i + 1, ]
  x00[(m - n + 1):m] <- c_init
  
  V00 <- matrix(0, m, m); V00[1:n, 1:n] <- diag(n)   # diffuse prior
  
  xtt <- NULL; Vtt <- list()
  
  for (t in 0:h) {
    
    ## 2.1  One-step prediction
    x01 <- F %*% x00
    V01 <- F %*% V00 %*% t(F) + Q
    
    ## 2.2  Conditioning on observed y_t (if any)
    yt  <- Y[T - h + t, ]
    ind <- which(!is.na(yt))
    if (length(ind) > 0) {
      
      Gt <- if (length(ind) == 1) matrix(G[ind, ], 1) else G[ind, ]
      Vy <- Gt %*% V01 %*% t(Gt)
      K  <- V01 %*% t(Gt) %*% solve(Vy)
      
      x00 <- x01 + K %*% (yt[ind] - Gt %*% x01)
      V00 <- V01 - K %*% Gt %*% V01
    } else {                      # no measurement update
      x00 <- x01;  V00 <- V01
    }
    
    if (t >= 1) {
      xtt <- rbind(xtt, t(x00))
      Vtt[[t]] <- V00
    }
  }
  
  list(xtt = xtt, Vtt = Vtt)
}

# ---------------------------------------------------------------------------
# 3.  Carter–Kohn smoother – single draw of x_{t|T}
# ---------------------------------------------------------------------------
smoother <- function(F, Q, xtt, Vtt, Y, h) {
  
  T <- nrow(Y);  n <- ncol(Y);  m <- ncol(xtt)
  xtT <- NULL                       # will accumulate backwards
  
  ## 3.1  Draw final-period state -------------------------------------------
  yt      <- Y[T, ]
  ind_obs <- which(!is.na(yt)); ind_mis <- which(is.na(yt))
  k       <- length(ind_mis)
  
  m0 <- matrix(xtt[h, ind_mis], k, 1)
  V0 <- Vtt[[h]][ind_mis, ind_mis]
  x_draw <- adjust_draw(yt,
                        m0 + t(chol(V0)) %*% rnorm(k))   # fill observed vars
  xtT <- rbind(xtT, t(x_draw))
  
  ## 3.2  Backward recursion -------------------------------------------------
  for (t in 1:(h - 1)) {
    
    yt       <- Y[T - t, ]
    ind_obs  <- which(!is.na(yt)); ind_mis <- which(is.na(yt))
    k        <- length(ind_mis)
    
    if (k == 0) {
      xtT <- rbind(xtT, t(yt));  next
    }
    
    Ft  <- if (k == 1) matrix(F[ind_mis, ], 1) else F[ind_mis, ]
    x00 <- matrix(xtt[h - t, ], m, 1);  V00 <- Vtt[[h - t]]
    x01 <- Ft %*% x00;  V01 <- Ft %*% V00 %*% t(Ft) + Q[ind_mis, ind_mis]
    
    x_temp  <- matrix(x_draw[ind_mis], k, 1)
    m0 <- x00 + V00 %*% t(Ft) %*% solve(V01) %*% (x_temp - x01)
    V0 <- V00 - V00 %*% t(Ft) %*% solve(V01) %*% Ft %*% V00
    
    x_draw <- adjust_draw(yt,
                          m0[ind_mis] + t(chol(V0[ind_mis, ind_mis])) %*% rnorm(k))
    xtT <- rbind(xtT, t(x_draw))
  }
  
  xtT[seq(nrow(xtT), 1, -1), ]      # reverse to chronological order
}

# Helper – merge simulated draws with observed y_t where available
adjust_draw <- function(yt, x0) {
  ind_obs <- which(!is.na(yt)); ind_mis <- which(is.na(yt))
  out <- matrix(NA_real_, length(yt), 1)
  out[ind_obs] <- yt[ind_obs];  out[ind_mis] <- x0
  out
}

# ---------------------------------------------------------------------------
# 4.  Unconditional forecast (all draws)                                     #
# ---------------------------------------------------------------------------
generate_unconditional <- function(data, df, res, p, h, n_sim, var_names0) {
  
  Y        <- rbind(as.matrix(data), matrix(NA, h, ncol(data)))
  end_ind  <- nrow(data)
  forecasts<- array(0, c(h, ncol(Y), n_sim))
  
  pb <- txtProgressBar(min = 1, max = n_sim, style = 3)
  for (i in seq_len(n_sim)) {
    beta  <- res$beta[i, , ]; sigma <- res$sigma[i, , ]
    fcst_i<- generate_forecasts_one_batch(
      Y, df, beta, sigma, p, h, var_names0, end_ind,
      c_init = beta[1, ])
    forecasts[ , , i] <- fcst_i
    setTxtProgressBar(pb, i)
  }
  cat("\n")
  forecasts
}

# ---------------------------------------------------------------------------
# 5.  Scenario forecast (future path partly pinned)                          #
# ---------------------------------------------------------------------------
generate_one_scenario <- function(Y, end_ind, df, res, p, h, n_sim, var_names0,
                                  base, scenario_name, verbose = TRUE) {
  
  forecasts <- array(0, c(h, ncol(Y), n_sim))
  if (verbose) pb <- txtProgressBar(min = 1, max = n_sim, style = 3)
  
  for (i in seq_len(n_sim)) {
    beta <- res$beta[i, , ]; sigma <- res$sigma[i, , ]
    base0<- if (scenario_name == "base") 0 else base[ , , i]
    
    forecasts[ , , i] <- generate_forecasts_one_batch(
      Y, df, beta, sigma, p, h,
      var_names0, end_ind, c_init = beta[1, ],
      base0 = base0)
    if (verbose) setTxtProgressBar(pb, i)
  }
  if (verbose) cat("\nfinished scenario:", scenario_name, "\n")
  forecasts
}

# ---------------------------------------------------------------------------
# 6.  Core workhorse – one Kalman-smoother pass for a single draw            #
# ---------------------------------------------------------------------------
generate_forecasts_one_batch <- function(
    Y, df, beta, sigma, p, h,
    var_names0, end_ind, c_init,
    base0 = 0) {
  
  ss  <- state_space(beta, sigma, p)
  km  <- kalman(ss$F, ss$Q, ss$G, Y, c_init, p, h)
  xtT <- smoother(ss$F, ss$Q, km$xtt, km$Vtt, Y, h)
  
  forecasts <- xtT
  if (!identical(base0, 0))
    forecasts <- forecasts + base0      # offset for conditional paths
  
  forecasts                         # (h × n_vars) matrix
}

# ---------------------------------------------------------------------------
# 7.  Attach forecasts to history in raw units                               #
# ---------------------------------------------------------------------------
combine_fcst_with_history_raw <- function(fcst_array,
                                          df_transformed,
                                          h, trans,
                                          var_names, var_names0,
                                          center = c("median", "mean")) {
  
  center <- match.arg(center)
  fcst_summary <- if (center == "median")
    apply(fcst_array, c(1, 2), median)
  else
    apply(fcst_array, c(1, 2), mean)
  
  # --- 7.1  Helper to undo model transforms -------------------------------
  revert_transforms_df <- function(df_in, var_names, var_names0, trans) {
    out <- df_in
    for (j in seq_along(var_names)) {
      v <- var_names[j]
      if (trans[j] == "log100") out[[v]] <- exp(out[[v]] / 100)
      if (trans[j] == "raw")    out[[v]] <- out[[v]] / 100
      if (var_names0[j] == "Potential GDP")
        out[[v]] <- out[[v]] / 10
    }
    out
  }
  
  # --- 7.2  Back-transform history & forecasts -----------------------------
  df_hist_raw <- revert_transforms_df(df_transformed, var_names, var_names0, trans)
  
  df_fcst <- as.data.frame(fcst_summary);  names(df_fcst) <- var_names
  df_fcst_raw <- revert_transforms_df(df_fcst, var_names, var_names0, trans)
  
  future_dates <- seq(tail(df_hist_raw$yq, 1) + 0.25,
                      by   = 0.25, length.out = h)
  df_fcst_raw$yq <- future_dates
  
  bind_rows(df_hist_raw, df_fcst_raw)
}

# ---------------------------------------------------------------------------
# 8.  Unused functions to potentially add back later                                      #
# ---------------------------------------------------------------------------
#  • entropic_tilting_pce()  – importance-sample draws to hit a PCE target
#  • check_explosive() / check_stability()  – eigenvalue diagnostics
#    (left untouched; original internal logic retained)

entropic_tilting_pce=function(forecast,var_names0){
  ind=which(var_names0=='PCE Price Index')
  h=12
  gbar=2
  Y=forecast[h,ind,]-forecast[(h-4),ind,]
  
  get_gamma=function(gamma,Y,gbar){
    gbar0=rep(gbar,length(Y))
    y=sum(exp(gamma*(Y-gbar0)))
    return(y)
  }
  
  optimizer=optim(1.1, get_gamma, Y=Y, gbar=gbar, method="BFGS")
  gamma=optimizer$par
  wstar=exp(gamma*Y)/sum(exp(gamma*Y))
  n_sim=dim(forecast)[3]
  index=sample(c(1:n_sim),(n_sim),replace=TRUE,prob=wstar)
  return(forecast[,,index])
}

check_explosive=function(res){
  beta=apply(res$beta,c(2,3),median)
  sigma=apply(res$sigma,c(2,3),mean)
  p=res$meta$lags
  ss=state_space(beta,sigma,p)
  lam_max=max(abs(eigen(ss$F)$values))
  return(lam_max>1.0001)
}

check_stability=function(res){
  is_stable=function(i,res){
    beta=res$beta[i,,]
    sigma=res$sigma[i,,]
    p=res$meta$lags
    ss=state_space(beta,sigma,p)
    lam_max=max(abs(eigen(ss$F)$values))
    outcome=ifelse(lam_max<=1.0001,1,0)
    return(outcome)
  }
  n_sim=dim(res$beta)[1]
  inds_stable=sapply(1:n_sim,is_stable,res)
  return(inds_stable)
}
