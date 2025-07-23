########################
## functions_forecasting_v6.R
########################

#-----------------------------------
# 1. building a state-space form from bvar coefficients
#-----------------------------------
state_space <- function(B, Sig, p) {
  # we want f, q, g, r in a companion form:
  #   y_t = g * x_t + v_t,  v_t ~ n(0, r)
  #   x_t = f * x_{t-1} + u_t,  u_t ~ n(0, q)
  #
  # b is (p+1) x n, with b[1,] as intercept
  # sig is n x n
  # f is (n*(p+1)) x (n*(p+1))
  #
  # we place the top block of f with the lag coefficients,
  # and the bottom block is mostly identity to move states forward.
  
  n   <- dim(B)[2]        # number of variables
  np1 <- n * (p + 1)      # companion state dimension
  np  <- n * p
  
  # r = 0 since measurement eq is y_t = x_t[1:n] + 0
  R <- matrix(0, n, n)
  
  # g picks off the first n elements of x_t to interpret as y_t
  G <- matrix(0, n, np1)
  G[, 1:n] <- diag(n)
  
  # q is (np1 x np1), with sig in the top-left
  Q <- matrix(0, np1, np1)
  Q[1:n, 1:n] <- Sig
  
  # build f
  F <- matrix(0, np1, np1)
  # b[2:(np+1), ] is the stacked lag coefficients
  F[1:n, 1:np] <- t(B[2:(np + 1), ])
  # identity to move intercept forward
  F[1:n, (np + 1):np1] <- diag(n)
  # if p>1, shift older lags
  if (p > 1) {
    F[(n + 1):np, 1:(n * (p - 1))] <- diag(n * (p - 1))
  }
  # the last n just carry forward
  F[(np + 1):np1, (np + 1):np1] <- diag(n)
  
  return(list(F = F, Q = Q, G = G, R = R))
}

#-----------------------------------
# 2. kalman filter for forecast extension
#-----------------------------------
kalman <- function(F, Q, G, Y, c_init, p, h) {
  # we do a forward pass to get x_{t|t} for t = T+1 ... T+h
  # c_init is the intercept row (b[1, ])
  # p is number of lags
  # y is dimension (t x n); we add rows of na for future quarters
  #
  # returns:
  #   xtt -> the sequence of filtered states for t=1..h
  #   vtt -> the sequence of var-cov for each t
  
  T <- nrow(Y)   # total time steps in y
  n <- nrow(G)
  m <- nrow(F)   # dimension of x_t
  
  # we build the initial state x00 from the last p obs plus intercept
  # for the new forecast steps
  x00 <- matrix(0, m, 1)
  for (i in 1:p) {
    # fill in the older obs in reverse
    x00[(1 + n * (i - 1)):(n * i)] <- Y[T - h - i + 1, ]
  }
  # place the intercept in the last n elements
  x00[(m - n + 1):m] <- c_init
  
  # build initial variance
  V00 <- matrix(0, m, m)
  # set the intercept part to identity
  V00[1:n, 1:n] <- diag(n)
  
  # placeholders
  xtt <- c()       # will store x_{t|t} for t=1..h
  Vtt <- list()
  
  for (t in 0:h) {
    # one-step predict
    x01 <- F %*% x00
    V01 <- F %*% V00 %*% t(F) + Q
    
    # next we update with Y at time T-h+t
    # note that T-h is the last historical point minus the forecast horizon
    # so T-h+t moves from T-h up to T
    yt  <- Y[T - h + t, ]
    ind <- which(!is.na(yt))
    
    if (length(ind) == 0) {
      # means all y are na, so no update
      # store predicted state
      if (t >= 1) {
        xtt <- rbind(xtt, t(x01))
        Vtt[[t]] <- V01
      }
      x00 <- x01
      V00 <- V01
      next
    } else {
      # measurement update
      Gt <- G[ind, ]  # partial pick if some subset is non-na
      if (length(ind) == 1) {
        Gt <- matrix(G[ind, ], nrow = 1)
      }
      
      Vy <- Gt %*% V01 %*% t(Gt)
      K  <- V01 %*% t(Gt) %*% solve(Vy)
      x00 <- x01 + K %*% (yt[ind] - Gt %*% x01)
      V00 <- V01 - K %*% Gt %*% V01
      
      # store after update
      if (t >= 1) {
        xtt <- rbind(xtt, t(x00))
        Vtt[[t]] <- V00
      }
    }
  }
  
  return(list(xtt = xtt, Vtt = Vtt))
}

#-----------------------------------
# 3. smoother (carter-kohn approach)
#-----------------------------------
smoother <- function(F, Q, xtt, Vtt, Y, h) {
  # we do a backward pass starting from t = T
  # T = last row of Y, h steps of forecast
  # xtt, Vtt are from the kalman filter forward pass
  #
  # returns draws of x_{t|T}, but we take just one draw
  # we keep a median path or single draw as needed
  
  T <- nrow(Y)
  n <- ncol(Y)
  m <- ncol(xtt)
  xtT <- c()
  
  # final period
  yt    <- Y[T, ]
  ind   <- which(!is.na(yt))
  ind_na<- which(is.na(yt))
  k     <- length(ind_na)
  
  # step 1: draw x_{T|T}
  m0 <- matrix(xtt[h, ind_na], k, 1)
  V0 <- Vtt[[h]][ind_na, ind_na]
  # random draw from n(m0, v0)
  x_draw0 <- m0 + t(chol(V0)) %*% rnorm(k)
  x_draw  <- adjust_draw(yt, x_draw0)  # fill in known obs
  xtT     <- rbind(xtT, t(x_draw))
  
  # now go backward t=1..(h-1)
  for (t in 1:(h - 1)) {
    yt  <- Y[T - t, ]
    ind <- which(!is.na(yt))
    ind_na <- which(is.na(yt))
    k <- length(ind_na)
    
    if (k == 0) {
      # no missing, we store actual
      xtT <- rbind(xtT, t(yt))
      next
    }
    
    # partial f for these missing vars
    if (k > 1) {
      Ft <- F[ind_na, ]
    } else {
      Ft <- matrix(F[ind_na, ], nrow = 1)
    }
    
    x00 <- matrix(xtt[h - t, ], m, 1)
    x01 <- Ft %*% x00
    V00 <- Vtt[[h - t]]
    
    # q subset
    qt  <- Q[ind_na, ind_na]
    V01 <- Ft %*% V00 %*% t(Ft) + qt
    
    x_temp <- matrix(x_draw[ind_na], k, 1)
    m0     <- x00 + V00 %*% t(Ft) %*% solve(V01) %*% (x_temp - x01)
    V0     <- V00 - V00 %*% t(Ft) %*% solve(V01) %*% Ft %*% V00
    
    x_draw0 <- m0[ind_na] + t(chol(V0[ind_na, ind_na])) %*% rnorm(k)
    x_draw  <- adjust_draw(yt, x_draw0)
    xtT     <- rbind(xtT, t(x_draw))
  }
  
  # reverse order so we get x_{T-h+1|T} ... x_{T|T}
  return(xtT[seq(nrow(xtT), 1, -1), ])
}

#-----------------------------------
# 4. helper for smoother draws
#-----------------------------------
adjust_draw <- function(yt, x0) {
  # we take a draw for the missing variables, but if a var is known (not na),
  # we fill it in from yt
  # so we combine them in a single vector
  ind     <- which(!is.na(yt))
  ind_na  <- which(is.na(yt))
  x1      <- matrix(NA, length(yt), 1)
  x1[ind] <- yt[ind]
  x1[ind_na] <- x0
  return(x1)
}

#-----------------------------------
# 5. generate_unconditional forecast
#-----------------------------------
generate_unconditional <- function(data, df, res, p, h, n_sim, var_names0) {
  # this function creates an unconditional forecast by:
  #   1) building y with future n/a
  #   2) for each draw, do a forward/backward pass
  #
  # returns a 3d array: [h x n_vars x n_sim]
  #
  Y       <- rbind(as.matrix(data), matrix(NA, h, ncol(data)))
  end_ind <- nrow(data)
  forecasts <- array(0, c(h, ncol(Y), n_sim))
  
  # we get the intercept from the first row of each draw
  # and do the full procedure
  pb <- txtProgressBar(min = 1, max = n_sim, style = 3)
  for (i in seq_len(n_sim)) {
    beta  <- res$beta[i, , ]
    sigma <- res$sigma[i, , ]
    c_init<- beta[1, ]  # intercept
    
    fcst_i <- generate_forecasts_one_batch(Y, df, beta, sigma, p, h, var_names0, end_ind, c_init)
    forecasts[,, i] <- fcst_i
    setTxtProgressBar(pb, i)
  }
  cat("\n")
  
  return(forecasts)
}

#-----------------------------------
# 6. generate_one_scenario forecast
#-----------------------------------
generate_one_scenario <- function(Y, end_ind, df, res, p, h, n_sim, var_names0,
                                  base, scenario_name, verbose = TRUE) {
  # this function pins or shifts the future in some way
  # but we let external functions (in functions_conditions_v6.R) do the pinning.
  # base is the unconditional scenario if scenario_name != 'base'
  #
  # returns a 3d array: [h x n_vars x n_sim]
  
  forecasts <- array(0, c(h, ncol(Y), n_sim))
  if (verbose) pb <- txtProgressBar(min = 1, max = n_sim, style = 3)
  
  for (i in seq_len(n_sim)) {
    beta   <- res$beta[i, , ]
    sigma  <- res$sigma[i, , ]
    c_init <- beta[1, ]
    
    # if scenario is 'base', we set base0 = 0
    # otherwise base is the unconditional forecast for that draw
    if (scenario_name == "base") {
      base0 <- 0
    } else {
      base0 <- base[,, i]
    }
    
    fcst_i <- generate_forecasts_one_batch(Y, df, beta, sigma, p, h, var_names0, end_ind, c_init,
                                           base0 = base0)
    forecasts[,, i] <- fcst_i
    
    if (verbose) setTxtProgressBar(pb, i)
  }
  if (verbose) cat("\nfinished scenario:", scenario_name, "\n")
  
  return(forecasts)
}

# -------------------------------------------------------------------
#  7. generate_forecasts_one_batch  – draw‑by‑draw Kalman smoother forecast
# -------------------------------------------------------------------
generate_forecasts_one_batch <- function(
    Y, df, beta, sigma, p, h,
    var_names0, end_ind, c_init,
    base0 = 0            # offset path (median of baseline if supplied)
){
  # Build state‑space representation
  ss  <- state_space(beta, sigma, p)
  
  # Forward filter + backward smoother
  km  <- kalman(ss$F, ss$Q, ss$G, Y, c_init, p, h)
  xtT <- smoother(ss$F, ss$Q, km$xtt, km$Vtt, Y, h)
  
  # Extract the top‑n elements (forecasts for original y_t)
  forecasts <- xtT      # dimensions: h × n_vars
  
  # Apply base‑path offset if provided (used by scenario functions)
  if (!identical(base0, 0)) {
    forecasts <- forecasts + base0
  }
  
  return(forecasts)
}


#-----------------------------------
# 8. combine_fcst_with_history_raw
#-----------------------------------
combine_fcst_with_history_raw <- function(fcst_array, df_transformed, h, trans, var_names, var_names0,
                                          center = c("median", "mean")) {
  # this function merges the forecast array with the historical df,
  # reverting the log/100 scaling back to real-world levels.
  #
  # fcst_array is [h x n_vars x n_sim]
  # df_transformed is the historical data in transformed space (e.g. log, raw)
  # we pick median or mean across the n_sim dimension
  # then we tack it on to the df (by creating future yq)
  
  center <- match.arg(center)
  if (center == "median") {
    fcst_summary <- apply(fcst_array, c(1, 2), median)
  } else {
    fcst_summary <- apply(fcst_array, c(1, 2), mean)
  }
  
  # revert the transforms
  revert_transforms_df <- function(df_in, var_names, var_names0, trans) {
    df_out <- df_in
    for (j in seq_along(var_names)) {
      col_j <- var_names[j]
      if (trans[j] == "log100") {
        df_out[[col_j]] <- exp(df_out[[col_j]] / 100)
      } else if (trans[j] == "raw") {
        df_out[[col_j]] <- df_out[[col_j]] / 100
      }
      # potential gdp division by 10
      if (var_names0[j] == "Potential GDP") {
        df_out[[col_j]] <- df_out[[col_j]] / 10
      }
    }
    return(df_out)
  }
  
  # revert historical
  df_history_raw <- revert_transforms_df(df_transformed, var_names, var_names0, trans)
  
  # revert forecast summary
  df_fcst <- as.data.frame(fcst_summary)
  colnames(df_fcst) <- var_names
  df_fcst_raw <- revert_transforms_df(df_fcst, var_names, var_names0, trans)
  
  # build future yq
  last_obs <- tail(df_history_raw$yq, 1)
  future_dates <- seq(as.yearqtr(last_obs) + 0.25, by = 0.25, length.out = h)
  df_fcst_raw$yq <- future_dates
  
  # combine
  df_combined <- bind_rows(df_history_raw, df_fcst_raw)
  return(df_combined)
}

#-----------------------------------
# Unused functions to potentially add back later 
#-----------------------------------
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


