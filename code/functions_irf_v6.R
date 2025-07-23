########################
## functions_irf_v6.R
########################

#----------------------------------------------
# generate_irf: structural impulse responses
#----------------------------------------------
#  res       -> bvar results (from estimate_bvar)
#  p         -> number of lags used in estimation
#  horizon   -> number of periods for the IRF
#  shock_var -> index or name of the variable that receives the unit shock
#  ident     -> identification scheme (only 'chol' currently)
#
#  returns a 3d array: [horizon+1 x n_vars x n_draws]
#----------------------------------------------

generate_irf <- function(res, p, horizon = 12, shock_var, ident = c("chol"),
                         shock_size = 1) {
  ident <- match.arg(ident)
  draws <- dim(res$beta)[1]
  n     <- dim(res$beta)[3]
  irfs  <- array(0, c(horizon + 1, n, draws))

  # determine index of the shocked variable
  if (is.character(shock_var)) {
    if (!is.null(res$meta$variables)) {
      shock_idx <- match(shock_var, res$meta$variables)
    } else {
      shock_idx <- match(shock_var, colnames(res$beta[1, , ]))
    }
  } else {
    shock_idx <- shock_var
  }
  if (is.na(shock_idx) || shock_idx < 1 || shock_idx > n) {
    stop("Shock variable not found")
  }

  for (d in seq_len(draws)) {
    B     <- res$beta[d, , ]
    Sigma <- res$sigma[d, , ]
    ss    <- state_space(B, Sigma, p)
    Fm    <- ss$F
    Gm    <- ss$G

    A0 <- switch(ident,
                 chol = t(chol(Sigma)),
                 diag(n))

    e0 <- rep(0, n)
    e0[shock_idx] <- shock_size
    shock <- A0 %*% e0

    np1 <- nrow(Fm)
    x   <- matrix(0, np1, horizon + 1)
    x[1:n, 1] <- shock

    for (h_i in 2:(horizon + 1)) {
      x[, h_i] <- Fm %*% x[, h_i - 1]
    }
    y <- Gm %*% x
    irfs[, , d] <- t(y)
  }

  return(irfs)
}

#----------------------------------------------
# plot_irf_core: convenience plotting
#----------------------------------------------
#  irf_array -> output from generate_irf
#  var_names0-> variable names (in estimation order)
#  shock_name-> name for plot title
#  center    -> 'median' or 'mean' across draws
#----------------------------------------------

plot_irf_core <- function(irf_array, var_names0, shock_name = "", 
                          center = c("median", "mean")) {
  center <- match.arg(center)
  if (center == "median") {
    irf <- apply(irf_array, c(1, 2), median)
  } else {
    irf <- apply(irf_array, c(1, 2), mean)
  }

  df <- as.data.frame(irf)
  colnames(df) <- var_names0
  df$horizon <- 0:(nrow(df) - 1)
  df_long <- df %>%
    pivot_longer(-horizon, names_to = "variable", values_to = "value")

  ggplot(df_long, aes(x = horizon, y = value)) +
    geom_line(color = "blue") +
    facet_wrap(~variable, scales = "free_y") +
    labs(title = paste0("IRF to shock in ", shock_name),
         x = "Quarter", y = "Response") +
    theme_bw()
}
