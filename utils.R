design_matrix <- function(model, var = NULL, val = NULL) {
 # Returns the design matrix. If var and val are provided, we make an
 # intervention beforehand.
 dat <- model$data
 if (!is.null(var) & !is.null(val)) {
  if (is.numeric(dat[[var]])) {
   dat <- mutate(dat, !!var := val)
  }
  else if (is.character(dat[[var]]) | is.factor(dat[[var]])) {
   dat <- mutate(dat, !!sym(var) := factor(val, levels = levels(dat[[var]])))
  }
  else {
   stop("Variable must be a character, factor or numeric.")
  }
 }
 res <- model.matrix(formula(model), data = dat)
 attr(res, "assign") <- NULL
 attr(res, "contrasts") <- NULL
 return(res)
}

dmu_fct <- function(model) {
 # The derivative dmu/deta
 # Same as make.link(model$family$link)$mu.eta
 switch(model$family$link,
        "identity" = {\(eta) 1},
        "log" = {\(eta) exp(eta)},
        "logit" = {\(eta) 1/(2 + exp(eta) + exp(-eta))},
        "inverse" = {\(eta) -1/(eta^2)},
        "probit" = {\(eta) dnorm(eta)},
        "cloglog" = {\(eta) exp(eta)*exp(-exp(eta))}
 )
}

d2mu_fct <- function(model) {
 # The second order derivative d2mu/d2eta
 switch(model$family$link,
        "identity" = {\(eta) 0},
        "log" = {\(eta) exp(eta)},
        "logit" = {\(eta) 2*exp(3*eta)/(1 + exp(eta))^3 + 1/(exp(-eta) + 1) - 3*exp(2*eta)/(1 + exp(eta))^2},
        "inverse" = {\(eta) 2/(eta^3)},
        "probit" = {\(eta) -eta*dnorm(eta)},
        "cloglog" = {\(eta) exp(eta)*exp(-exp(eta))*(1 - exp(eta))}
 )
}

v_fct <- function(model) {
 # The variance function
 # Same as model$family$variance
 switch(model$family$family,
        "gaussian" = {\(mu) 1},
        "binomial" = {\(mu) mu*(1 - mu)},
        "Gamma" = {\(mu) mu^2}
 )
}

dv_fct <- function(model) {
 # The derivative dv/dmu
 switch(model$family$family,
        "gaussian" = {\(mu) 0},
        "binomial" = {\(mu) 1 - 2*mu},
        "Gamma" = {\(mu) 2*mu}
 )
}

M_mat <- function(model, oim = TRUE) { #oim: use observed information matrix, else use expected
 link <- model$family$linkfun
 beta <- model$coefficients
 y <- model$y
 mu <- model$fitted.values
 eta <- link(mu)
 dmu <- purrr::map_vec(eta, dmu_fct(model))
 d2mu <- purrr::map_vec(eta, d2mu_fct(model))
 v <- purrr::map_vec(mu, v_fct(model))
 dv <- purrr::map_vec(mu, dv_fct(model))
 X_obs <- design_matrix(model)
 n_obs <- nobs(model)
 M_middle <- -dmu*dmu/v
 if (oim) {
  M_middle <- M_middle + (d2mu/v - dv*dmu*dmu/(v*v))*(y - mu)
 }
 return(t(X_obs) %*% diag(M_middle) %*% X_obs / n_obs) # Could also do t(X_bs*M_middle) %*% X_obs / n_obs without diagnonalizing M_middle. Maybe that's faster.
}

sigma_mat <- function(model) {
 link <- model$family$linkfun
 beta <- model$coefficients
 y <- model$model[, 1]
 mu <- model$fitted.values
 eta <- link(mu)
 dmu <- purrr::map_vec(eta, dmu_fct(model))
 d2mu <- purrr::map_vec(eta, d2mu_fct(model))
 v <- purrr::map_vec(mu, v_fct(model))
 dv <- purrr::map_vec(mu, dv_fct(model))
 X_obs <- design_matrix(model)
 n_obs <- nobs(model)
 sigma_middle <- (dmu*(y - mu)/v)^2
 return(t(X_obs) %*% diag(sigma_middle) %*% X_obs / n_obs)
}

beta_infl_fct <- function(model, oim = TRUE) {
 link <- model$family$linkfun
 beta <- model$coefficients
 y <- model$model[, 1]
 mu <- model$fitted.values
 eta <- link(mu)
 dmu <- purrr::map_vec(eta, dmu_fct(model))
 v <- purrr::map_vec(mu, v_fct(model))
 X_obs <- design_matrix(model)
 M_inv <- solve(M_mat(model, oim))
 return(-((X_obs*dmu/v) %*% M_inv)*(y - mu))
}

po_var <- function(model, var, vals, oim = TRUE) {
 n_obs <- nobs(model)
 beta <- model$coefficients
 beta_infl <- beta_infl_fct(model, oim)
 
 X_int <- purrr::map2(
  rep(var, length(vals)),
  vals,
  \(var, val) design_matrix(model, var, val)
 )
 
 po <- list()
 dmu_int <- list()
 for (i in 1:length(vals)) {
  eta <- X_int[[i]] %*% beta
  po[[i]] <- model$family$linkinv(eta)
  dmu_int[[i]] <- map_vec(eta, dmu_fct(model))
 }
 
 po_mean <- vector()
 po_mean_infl <- matrix(data = NA, nrow = n_obs, ncol = length(vals))
 po_mean_infl_pt1 <- matrix(data = NA, nrow = n_obs, ncol = length(vals))
 po_mean_infl_pt2 <- matrix(data = NA, nrow = n_obs, ncol = length(vals))
 
 for (i in 1:length(vals)) {
  po_mean[i] <- mean(po[[i]])
  po_mean_infl_pt1[, i] <- po[[i]] - po_mean[[i]]
  po_mean_infl_pt2[, i] <- beta_infl %*% colMeans(X_int[[i]]*dmu_int[[i]])
 }
 
 po_var_no_cross <- (t(po_mean_infl_pt1) %*% po_mean_infl_pt1 + t(po_mean_infl_pt2) %*% po_mean_infl_pt2)/n_obs
 po_var_pt1 <- (t(po_mean_infl_pt1) %*% po_mean_infl_pt1)/n_obs 
 po_var_pt2 <- (t(po_mean_infl_pt2) %*% po_mean_infl_pt2)/n_obs 
 po_var_pt1pt2 <- (t(po_mean_infl_pt1) %*% po_mean_infl_pt2)/n_obs + (t(po_mean_infl_pt2) %*% po_mean_infl_pt1)/n_obs
 
 list(
  po_mean = po_mean,
  po_var_pt1 = po_var_pt1,
  po_var_pt2 = po_var_pt2,
  po_var_pt1pt2 = po_var_pt1pt2,
  po_var_total = po_var_pt1 + po_var_pt2 + po_var_pt1pt2,
  n_obs = n_obs
 )
}

po_se <- function(model, var, vals, oim = FALSE) {
 po_res <- po_var(model, var, vals, oim)
 n_obs <- po_res$n_obs
 po_se <- sqrt(diag(po_res$po_var_total)/n_obs)
 po_se_fixed <- sqrt(diag(po_res$po_var_pt2)/n_obs)
 tibble(
  {{var}} := vals,
  pom = po_res$po_mean,
  se = po_se,
  se_fixed = po_se_fixed
 )
}

