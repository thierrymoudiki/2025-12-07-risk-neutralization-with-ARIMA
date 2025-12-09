## ----setup, include=FALSE---------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----1-simulate-SVJD, cache=FALSE-------------------------------------------------------------
# ARIMA-Black-Scholes: Semi-Parametric Risk-Neutral Pricing
# T. Moudiki
# 2025-12-04

library(esgtoolkit)
library(forecast)
library(ahead)

# =============================================================================
# 1 - SVJD SIMULATION (Physical Measure)
# =============================================================================

set.seed(123)
n <- 100L
h <- 5
freq <- "daily"
r <- 0.05
maturity <- 5
S0 <- 100
mu <- 0.08
sigma <- 0.04

# Simulate under physical measure with stochastic volatility and jumps
sim_GBM <- esgtoolkit::simdiff(n=n, horizon = h, frequency = freq, x0=S0, theta1 = mu, theta2 = sigma)
sim_SVJD <- esgtoolkit::rsvjd(n=n, r0=mu)
sim_Heston <- esgtoolkit::rsvjd(n=n, r0=mu, 
                                lambda = 0,
                                mu_J = 0,
                                sigma_J = 0)

cat("Simulation dimensions:\n")
cat("Start:", start(sim_SVJD), "\n")
cat("End:", end(sim_SVJD), "\n")
cat("Paths:", ncol(sim_SVJD), "\n")
cat("Time steps:", nrow(sim_SVJD), "\n\n")

par(mfrow=c(1, 3))
# Plot historical (physical measure) paths
esgtoolkit::esgplotbands(sim_GBM, main="GBM Paths under the Physical Measure", 
                         xlab="Time",
                         ylab="Stock prices")
esgtoolkit::esgplotbands(sim_SVJD, main="SVJD Paths under the Physical Measure", 
                         xlab="Time",
                         ylab="Stock prices")
esgtoolkit::esgplotbands(sim_Heston, main="Heston Paths under the Physical Measure", 
                         xlab="Time",
                         ylab="Stock prices")                         

# Summary statistics
cat("Physical measure statistics (GBM):\n")
terminal_prices_physical_GBM <- sim_GBM[nrow(sim_GBM), ]
cat("Mean terminal price:", mean(terminal_prices_physical_GBM), "\n")
cat("Std terminal price:", sd(terminal_prices_physical_GBM), "\n")
cat("Expected under P:", S0 * exp(mu * maturity), "\n\n")
 
cat("Physical measure statistics (SVJD):\n")
terminal_prices_physical_SVJD <- sim_SVJD[nrow(sim_SVJD), ]
cat("Mean terminal price:", mean(terminal_prices_physical_SVJD), "\n")
cat("Std terminal price:", sd(terminal_prices_physical_SVJD), "\n")
cat("Expected under P:", S0 * exp(mu * maturity), "\n\n")

cat("Physical measure statistics (Heston):\n")
terminal_prices_physical_Heston <- sim_Heston[nrow(sim_Heston), ]
cat("Mean terminal price:", mean(terminal_prices_physical_Heston), "\n")
cat("Std terminal price:", sd(terminal_prices_physical_Heston), "\n")
cat("Expected under P:", S0 * exp(mu * maturity), "\n\n")


## ----2-compute-discounted, cache=FALSE--------------------------------------------------------
# =============================================================================
# 2 - COMPUTE DISCOUNTED PRICES (Transform to Martingale Domain)
# =============================================================================

discounted_prices_GBM <- esgtoolkit::esgdiscountfactor(r=r, X=sim_GBM)
discounted_prices_SVJD <- esgtoolkit::esgdiscountfactor(r=r, X=sim_SVJD)
discounted_prices_Heston <- esgtoolkit::esgdiscountfactor(r=r, X=sim_Heston)
martingale_diff_GBM <- discounted_prices_GBM - S0 
martingale_diff_SVJD <- discounted_prices_SVJD - S0 
martingale_diff_Heston <- discounted_prices_Heston - S0 

cat("Martingale differences dimensions (GBM):", dim(martingale_diff_GBM), "\n")
cat("Mean martingale diff (should be ≠ 0 under P):\n")
print(t.test(rowMeans(martingale_diff_GBM)))

cat("\nMartingale differences dimensions (SVJD):", dim(martingale_diff_SVJD), "\n")
cat("Mean martingale diff (should be ≠ 0 under P):\n")
print(t.test(rowMeans(martingale_diff_SVJD)))

cat("\nMartingale differences dimensions (Heston):", dim(martingale_diff_Heston), "\n")
cat("Mean martingale diff (should be ≠ 0 under P):\n")
print(t.test(rowMeans(martingale_diff_Heston)))

# =============================================================================
# 3 - VISUALIZE RISK PREMIUM
# =============================================================================

par(mfrow=c(2,2))

matplot(discounted_prices_GBM, type='l', col=rgb(0,0,1,0.3),
        main="Discounted Stock Prices (Physical Measure - GBM)",
        ylab="exp(-rt) * S_t", xlab="Time step")
abline(h=S0, col='red', lwd=2, lty=2)

matplot(discounted_prices_SVJD, type='l', col=rgb(0,0,1,0.3),
        main="Discounted Stock Prices (Physical Measure - SVJD)",
        ylab="exp(-rt) * S_t", xlab="Time step")
abline(h=S0, col='red', lwd=2, lty=2)

matplot(discounted_prices_Heston, type='l', col=rgb(0,0,1,0.3),
        main="Discounted Stock Prices (Physical Measure - Heston)",
        ylab="exp(-rt) * S_t", xlab="Time step")
abline(h=S0, col='red', lwd=2, lty=2)

par(mfrow=c(1,1))

mean_disc_path_GBM <- rowMeans(discounted_prices_GBM)
times_plot <- as.numeric(time(discounted_prices_GBM)) 
plot(times_plot, mean_disc_path_GBM, type='l', lwd=2, col='blue',
     main="Risk Premium in Discounted Prices (GBM)",
     xlab="Time (years)", ylab="E[exp(-rt)*S_t]")
abline(h=S0, col='red', lwd=2, lty=2)
lines(times_plot, S0 * exp((mu-r)*times_plot), col='green', lwd=2, lty=3)
legend("topleft",
       legend=c("Empirical mean", "S0", "Theoretical (μ-r drift)"),
       col=c('blue','red','green'), lty=c(1,2,3), lwd=2)

mean_disc_path_SVJD <- rowMeans(discounted_prices_SVJD)
times_plot <- as.numeric(time(discounted_prices_SVJD)) 
plot(times_plot, mean_disc_path_SVJD, type='l', lwd=2, col='blue',
     main="Risk Premium in Discounted Prices (SVJD)",
     xlab="Time (years)", ylab="E[exp(-rt)*S_t]")
abline(h=S0, col='red', lwd=2, lty=2)
lines(times_plot, S0 * exp((mu-r)*times_plot), col='green', lwd=2, lty=3)
legend("topleft",
       legend=c("Empirical mean", "S0", "Theoretical (μ-r drift)"),
       col=c('blue','red','green'), lty=c(1,2,3), lwd=2)

mean_disc_path_Heston <- rowMeans(discounted_prices_Heston)
times_plot <- as.numeric(time(discounted_prices_Heston)) 
plot(times_plot, mean_disc_path_Heston, type='l', lwd=2, col='blue',
     main="Risk Premium in Discounted Prices (Heston)",
     xlab="Time (years)", ylab="E[exp(-rt)*S_t]")
abline(h=S0, col='red', lwd=2, lty=2)
lines(times_plot, S0 * exp((mu-r)*times_plot), col='green', lwd=2, lty=3)
legend("topleft",
       legend=c("Empirical mean", "S0", "Theoretical (μ-r drift)"),
       col=c('blue','red','green'), lty=c(1,2,3), lwd=2)     


## ----3-fit-ARIMA, cache=FALSE, eval=TRUE------------------------------------------------------
# =============================================================================
# 4 - FIT ARIMA MODELS TO EXTRACT RISK PREMIUM
# =============================================================================

n_periods <- nrow(martingale_diff_GBM)
n_paths <- ncol(martingale_diff_GBM)

martingale_increments_GBM <- diff(martingale_diff_GBM)
martingale_increments_SVJD <- diff(martingale_diff_SVJD)
martingale_increments_Heston <- diff(martingale_diff_Heston)

# Initialize storage arrays
arima_residuals_GBM <- array(NA, dim = c(nrow(martingale_increments_GBM), n_paths))
centered_arima_residuals_GBM <- array(NA, dim = c(nrow(martingale_increments_GBM), n_paths))
means_arima_residuals_GBM <- rep(NA, n_paths)
arima_models_GBM <- list()

arima_residuals_SVJD <- array(NA, dim = c(nrow(martingale_increments_SVJD), n_paths))
centered_arima_residuals_SVJD <- array(NA, dim = c(nrow(martingale_increments_SVJD), n_paths))
means_arima_residuals_SVJD <- rep(NA, n_paths)
arima_models_SVJD <- list()

arima_residuals_Heston <- array(NA, dim = c(nrow(martingale_increments_Heston), n_paths))
centered_arima_residuals_Heston <- array(NA, dim = c(nrow(martingale_increments_Heston), n_paths))
means_arima_residuals_Heston <- rep(NA, n_paths)
arima_models_Heston <- list()

# Fit ARIMA models to GBM
cat("Fitting ARIMA models to", n_paths, "GBM paths...\n")
for (i in 1:n_paths) {
  y <- as.numeric(martingale_increments_GBM[, i])
  fit <- forecast::thetaf(y)
  arima_models_GBM[[i]] <- fit
  
  res <- as.numeric(residuals(fit))
  arima_residuals_GBM[, i] <- res
  
  centre_arima_residuals <- scale(res, center = TRUE, scale = FALSE)
  means_arima_residuals_GBM[i] <- attr(centre_arima_residuals, "scaled:center")
  centered_arima_residuals_GBM[, i] <- centre_arima_residuals[,1]
}

# Fit ARIMA models to SVJD
cat("Fitting ARIMA models to", n_paths, "SVJD paths...\n")
for (i in 1:n_paths) {
  y <- as.numeric(martingale_increments_SVJD[, i])
  fit <- forecast::thetaf(y)
  arima_models_SVJD[[i]] <- fit
  
  res <- as.numeric(residuals(fit))
  arima_residuals_SVJD[, i] <- res
  
  centre_arima_residuals <- scale(res, center = TRUE, scale = FALSE)
  means_arima_residuals_SVJD[i] <- attr(centre_arima_residuals, "scaled:center")
  centered_arima_residuals_SVJD[, i] <- centre_arima_residuals[,1]
}

# Fit ARIMA models to Heston
cat("Fitting ARIMA models to", n_paths, "Heston paths...\n")
for (i in 1:n_paths) {
  y <- as.numeric(martingale_increments_Heston[, i])
  fit <- forecast::thetaf(y)
  arima_models_Heston[[i]] <- fit
  
  res <- as.numeric(residuals(fit))
  arima_residuals_Heston[, i] <- res
  
  centre_arima_residuals <- scale(res, center = TRUE, scale = FALSE)
  means_arima_residuals_Heston[i] <- attr(centre_arima_residuals, "scaled:center")
  centered_arima_residuals_Heston[, i] <- centre_arima_residuals[,1]
}

# Box-Ljung tests
pvalues_GBM <- sapply(1:ncol(centered_arima_residuals_GBM), 
                  function(i) Box.test(centered_arima_residuals_GBM[,i])$p.value)
cat("\nBox-Ljung test p-values (GBM):\n")
cat("Mean p-value:", mean(pvalues_GBM), "\n")
cat("Proportion > 0.05:", mean(pvalues_GBM > 0.05), "\n")

pvalues_SVJD <- sapply(1:ncol(centered_arima_residuals_SVJD), 
                  function(i) Box.test(centered_arima_residuals_SVJD[,i])$p.value)
cat("\nBox-Ljung test p-values (SVJD):\n")
cat("Mean p-value:", mean(pvalues_SVJD), "\n")
cat("Proportion > 0.05:", mean(pvalues_SVJD > 0.05), "\n")

pvalues_Heston <- sapply(1:ncol(centered_arima_residuals_Heston), 
                  function(i) Box.test(centered_arima_residuals_Heston[,i])$p.value)
cat("\nBox-Ljung test p-values (Heston):\n")
cat("Mean p-value:", mean(pvalues_Heston), "\n")
cat("Proportion > 0.05:", mean(pvalues_Heston > 0.05), "\n")

par(mfrow=c(1,3))
hist(pvalues_GBM, breaks=20, col='lightgreen',
     main="Box-Ljung P-values (GBM)",
     xlab="P-value")
abline(v=0.05, col='red', lwd=2, lty=2)

hist(pvalues_SVJD, breaks=20, col='lightblue',
     main="Box-Ljung P-values (SVJD)",
     xlab="P-value")
abline(v=0.05, col='red', lwd=2, lty=2)

hist(pvalues_Heston, breaks=20, col='lightcoral',
     main="Box-Ljung P-values (Heston)",
     xlab="P-value")
abline(v=0.05, col='red', lwd=2, lty=2)
par(mfrow=c(1,1))


## ----4-generate-rn-paths, cache=FALSE, eval=TRUE----------------------------------------------
# =============================================================================
# 5 - GENERATE RISK-NEUTRAL PATHS
# =============================================================================

bootstrap_independent <- function(x, R = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- length(x)
  
  # Each column is an independent bootstrap sample
  replicate(R, sample(x, size = n, replace = TRUE))
}

cat("\n\nGenerating risk-neutral paths from ALL historical paths...\n")
n_sim_per_path <- 20  # Generate 10 paths per historical path
times <- seq(0, maturity, length.out = n_periods)
discount_factor <- exp(r * times)

# Storage for all risk-neutral paths
all_S_tilde_GBM <- list()
all_S_tilde_SVJD <- list()
all_S_tilde_Heston <- list()

# Generate GBM risk-neutral paths
for (i in 1:n_paths) {
  resampled_residuals <- bootstrap_independent(centered_arima_residuals_GBM[, i], 
                                               R = n_sim_per_path)
  fit <- arima_models_GBM[[i]]
  fitted_increments <- as.numeric(fitted(fit))
  
  discounted_path <- matrix(0, nrow = n_periods, ncol = n_sim_per_path)
  discounted_path[1, ] <- S0
  
  increments <- matrix(scale(fitted_increments, center = TRUE, 
                             scale = FALSE)[,1], nrow = n_periods - 1, ncol = n_sim_per_path) + 
                resampled_residuals[1:(n_periods - 1), ]
  
  discounted_path[-1, ] <- S0 + apply(increments, 2, cumsum)
  S_tilde_price <- discounted_path * discount_factor
  all_S_tilde_GBM[[i]] <- S_tilde_price
}

# Generate SVJD risk-neutral paths
for (i in 1:n_paths) {
  resampled_residuals <- bootstrap_independent(centered_arima_residuals_SVJD[, i], 
                                               R = n_sim_per_path)
  fit <- arima_models_SVJD[[i]]
  fitted_increments <- as.numeric(fitted(fit))
  
  discounted_path <- matrix(0, nrow = n_periods, ncol = n_sim_per_path)
  discounted_path[1, ] <- S0
  
  increments <- matrix(scale(fitted_increments, center = TRUE, 
                             scale = FALSE)[,1], nrow = n_periods - 1, ncol = n_sim_per_path) + 
                resampled_residuals[1:(n_periods - 1), ]
  
  discounted_path[-1, ] <- S0 + apply(increments, 2, cumsum)
  S_tilde_price <- discounted_path * discount_factor
  all_S_tilde_SVJD[[i]] <- S_tilde_price
}

# Generate Heston risk-neutral paths
for (i in 1:n_paths) {
  resampled_residuals <- bootstrap_independent(centered_arima_residuals_Heston[, i], 
                                               R = n_sim_per_path)
  fit <- arima_models_Heston[[i]]
  fitted_increments <- as.numeric(fitted(fit))
  
  discounted_path <- matrix(0, nrow = n_periods, ncol = n_sim_per_path)
  discounted_path[1, ] <- S0
  
  increments <- matrix(scale(fitted_increments, center = TRUE, 
                             scale = FALSE)[,1], nrow = n_periods - 1, ncol = n_sim_per_path) + 
                resampled_residuals[1:(n_periods - 1), ]
  
  discounted_path[-1, ] <- S0 + apply(increments, 2, cumsum)
  S_tilde_price <- discounted_path * discount_factor
  all_S_tilde_Heston[[i]] <- S_tilde_price
}

# Combine all paths
S_tilde_combined_GBM <- do.call(cbind, all_S_tilde_GBM)
cat("Total GBM risk-neutral paths generated:", ncol(S_tilde_combined_GBM), "\n")

S_tilde_combined_SVJD <- do.call(cbind, all_S_tilde_SVJD)
cat("Total SVJD risk-neutral paths generated:", ncol(S_tilde_combined_SVJD), "\n")

S_tilde_combined_Heston <- do.call(cbind, all_S_tilde_Heston)
cat("Total Heston risk-neutral paths generated:", ncol(S_tilde_combined_Heston), "\n\n")

# Convert to time series
S_tilde_ts_GBM <- ts(S_tilde_combined_GBM, start = start(sim_GBM), 
                     frequency = frequency(sim_GBM))
S_tilde_ts_SVJD <- ts(S_tilde_combined_SVJD, start = start(sim_SVJD), 
                      frequency = frequency(sim_SVJD))
S_tilde_ts_Heston <- ts(S_tilde_combined_Heston, start = start(sim_Heston), 
                        frequency = frequency(sim_Heston))

# Visualize risk-neutral paths
par(mfrow=c(1,3))
esgtoolkit::esgplotbands(S_tilde_ts_GBM, 
                        main="Risk-Neutral Paths - GBM")
esgtoolkit::esgplotbands(S_tilde_ts_SVJD, 
                        main="Risk-Neutral Paths - SVJD")
esgtoolkit::esgplotbands(S_tilde_ts_Heston, 
                        main="Risk-Neutral Paths - Heston")

# Sample plots
par(mfrow=c(1,3))
matplot(S_tilde_combined_GBM[, sample(ncol(S_tilde_combined_GBM), 200)], 
        type='l', col=rgb(0,0,1,0.1),
        main="Risk-Neutral Paths (GBM)",
        xlab="Time step", ylab="Stock Price")
lines(rowMeans(S_tilde_combined_GBM), col='red', lwd=3)
abline(h=S0, col='green', lwd=2, lty=2)

matplot(S_tilde_combined_SVJD[, sample(ncol(S_tilde_combined_SVJD), 200)], 
        type='l', col=rgb(0,0,1,0.1),
        main="Risk-Neutral Paths (SVJD)",
        xlab="Time step", ylab="Stock Price")
lines(rowMeans(S_tilde_combined_SVJD), col='red', lwd=3)
abline(h=S0, col='green', lwd=2, lty=2)

matplot(S_tilde_combined_Heston[, sample(ncol(S_tilde_combined_Heston), 200)], 
        type='l', col=rgb(0,0,1,0.1),
        main="Risk-Neutral Paths (Heston)",
        xlab="Time step", ylab="Stock Price")
lines(rowMeans(S_tilde_combined_Heston), col='red', lwd=3)
abline(h=S0, col='green', lwd=2, lty=2)
par(mfrow=c(1,1))


## ----5-rn-verif, cache=FALSE, eval=TRUE-------------------------------------------------------
# =============================================================================
# 6 - VERIFY RISK-NEUTRAL PROPERTY
# =============================================================================

cat("\n=== RISK-NEUTRAL VERIFICATION ===\n\n")

terminal_prices_rn_GBM <- S_tilde_combined_GBM[n_periods, ]
terminal_prices_rn_SVJD <- S_tilde_combined_SVJD[n_periods, ]
terminal_prices_rn_Heston <- S_tilde_combined_Heston[n_periods, ]
capitalized_stock_price <- S0 * exp(r * maturity)

cat("GBM Risk-Neutral Verification:\n")
cat("Expected terminal price (Q):", capitalized_stock_price, "\n")
cat("Empirical mean:", mean(terminal_prices_rn_GBM), "\n")
cat("Difference:", mean(terminal_prices_rn_GBM) - capitalized_stock_price, "\n")
print(t.test(terminal_prices_rn_GBM - capitalized_stock_price))

cat("\nSVJD Risk-Neutral Verification:\n")
cat("Expected terminal price (Q):", capitalized_stock_price, "\n")
cat("Empirical mean:", mean(terminal_prices_rn_SVJD), "\n")
cat("Difference:", mean(terminal_prices_rn_SVJD) - capitalized_stock_price, "\n")
print(t.test(terminal_prices_rn_SVJD - capitalized_stock_price))

cat("\nHeston Risk-Neutral Verification:\n")
cat("Expected terminal price (Q):", capitalized_stock_price, "\n")
cat("Empirical mean:", mean(terminal_prices_rn_Heston), "\n")
cat("Difference:", mean(terminal_prices_rn_Heston) - capitalized_stock_price, "\n")
print(t.test(terminal_prices_rn_Heston - capitalized_stock_price))

# Visualization comparison
par(mfrow=c(3, 2))
hist(terminal_prices_physical_GBM, breaks=30, col=rgb(1,0,0,0.5),
     main="Terminal Prices: Physical (GBM)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_physical_GBM), col='red', lwd=2)
abline(v=S0*exp(mu*maturity), col='blue', lwd=2, lty=2)

hist(terminal_prices_rn_GBM, breaks=30, col=rgb(0,0,1,0.5),
     main="Terminal Prices: Risk-Neutral (GBM)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_rn_GBM), col='blue', lwd=2)
abline(v=S0*exp(r*maturity), col='red', lwd=2, lty=2)

hist(terminal_prices_physical_SVJD, breaks=30, col=rgb(1,0,0,0.5),
     main="Terminal Prices: Physical (SVJD)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_physical_SVJD), col='red', lwd=2)
abline(v=S0*exp(mu*maturity), col='blue', lwd=2, lty=2)

hist(terminal_prices_rn_SVJD, breaks=30, col=rgb(0,0,1,0.5),
     main="Terminal Prices: Risk-Neutral (SVJD)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_rn_SVJD), col='blue', lwd=2)
abline(v=S0*exp(r*maturity), col='red', lwd=2, lty=2)

hist(terminal_prices_physical_Heston, breaks=30, col=rgb(1,0,0,0.5),
     main="Terminal Prices: Physical (Heston)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_physical_Heston), col='red', lwd=2)
abline(v=S0*exp(mu*maturity), col='blue', lwd=2, lty=2)

hist(terminal_prices_rn_Heston, breaks=30, col=rgb(0,0,1,0.5),
     main="Terminal Prices: Risk-Neutral (Heston)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_rn_Heston), col='blue', lwd=2)
abline(v=S0*exp(r*maturity), col='red', lwd=2, lty=2)
par(mfrow=c(1,1))


## ----6-option-pricing, cache=FALSE, eval=TRUE-------------------------------------------------
# =============================================================================
# 7 - OPTION PRICING
# =============================================================================

cat("\n=== OPTION PRICING ===\n\n")

bs_price <- function(S, K, r, sigma, T, q = 0) {

  d1 <- (log(S / K) + (r - q + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)

  call <- S * exp(-q * T) * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  put  <- K * exp(-r * T) * pnorm(-d2) - S * exp(-q * T) * pnorm(-d1)

  list(call = call, put = put)
}


strikes <- seq(80, 160, by=10)
d_f <- exp(-r * maturity)

# Function to price options
price_options <- function(terminal_prices, strikes, discount_factor) {
  n_strikes <- length(strikes)
  call_prices <- numeric(n_strikes)
  bs_call_prices <- numeric(n_strikes)
  put_prices <- numeric(n_strikes)
  bs_put_prices <- numeric(n_strikes)
  call_se <- numeric(n_strikes)
  put_se <- numeric(n_strikes)
  
  for (k in 1:n_strikes) {
    K <- strikes[k]
    call_payoffs <- pmax(terminal_prices - K, 0)
    call_prices[k] <- mean(call_payoffs) * discount_factor
    call_se[k] <- sd(call_payoffs) / sqrt(length(call_payoffs)) * discount_factor
    
    put_payoffs <- pmax(K - terminal_prices, 0)
    put_prices[k] <- mean(put_payoffs) * discount_factor
    put_se[k] <- sd(put_payoffs) / sqrt(length(put_payoffs)) * discount_factor
    
    bs_prices <- bs_price(S0, K, r, sigma=sigma, T=5, q = 0)
    bs_call_prices[k] <- bs_prices$call
    bs_put_prices[k] <- bs_prices$put
  }
  
  list(call_prices = call_prices, put_prices = put_prices,
       bs_call_prices = bs_call_prices, bs_put_prices = bs_put_prices,
       call_se = call_se, put_se = put_se)
}

# Price options for all models
options_GBM <- price_options(terminal_prices_rn_GBM, strikes, d_f)
options_SVJD <- price_options(terminal_prices_rn_SVJD, strikes, d_f)
options_Heston <- price_options(terminal_prices_rn_Heston, strikes, d_f)


## ----7-recap-options, cache=FALSE, eval=TRUE--------------------------------------------------
kableExtra::kable(as.data.frame(options_GBM))
kableExtra::kable(as.data.frame(options_Heston))
kableExtra::kable(as.data.frame(options_SVJD))


## ----3-fit-ARIMA-2, cache=FALSE, eval=TRUE----------------------------------------------------
# =============================================================================
# 4 - FIT ARIMA MODELS TO EXTRACT RISK PREMIUM
# =============================================================================

n_periods <- nrow(martingale_diff_GBM)
n_paths <- ncol(martingale_diff_GBM)

martingale_increments_GBM <- diff(martingale_diff_GBM)
martingale_increments_SVJD <- diff(martingale_diff_SVJD)
martingale_increments_Heston <- diff(martingale_diff_Heston)

# Initialize storage arrays
arima_residuals_GBM <- array(NA, dim = c(nrow(martingale_increments_GBM), n_paths))
centered_arima_residuals_GBM <- array(NA, dim = c(nrow(martingale_increments_GBM), n_paths))
means_arima_residuals_GBM <- rep(NA, n_paths)
arima_models_GBM <- list()

arima_residuals_SVJD <- array(NA, dim = c(nrow(martingale_increments_SVJD), n_paths))
centered_arima_residuals_SVJD <- array(NA, dim = c(nrow(martingale_increments_SVJD), n_paths))
means_arima_residuals_SVJD <- rep(NA, n_paths)
arima_models_SVJD <- list()

arima_residuals_Heston <- array(NA, dim = c(nrow(martingale_increments_Heston), n_paths))
centered_arima_residuals_Heston <- array(NA, dim = c(nrow(martingale_increments_Heston), n_paths))
means_arima_residuals_Heston <- rep(NA, n_paths)
arima_models_Heston <- list()

# Fit ARIMA models to GBM
cat("Fitting ARIMA models to", n_paths, "GBM paths...\n")
for (i in 1:n_paths) {
  y <- as.numeric(martingale_increments_GBM[, i])
  fit <- forecast::ets(y)
  arima_models_GBM[[i]] <- fit
  
  res <- as.numeric(residuals(fit))
  arima_residuals_GBM[, i] <- res
  
  centre_arima_residuals <- scale(res, center = TRUE, scale = FALSE)
  means_arima_residuals_GBM[i] <- attr(centre_arima_residuals, "scaled:center")
  centered_arima_residuals_GBM[, i] <- centre_arima_residuals[,1]
}

# Fit ARIMA models to SVJD
cat("Fitting ARIMA models to", n_paths, "SVJD paths...\n")
for (i in 1:n_paths) {
  y <- as.numeric(martingale_increments_SVJD[, i])
  fit <- forecast::ets(y)
  arima_models_SVJD[[i]] <- fit
  
  res <- as.numeric(residuals(fit))
  arima_residuals_SVJD[, i] <- res
  
  centre_arima_residuals <- scale(res, center = TRUE, scale = FALSE)
  means_arima_residuals_SVJD[i] <- attr(centre_arima_residuals, "scaled:center")
  centered_arima_residuals_SVJD[, i] <- centre_arima_residuals[,1]
}

# Fit ARIMA models to Heston
cat("Fitting ARIMA models to", n_paths, "Heston paths...\n")
for (i in 1:n_paths) {
  y <- as.numeric(martingale_increments_Heston[, i])
  fit <- forecast::thetaf(y)
  arima_models_Heston[[i]] <- fit
  
  res <- as.numeric(residuals(fit))
  arima_residuals_Heston[, i] <- res
  
  centre_arima_residuals <- scale(res, center = TRUE, scale = FALSE)
  means_arima_residuals_Heston[i] <- attr(centre_arima_residuals, "scaled:center")
  centered_arima_residuals_Heston[, i] <- centre_arima_residuals[,1]
}

# Box-Ljung tests
pvalues_GBM <- sapply(1:ncol(centered_arima_residuals_GBM), 
                  function(i) Box.test(centered_arima_residuals_GBM[,i])$p.value)
cat("\nBox-Ljung test p-values (GBM):\n")
cat("Mean p-value:", mean(pvalues_GBM), "\n")
cat("Proportion > 0.05:", mean(pvalues_GBM > 0.05), "\n")

pvalues_SVJD <- sapply(1:ncol(centered_arima_residuals_SVJD), 
                  function(i) Box.test(centered_arima_residuals_SVJD[,i])$p.value)
cat("\nBox-Ljung test p-values (SVJD):\n")
cat("Mean p-value:", mean(pvalues_SVJD), "\n")
cat("Proportion > 0.05:", mean(pvalues_SVJD > 0.05), "\n")

pvalues_Heston <- sapply(1:ncol(centered_arima_residuals_Heston), 
                  function(i) Box.test(centered_arima_residuals_Heston[,i])$p.value)
cat("\nBox-Ljung test p-values (Heston):\n")
cat("Mean p-value:", mean(pvalues_Heston), "\n")
cat("Proportion > 0.05:", mean(pvalues_Heston > 0.05), "\n")

par(mfrow=c(1,3))
hist(pvalues_GBM, breaks=20, col='lightgreen',
     main="Box-Ljung P-values (GBM)",
     xlab="P-value")
abline(v=0.05, col='red', lwd=2, lty=2)

hist(pvalues_SVJD, breaks=20, col='lightblue',
     main="Box-Ljung P-values (SVJD)",
     xlab="P-value")
abline(v=0.05, col='red', lwd=2, lty=2)

hist(pvalues_Heston, breaks=20, col='lightcoral',
     main="Box-Ljung P-values (Heston)",
     xlab="P-value")
abline(v=0.05, col='red', lwd=2, lty=2)
par(mfrow=c(1,1))


## ----4-generate-rn-paths-2, cache=FALSE, eval=TRUE--------------------------------------------
# =============================================================================
# 5 - GENERATE RISK-NEUTRAL PATHS
# =============================================================================

bootstrap_independent <- function(x, R = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- length(x)
  
  # Each column is an independent bootstrap sample
  replicate(R, sample(x, size = n, replace = TRUE))
}

cat("\n\nGenerating risk-neutral paths from ALL historical paths...\n")
n_sim_per_path <- 20  # Generate 10 paths per historical path
times <- seq(0, maturity, length.out = n_periods)
discount_factor <- exp(r * times)

# Storage for all risk-neutral paths
all_S_tilde_GBM <- list()
all_S_tilde_SVJD <- list()
all_S_tilde_Heston <- list()

# Generate GBM risk-neutral paths
for (i in 1:n_paths) {
  resampled_residuals <- bootstrap_independent(centered_arima_residuals_GBM[, i], 
                                               R = n_sim_per_path)
  fit <- arima_models_GBM[[i]]
  fitted_increments <- as.numeric(fitted(fit))
  
  discounted_path <- matrix(0, nrow = n_periods, ncol = n_sim_per_path)
  discounted_path[1, ] <- S0
  
  increments <- matrix(scale(fitted_increments, center = TRUE, 
                             scale = FALSE)[,1], nrow = n_periods - 1, ncol = n_sim_per_path) + 
                resampled_residuals[1:(n_periods - 1), ]
  
  discounted_path[-1, ] <- S0 + apply(increments, 2, cumsum)
  S_tilde_price <- discounted_path * discount_factor
  all_S_tilde_GBM[[i]] <- S_tilde_price
}

# Generate SVJD risk-neutral paths
for (i in 1:n_paths) {
  resampled_residuals <- bootstrap_independent(centered_arima_residuals_SVJD[, i], 
                                               R = n_sim_per_path)
  fit <- arima_models_SVJD[[i]]
  fitted_increments <- as.numeric(fitted(fit))
  
  discounted_path <- matrix(0, nrow = n_periods, ncol = n_sim_per_path)
  discounted_path[1, ] <- S0
  
  increments <- matrix(scale(fitted_increments, center = TRUE, 
                             scale = FALSE)[,1], nrow = n_periods - 1, ncol = n_sim_per_path) + 
                resampled_residuals[1:(n_periods - 1), ]
  
  discounted_path[-1, ] <- S0 + apply(increments, 2, cumsum)
  S_tilde_price <- discounted_path * discount_factor
  all_S_tilde_SVJD[[i]] <- S_tilde_price
}

# Generate Heston risk-neutral paths
for (i in 1:n_paths) {
  resampled_residuals <- bootstrap_independent(centered_arima_residuals_Heston[, i], 
                                               R = n_sim_per_path)
  fit <- arima_models_Heston[[i]]
  fitted_increments <- as.numeric(fitted(fit))
  
  discounted_path <- matrix(0, nrow = n_periods, ncol = n_sim_per_path)
  discounted_path[1, ] <- S0
  
  increments <- matrix(scale(fitted_increments, center = TRUE, 
                             scale = FALSE)[,1], nrow = n_periods - 1, ncol = n_sim_per_path) + 
                resampled_residuals[1:(n_periods - 1), ]
  
  discounted_path[-1, ] <- S0 + apply(increments, 2, cumsum)
  S_tilde_price <- discounted_path * discount_factor
  all_S_tilde_Heston[[i]] <- S_tilde_price
}

# Combine all paths
S_tilde_combined_GBM <- do.call(cbind, all_S_tilde_GBM)
cat("Total GBM risk-neutral paths generated:", ncol(S_tilde_combined_GBM), "\n")

S_tilde_combined_SVJD <- do.call(cbind, all_S_tilde_SVJD)
cat("Total SVJD risk-neutral paths generated:", ncol(S_tilde_combined_SVJD), "\n")

S_tilde_combined_Heston <- do.call(cbind, all_S_tilde_Heston)
cat("Total Heston risk-neutral paths generated:", ncol(S_tilde_combined_Heston), "\n\n")

# Convert to time series
S_tilde_ts_GBM <- ts(S_tilde_combined_GBM, start = start(sim_GBM), 
                     frequency = frequency(sim_GBM))
S_tilde_ts_SVJD <- ts(S_tilde_combined_SVJD, start = start(sim_SVJD), 
                      frequency = frequency(sim_SVJD))
S_tilde_ts_Heston <- ts(S_tilde_combined_Heston, start = start(sim_Heston), 
                        frequency = frequency(sim_Heston))

# Visualize risk-neutral paths
par(mfrow=c(1,3))
esgtoolkit::esgplotbands(S_tilde_ts_GBM, 
                        main="Risk-Neutral Paths - GBM")
esgtoolkit::esgplotbands(S_tilde_ts_SVJD, 
                        main="Risk-Neutral Paths - SVJD")
esgtoolkit::esgplotbands(S_tilde_ts_Heston, 
                        main="Risk-Neutral Paths - Heston")

# Sample plots
par(mfrow=c(1,3))
matplot(S_tilde_combined_GBM[, sample(ncol(S_tilde_combined_GBM), 200)], 
        type='l', col=rgb(0,0,1,0.1),
        main="Risk-Neutral Paths (GBM)",
        xlab="Time step", ylab="Stock Price")
lines(rowMeans(S_tilde_combined_GBM), col='red', lwd=3)
abline(h=S0, col='green', lwd=2, lty=2)

matplot(S_tilde_combined_SVJD[, sample(ncol(S_tilde_combined_SVJD), 200)], 
        type='l', col=rgb(0,0,1,0.1),
        main="Risk-Neutral Paths (SVJD)",
        xlab="Time step", ylab="Stock Price")
lines(rowMeans(S_tilde_combined_SVJD), col='red', lwd=3)
abline(h=S0, col='green', lwd=2, lty=2)

matplot(S_tilde_combined_Heston[, sample(ncol(S_tilde_combined_Heston), 200)], 
        type='l', col=rgb(0,0,1,0.1),
        main="Risk-Neutral Paths (Heston)",
        xlab="Time step", ylab="Stock Price")
lines(rowMeans(S_tilde_combined_Heston), col='red', lwd=3)
abline(h=S0, col='green', lwd=2, lty=2)
par(mfrow=c(1,1))


## ----5-rn-verif-2, cache=FALSE, eval=TRUE-----------------------------------------------------
# =============================================================================
# 6 - VERIFY RISK-NEUTRAL PROPERTY
# =============================================================================

cat("\n=== RISK-NEUTRAL VERIFICATION ===\n\n")

terminal_prices_rn_GBM <- S_tilde_combined_GBM[n_periods, ]
terminal_prices_rn_SVJD <- S_tilde_combined_SVJD[n_periods, ]
terminal_prices_rn_Heston <- S_tilde_combined_Heston[n_periods, ]
capitalized_stock_price <- S0 * exp(r * maturity)

cat("GBM Risk-Neutral Verification:\n")
cat("Expected terminal price (Q):", capitalized_stock_price, "\n")
cat("Empirical mean:", mean(terminal_prices_rn_GBM), "\n")
cat("Difference:", mean(terminal_prices_rn_GBM) - capitalized_stock_price, "\n")
print(t.test(terminal_prices_rn_GBM - capitalized_stock_price))

cat("\nSVJD Risk-Neutral Verification:\n")
cat("Expected terminal price (Q):", capitalized_stock_price, "\n")
cat("Empirical mean:", mean(terminal_prices_rn_SVJD), "\n")
cat("Difference:", mean(terminal_prices_rn_SVJD) - capitalized_stock_price, "\n")
print(t.test(terminal_prices_rn_SVJD - capitalized_stock_price))

cat("\nHeston Risk-Neutral Verification:\n")
cat("Expected terminal price (Q):", capitalized_stock_price, "\n")
cat("Empirical mean:", mean(terminal_prices_rn_Heston), "\n")
cat("Difference:", mean(terminal_prices_rn_Heston) - capitalized_stock_price, "\n")
print(t.test(terminal_prices_rn_Heston - capitalized_stock_price))

# Visualization comparison
par(mfrow=c(3, 2))
hist(terminal_prices_physical_GBM, breaks=30, col=rgb(1,0,0,0.5),
     main="Terminal Prices: Physical (GBM)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_physical_GBM), col='red', lwd=2)
abline(v=S0*exp(mu*maturity), col='blue', lwd=2, lty=2)

hist(terminal_prices_rn_GBM, breaks=30, col=rgb(0,0,1,0.5),
     main="Terminal Prices: Risk-Neutral (GBM)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_rn_GBM), col='blue', lwd=2)
abline(v=S0*exp(r*maturity), col='red', lwd=2, lty=2)

hist(terminal_prices_physical_SVJD, breaks=30, col=rgb(1,0,0,0.5),
     main="Terminal Prices: Physical (SVJD)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_physical_SVJD), col='red', lwd=2)
abline(v=S0*exp(mu*maturity), col='blue', lwd=2, lty=2)

hist(terminal_prices_rn_SVJD, breaks=30, col=rgb(0,0,1,0.5),
     main="Terminal Prices: Risk-Neutral (SVJD)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_rn_SVJD), col='blue', lwd=2)
abline(v=S0*exp(r*maturity), col='red', lwd=2, lty=2)

hist(terminal_prices_physical_Heston, breaks=30, col=rgb(1,0,0,0.5),
     main="Terminal Prices: Physical (Heston)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_physical_Heston), col='red', lwd=2)
abline(v=S0*exp(mu*maturity), col='blue', lwd=2, lty=2)

hist(terminal_prices_rn_Heston, breaks=30, col=rgb(0,0,1,0.5),
     main="Terminal Prices: Risk-Neutral (Heston)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_rn_Heston), col='blue', lwd=2)
abline(v=S0*exp(r*maturity), col='red', lwd=2, lty=2)
par(mfrow=c(1,1))


## ----6-option-pricing-2, cache=FALSE, eval=TRUE-----------------------------------------------
# =============================================================================
# 7 - OPTION PRICING
# =============================================================================

cat("\n=== OPTION PRICING ===\n\n")

bs_price <- function(S, K, r, sigma, T, q = 0) {

  d1 <- (log(S / K) + (r - q + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)

  call <- S * exp(-q * T) * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  put  <- K * exp(-r * T) * pnorm(-d2) - S * exp(-q * T) * pnorm(-d1)

  list(call = call, put = put)
}


strikes <- seq(80, 160, by=10)
d_f <- exp(-r * maturity)

# Function to price options
price_options <- function(terminal_prices, strikes, discount_factor) {
  n_strikes <- length(strikes)
  call_prices <- numeric(n_strikes)
  bs_call_prices <- numeric(n_strikes)
  put_prices <- numeric(n_strikes)
  bs_put_prices <- numeric(n_strikes)
  call_se <- numeric(n_strikes)
  put_se <- numeric(n_strikes)
  
  for (k in 1:n_strikes) {
    K <- strikes[k]
    call_payoffs <- pmax(terminal_prices - K, 0)
    call_prices[k] <- mean(call_payoffs) * discount_factor
    call_se[k] <- sd(call_payoffs) / sqrt(length(call_payoffs)) * discount_factor
    
    put_payoffs <- pmax(K - terminal_prices, 0)
    put_prices[k] <- mean(put_payoffs) * discount_factor
    put_se[k] <- sd(put_payoffs) / sqrt(length(put_payoffs)) * discount_factor
    
    bs_prices <- bs_price(S0, K, r, sigma=sigma, T=5, q = 0)
    bs_call_prices[k] <- bs_prices$call
    bs_put_prices[k] <- bs_prices$put
  }
  
  list(call_prices = call_prices, put_prices = put_prices,
       bs_call_prices = bs_call_prices, bs_put_prices = bs_put_prices,
       call_se = call_se, put_se = put_se)
}

# Price options for all models
options_GBM <- price_options(terminal_prices_rn_GBM, strikes, d_f)
options_SVJD <- price_options(terminal_prices_rn_SVJD, strikes, d_f)
options_Heston <- price_options(terminal_prices_rn_Heston, strikes, d_f)


## ----7-recap-options-2, cache=FALSE, eval=TRUE------------------------------------------------
kableExtra::kable(as.data.frame(options_GBM))
kableExtra::kable(as.data.frame(options_Heston))
kableExtra::kable(as.data.frame(options_SVJD))


## ----3-fit-ARIMA-3, cache=FALSE, eval=TRUE----------------------------------------------------
# =============================================================================
# 4 - FIT ARIMA MODELS TO EXTRACT RISK PREMIUM
# =============================================================================

n_periods <- nrow(martingale_diff_GBM)
n_paths <- ncol(martingale_diff_GBM)

martingale_increments_GBM <- diff(martingale_diff_GBM)
martingale_increments_SVJD <- diff(martingale_diff_SVJD)
martingale_increments_Heston <- diff(martingale_diff_Heston)

# Initialize storage arrays
arima_residuals_GBM <- array(NA, dim = c(nrow(martingale_increments_GBM), n_paths))
centered_arima_residuals_GBM <- array(NA, dim = c(nrow(martingale_increments_GBM), n_paths))
means_arima_residuals_GBM <- rep(NA, n_paths)
arima_models_GBM <- list()

arima_residuals_SVJD <- array(NA, dim = c(nrow(martingale_increments_SVJD), n_paths))
centered_arima_residuals_SVJD <- array(NA, dim = c(nrow(martingale_increments_SVJD), n_paths))
means_arima_residuals_SVJD <- rep(NA, n_paths)
arima_models_SVJD <- list()

arima_residuals_Heston <- array(NA, dim = c(nrow(martingale_increments_Heston), n_paths))
centered_arima_residuals_Heston <- array(NA, dim = c(nrow(martingale_increments_Heston), n_paths))
means_arima_residuals_Heston <- rep(NA, n_paths)
arima_models_Heston <- list()

# Fit ARIMA models to GBM
cat("Fitting ARIMA models to", n_paths, "GBM paths...\n")
for (i in 1:n_paths) {
  y <- as.numeric(martingale_increments_GBM[, i])
  fit <- forecast::auto.arima(y)
  arima_models_GBM[[i]] <- fit
  
  res <- as.numeric(residuals(fit))
  arima_residuals_GBM[, i] <- res
  
  centre_arima_residuals <- scale(res, center = TRUE, scale = FALSE)
  means_arima_residuals_GBM[i] <- attr(centre_arima_residuals, "scaled:center")
  centered_arima_residuals_GBM[, i] <- centre_arima_residuals[,1]
}

# Fit ARIMA models to SVJD
cat("Fitting ARIMA models to", n_paths, "SVJD paths...\n")
for (i in 1:n_paths) {
  y <- as.numeric(martingale_increments_SVJD[, i])
  fit <- forecast::auto.arima(y)
  arima_models_SVJD[[i]] <- fit
  
  res <- as.numeric(residuals(fit))
  arima_residuals_SVJD[, i] <- res
  
  centre_arima_residuals <- scale(res, center = TRUE, scale = FALSE)
  means_arima_residuals_SVJD[i] <- attr(centre_arima_residuals, "scaled:center")
  centered_arima_residuals_SVJD[, i] <- centre_arima_residuals[,1]
}

# Fit ARIMA models to Heston
cat("Fitting ARIMA models to", n_paths, "Heston paths...\n")
for (i in 1:n_paths) {
  y <- as.numeric(martingale_increments_Heston[, i])
  fit <- forecast::thetaf(y)
  arima_models_Heston[[i]] <- fit
  
  res <- as.numeric(residuals(fit))
  arima_residuals_Heston[, i] <- res
  
  centre_arima_residuals <- scale(res, center = TRUE, scale = FALSE)
  means_arima_residuals_Heston[i] <- attr(centre_arima_residuals, "scaled:center")
  centered_arima_residuals_Heston[, i] <- centre_arima_residuals[,1]
}

# Box-Ljung tests
pvalues_GBM <- sapply(1:ncol(centered_arima_residuals_GBM), 
                  function(i) Box.test(centered_arima_residuals_GBM[,i])$p.value)
cat("\nBox-Ljung test p-values (GBM):\n")
cat("Mean p-value:", mean(pvalues_GBM), "\n")
cat("Proportion > 0.05:", mean(pvalues_GBM > 0.05), "\n")

pvalues_SVJD <- sapply(1:ncol(centered_arima_residuals_SVJD), 
                  function(i) Box.test(centered_arima_residuals_SVJD[,i])$p.value)
cat("\nBox-Ljung test p-values (SVJD):\n")
cat("Mean p-value:", mean(pvalues_SVJD), "\n")
cat("Proportion > 0.05:", mean(pvalues_SVJD > 0.05), "\n")

pvalues_Heston <- sapply(1:ncol(centered_arima_residuals_Heston), 
                  function(i) Box.test(centered_arima_residuals_Heston[,i])$p.value)
cat("\nBox-Ljung test p-values (Heston):\n")
cat("Mean p-value:", mean(pvalues_Heston), "\n")
cat("Proportion > 0.05:", mean(pvalues_Heston > 0.05), "\n")

par(mfrow=c(1,3))
hist(pvalues_GBM, breaks=20, col='lightgreen',
     main="Box-Ljung P-values (GBM)",
     xlab="P-value")
abline(v=0.05, col='red', lwd=2, lty=2)

hist(pvalues_SVJD, breaks=20, col='lightblue',
     main="Box-Ljung P-values (SVJD)",
     xlab="P-value")
abline(v=0.05, col='red', lwd=2, lty=2)

hist(pvalues_Heston, breaks=20, col='lightcoral',
     main="Box-Ljung P-values (Heston)",
     xlab="P-value")
abline(v=0.05, col='red', lwd=2, lty=2)
par(mfrow=c(1,1))


## ----4-generate-rn-paths-3, cache=FALSE, eval=TRUE--------------------------------------------
# =============================================================================
# 5 - GENERATE RISK-NEUTRAL PATHS
# =============================================================================

bootstrap_independent <- function(x, R = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- length(x)
  
  # Each column is an independent bootstrap sample
  replicate(R, sample(x, size = n, replace = TRUE))
}

cat("\n\nGenerating risk-neutral paths from ALL historical paths...\n")
n_sim_per_path <- 20  # Generate 10 paths per historical path
times <- seq(0, maturity, length.out = n_periods)
discount_factor <- exp(r * times)

# Storage for all risk-neutral paths
all_S_tilde_GBM <- list()
all_S_tilde_SVJD <- list()
all_S_tilde_Heston <- list()

# Generate GBM risk-neutral paths
for (i in 1:n_paths) {
  resampled_residuals <- bootstrap_independent(centered_arima_residuals_GBM[, i], 
                                               R = n_sim_per_path)
  fit <- arima_models_GBM[[i]]
  fitted_increments <- as.numeric(fitted(fit))
  
  discounted_path <- matrix(0, nrow = n_periods, ncol = n_sim_per_path)
  discounted_path[1, ] <- S0
  
  increments <- matrix(scale(fitted_increments, center = TRUE, 
                             scale = FALSE)[,1], nrow = n_periods - 1, ncol = n_sim_per_path) + 
                resampled_residuals[1:(n_periods - 1), ]
  
  discounted_path[-1, ] <- S0 + apply(increments, 2, cumsum)
  S_tilde_price <- discounted_path * discount_factor
  all_S_tilde_GBM[[i]] <- S_tilde_price
}

# Generate SVJD risk-neutral paths
for (i in 1:n_paths) {
  resampled_residuals <- bootstrap_independent(centered_arima_residuals_SVJD[, i], 
                                               R = n_sim_per_path)
  fit <- arima_models_SVJD[[i]]
  fitted_increments <- as.numeric(fitted(fit))
  
  discounted_path <- matrix(0, nrow = n_periods, ncol = n_sim_per_path)
  discounted_path[1, ] <- S0
  
  increments <- matrix(scale(fitted_increments, center = TRUE, 
                             scale = FALSE)[,1], nrow = n_periods - 1, ncol = n_sim_per_path) + 
                resampled_residuals[1:(n_periods - 1), ]
  
  discounted_path[-1, ] <- S0 + apply(increments, 2, cumsum)
  S_tilde_price <- discounted_path * discount_factor
  all_S_tilde_SVJD[[i]] <- S_tilde_price
}

# Generate Heston risk-neutral paths
for (i in 1:n_paths) {
  resampled_residuals <- bootstrap_independent(centered_arima_residuals_Heston[, i], 
                                               R = n_sim_per_path)
  fit <- arima_models_Heston[[i]]
  fitted_increments <- as.numeric(fitted(fit))
  
  discounted_path <- matrix(0, nrow = n_periods, ncol = n_sim_per_path)
  discounted_path[1, ] <- S0
  
  increments <- matrix(scale(fitted_increments, center = TRUE, 
                             scale = FALSE)[,1], nrow = n_periods - 1, ncol = n_sim_per_path) + 
                resampled_residuals[1:(n_periods - 1), ]
  
  discounted_path[-1, ] <- S0 + apply(increments, 2, cumsum)
  S_tilde_price <- discounted_path * discount_factor
  all_S_tilde_Heston[[i]] <- S_tilde_price
}

# Combine all paths
S_tilde_combined_GBM <- do.call(cbind, all_S_tilde_GBM)
cat("Total GBM risk-neutral paths generated:", ncol(S_tilde_combined_GBM), "\n")

S_tilde_combined_SVJD <- do.call(cbind, all_S_tilde_SVJD)
cat("Total SVJD risk-neutral paths generated:", ncol(S_tilde_combined_SVJD), "\n")

S_tilde_combined_Heston <- do.call(cbind, all_S_tilde_Heston)
cat("Total Heston risk-neutral paths generated:", ncol(S_tilde_combined_Heston), "\n\n")

# Convert to time series
S_tilde_ts_GBM <- ts(S_tilde_combined_GBM, start = start(sim_GBM), 
                     frequency = frequency(sim_GBM))
S_tilde_ts_SVJD <- ts(S_tilde_combined_SVJD, start = start(sim_SVJD), 
                      frequency = frequency(sim_SVJD))
S_tilde_ts_Heston <- ts(S_tilde_combined_Heston, start = start(sim_Heston), 
                        frequency = frequency(sim_Heston))

# Visualize risk-neutral paths
par(mfrow=c(1,3))
esgtoolkit::esgplotbands(S_tilde_ts_GBM, 
                        main="Risk-Neutral Paths - GBM")
esgtoolkit::esgplotbands(S_tilde_ts_SVJD, 
                        main="Risk-Neutral Paths - SVJD")
esgtoolkit::esgplotbands(S_tilde_ts_Heston, 
                        main="Risk-Neutral Paths - Heston")

# Sample plots
par(mfrow=c(1,3))
matplot(S_tilde_combined_GBM[, sample(ncol(S_tilde_combined_GBM), 200)], 
        type='l', col=rgb(0,0,1,0.1),
        main="Risk-Neutral Paths (GBM)",
        xlab="Time step", ylab="Stock Price")
lines(rowMeans(S_tilde_combined_GBM), col='red', lwd=3)
abline(h=S0, col='green', lwd=2, lty=2)

matplot(S_tilde_combined_SVJD[, sample(ncol(S_tilde_combined_SVJD), 200)], 
        type='l', col=rgb(0,0,1,0.1),
        main="Risk-Neutral Paths (SVJD)",
        xlab="Time step", ylab="Stock Price")
lines(rowMeans(S_tilde_combined_SVJD), col='red', lwd=3)
abline(h=S0, col='green', lwd=2, lty=2)

matplot(S_tilde_combined_Heston[, sample(ncol(S_tilde_combined_Heston), 200)], 
        type='l', col=rgb(0,0,1,0.1),
        main="Risk-Neutral Paths (Heston)",
        xlab="Time step", ylab="Stock Price")
lines(rowMeans(S_tilde_combined_Heston), col='red', lwd=3)
abline(h=S0, col='green', lwd=2, lty=2)
par(mfrow=c(1,1))


## ----5-rn-verif-3, cache=FALSE, eval=TRUE-----------------------------------------------------
# =============================================================================
# 6 - VERIFY RISK-NEUTRAL PROPERTY
# =============================================================================

cat("\n=== RISK-NEUTRAL VERIFICATION ===\n\n")

terminal_prices_rn_GBM <- S_tilde_combined_GBM[n_periods, ]
terminal_prices_rn_SVJD <- S_tilde_combined_SVJD[n_periods, ]
terminal_prices_rn_Heston <- S_tilde_combined_Heston[n_periods, ]
capitalized_stock_price <- S0 * exp(r * maturity)

cat("GBM Risk-Neutral Verification:\n")
cat("Expected terminal price (Q):", capitalized_stock_price, "\n")
cat("Empirical mean:", mean(terminal_prices_rn_GBM), "\n")
cat("Difference:", mean(terminal_prices_rn_GBM) - capitalized_stock_price, "\n")
print(t.test(terminal_prices_rn_GBM - capitalized_stock_price))

cat("\nSVJD Risk-Neutral Verification:\n")
cat("Expected terminal price (Q):", capitalized_stock_price, "\n")
cat("Empirical mean:", mean(terminal_prices_rn_SVJD), "\n")
cat("Difference:", mean(terminal_prices_rn_SVJD) - capitalized_stock_price, "\n")
print(t.test(terminal_prices_rn_SVJD - capitalized_stock_price))

cat("\nHeston Risk-Neutral Verification:\n")
cat("Expected terminal price (Q):", capitalized_stock_price, "\n")
cat("Empirical mean:", mean(terminal_prices_rn_Heston), "\n")
cat("Difference:", mean(terminal_prices_rn_Heston) - capitalized_stock_price, "\n")
print(t.test(terminal_prices_rn_Heston - capitalized_stock_price))

# Visualization comparison
par(mfrow=c(3, 2))
hist(terminal_prices_physical_GBM, breaks=30, col=rgb(1,0,0,0.5),
     main="Terminal Prices: Physical (GBM)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_physical_GBM), col='red', lwd=2)
abline(v=S0*exp(mu*maturity), col='blue', lwd=2, lty=2)

hist(terminal_prices_rn_GBM, breaks=30, col=rgb(0,0,1,0.5),
     main="Terminal Prices: Risk-Neutral (GBM)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_rn_GBM), col='blue', lwd=2)
abline(v=S0*exp(r*maturity), col='red', lwd=2, lty=2)

hist(terminal_prices_physical_SVJD, breaks=30, col=rgb(1,0,0,0.5),
     main="Terminal Prices: Physical (SVJD)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_physical_SVJD), col='red', lwd=2)
abline(v=S0*exp(mu*maturity), col='blue', lwd=2, lty=2)

hist(terminal_prices_rn_SVJD, breaks=30, col=rgb(0,0,1,0.5),
     main="Terminal Prices: Risk-Neutral (SVJD)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_rn_SVJD), col='blue', lwd=2)
abline(v=S0*exp(r*maturity), col='red', lwd=2, lty=2)

hist(terminal_prices_physical_Heston, breaks=30, col=rgb(1,0,0,0.5),
     main="Terminal Prices: Physical (Heston)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_physical_Heston), col='red', lwd=2)
abline(v=S0*exp(mu*maturity), col='blue', lwd=2, lty=2)

hist(terminal_prices_rn_Heston, breaks=30, col=rgb(0,0,1,0.5),
     main="Terminal Prices: Risk-Neutral (Heston)",
     xlab="Price", xlim=c(50, 300))
abline(v=mean(terminal_prices_rn_Heston), col='blue', lwd=2)
abline(v=S0*exp(r*maturity), col='red', lwd=2, lty=2)
par(mfrow=c(1,1))


## ----6-option-pricing-3, cache=FALSE, eval=TRUE-----------------------------------------------
# =============================================================================
# 7 - OPTION PRICING
# =============================================================================

cat("\n=== OPTION PRICING ===\n\n")

bs_price <- function(S, K, r, sigma, T, q = 0) {

  d1 <- (log(S / K) + (r - q + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)

  call <- S * exp(-q * T) * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  put  <- K * exp(-r * T) * pnorm(-d2) - S * exp(-q * T) * pnorm(-d1)

  list(call = call, put = put)
}


strikes <- seq(80, 160, by=10)
d_f <- exp(-r * maturity)

# Function to price options
price_options <- function(terminal_prices, strikes, discount_factor) {
  n_strikes <- length(strikes)
  call_prices <- numeric(n_strikes)
  bs_call_prices <- numeric(n_strikes)
  put_prices <- numeric(n_strikes)
  bs_put_prices <- numeric(n_strikes)
  call_se <- numeric(n_strikes)
  put_se <- numeric(n_strikes)
  
  for (k in 1:n_strikes) {
    K <- strikes[k]
    call_payoffs <- pmax(terminal_prices - K, 0)
    call_prices[k] <- mean(call_payoffs) * discount_factor
    call_se[k] <- sd(call_payoffs) / sqrt(length(call_payoffs)) * discount_factor
    
    put_payoffs <- pmax(K - terminal_prices, 0)
    put_prices[k] <- mean(put_payoffs) * discount_factor
    put_se[k] <- sd(put_payoffs) / sqrt(length(put_payoffs)) * discount_factor
    
    bs_prices <- bs_price(S0, K, r, sigma=sigma, T=5, q = 0)
    bs_call_prices[k] <- bs_prices$call
    bs_put_prices[k] <- bs_prices$put
  }
  
  list(call_prices = call_prices, put_prices = put_prices,
       bs_call_prices = bs_call_prices, bs_put_prices = bs_put_prices,
       call_se = call_se, put_se = put_se)
}

# Price options for all models
options_GBM <- price_options(terminal_prices_rn_GBM, strikes, d_f)
options_SVJD <- price_options(terminal_prices_rn_SVJD, strikes, d_f)
options_Heston <- price_options(terminal_prices_rn_Heston, strikes, d_f)


## ----7-recap-options-3, cache=FALSE, eval=TRUE------------------------------------------------
kableExtra::kable(as.data.frame(options_GBM))
kableExtra::kable(as.data.frame(options_Heston))
kableExtra::kable(as.data.frame(options_SVJD))

