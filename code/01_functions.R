# Useful functions that will be used over again

# Draw Gaussian zCDP noise
# Input: vector of statistics, sensitivity, rho
# Output: noise vector of the same dimension as the input vector
# To generate noisy answer add the vector and the noise vector
gaussian_mech_zCDP_vec <- function(vec, sensitivity, rho) {
  # Set sigma for the noise draw
  sigma <- sqrt((sensitivity ^ 2) / (2 * rho))
  
  # Draw independent noise for each entry in the input vector and return output
  return(sapply(vec, function(v)
    rnorm(1, 0, sigma)))
  
}

# Generate synthetic data of binary input data
# Input: Dataset with binary features
# Output: Synthetic dataset 
synthetic_data <- function(data, dp = FALSE, rho = NULL) {
  # Generate histogram of counts
  counts <-
    table(apply(data, 1, function(x)
      paste(x, collapse = " ")))
  
  # If dp = TRUE generate noisy counts
  
  if (dp) {
    # Draw noise from gaussian, distribute rho evenly
    noise <- gaussian_mech_zCDP_vec(counts, 1, rho / length(counts))
    
    # Add noise to counts
    counts <- counts + noise
    
    # Handle potential negative counts
    if (any(counts < 0)) {
      counts <- counts + abs(min(counts[counts < 0]))
    }
  }
  
  # Normalize counts to get probabilities
  prob <- counts / sum(counts)
  
  # Use probabilities to draw N times from a multinomial distribution
  synth <- t(rmultinom(n = nrow(data), 1, prob))
  
  # Format Histogram back to original format
  synth_df <-
    data.frame(do.call(rbind, stringr::str_split(colnames(synth)[synth %*% 1:length(counts)], " ")))
  synth_df <- data.frame(lapply(synth_df, as.numeric))
  colnames(synth_df) <- colnames(data)
  
  # Return synthetic data frame
  return(synth_df)
}

# Function to calculate logit regression on synthetic data for reuse in bootstrap
# Function is particular to logit experiments
# Input: Data, dp and rho
# Output: Regression coefficient
logit_function <- function(data, dp = F, rho = NULL) {
  # Generate synthetic copy of input data
  synth_df <- synthetic_data(data, dp = dp, rho = rho)
  
  # Run logit regression on synthetic copy
  synth_reg <-
    glm(y ~ x, data = synth_df, family = binomial(link = "logit"))
  
  # Pull out coefficient of interest from result object
  beta_hat <- synth_reg$coefficients[2]
  
  # Return coefficient
  return(beta_hat)
}

# Function to sample data for the logit experiments
# Input: betas = True coefficients, n_data = Number of observations
# Output: Data frame 
sample_logit_data <- function(betas, n_data) {
  # Sample independent variable from Bernoulli distribution with p = 0.8
  x <- rbinom(n_data, 1, 0.8)
  
  # Generate probabilities for y using inverse logit link function
  p <- 1 / (1 + exp(-cbind(1, x) %*% betas))
  
  # Sample y from Bernoulli based on p
  y <- rbinom(n_data, 1, p)
  
  # Put y and x together in a data.frame
  orig_data <- data.frame(y, x)
  
  # Return data.frame
  return(orig_data)
}


jackknife_logit <- function(data) {
  res <- NULL
  
  for(i in 1:nrow(data)){
    tmp_reg <-
      glm(y ~ x, data = data[-i,], family = binomial(link = "logit"))
    
    tmp_coef <- tmp_reg$coefficients[2]
    
    res <- c(res, tmp_coef)
  }
  return(res)
}

bc_ci <- function(bootstrap_distribution, initial_estimate, jackknife_distribution, alpha = 0.025) {
  
  z0 <- qnorm(sum(bootstrap_distribution < initial_estimate)/length(bootstrap_distribution))
  
  
  theta_hat_dot <- mean(jackknife_distribution)
  
  
  
  a_hat <- sum((theta_hat_dot - jackknife_distribution)^3) / 6*sum((theta_hat_dot - jackknife_distribution)^2)^(3/2)
  
  alpha1 <-
    pnorm(z0 + ((z0 + qnorm(alpha)) / (1 - a_hat * (z0 + qnorm(
      alpha
    )))))
  
  alpha2 <-
    pnorm(z0 + ((z0 + qnorm(1 - alpha)) / (1 - a_hat * (z0 + qnorm(
      1 - alpha
    )))))
  
  quantile(bootstrap_distribution, c(alpha1, alpha2))
  
}

# Function to plot Confidence Intervals of repeated experiments
# Calculates the empirical coverage 
plot_cis <-
  function(d,
           true_pop_mean,
           label = NULL,
           xlim = c(0, 1),
           xlab = "Quantity of Interest") {
    par(mar = c(5.1, 4.1, 5.1, 2.1))
    plot(
      1,
      1,
      ylim = c(1, nrow(d)),
      xlim = xlim,
      type = "n",
      bty = "n",
      yaxt = "n",
      ylab = "",
      xlab = xlab,
      main = ifelse(
        is.null(label),
        paste0(
          "95% Confidence Intervals \n Empirical Coverage: ",
          round(sum(
            apply(d, 1, function(x)
              true_pop_mean >= x[1] &
                true_pop_mean <= x[2])
          ) / nrow(d) * 100, 1),
          "%"
        ),
        paste0(
          label,
          "\n",
          "95% Confidence Intervals \n Empirical Coverage: ",
          round(sum(
            apply(d, 1, function(x)
              true_pop_mean >= x[1] &
                true_pop_mean <= x[2])
          ) / nrow(d) * 100, 1),
          "%"
          
        )
      )
    )
    segments(d[, 1], 1:nrow(d), d[, 2], col = viridis::viridis(3)[(apply(d, 1, function(x)
      true_pop_mean >= x[1] & true_pop_mean <= x[2]) * 1) + 1])
    abline(v = true_pop_mean,
           lty = "dashed",
           col = viridis::viridis(3)[3])
    
  }