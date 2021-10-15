# Code to run experiments

# Source other scripts (Assuming the folder with this file is the working directory)
source("00_setup.R")
source("01_functions.R")

# Set up parallel computing backend
# Detect number of available cores, leave one free
n_cores <- parallel::detectCores() - 1

# Create the cluster
my_cluster <- parallel::makeCluster(n_cores,
                                    type = "FORK")

# Register the cluster to be used by %dopar% for parallelization of for loops
doParallel::registerDoParallel(cl = my_cluster)

# Check if the cluster is registered
foreach::getDoParRegistered()

# Setup of true population

# Set true population parameters for the regression
betas <- c(0, 1.5)
n_data <- 2000

# Set rho parameter
rho <- 0.00337869

# Check achieved privacy in terms of delta and epsilon
delta <- 1 / (2 * n_data)

epsilon <- rho + 2 * sqrt(rho * log(1 / delta))

# Initialize empty objects to store results
betas_hat <- NULL
dp_betas_hat <- NULL
bootstrap_dists <- list()

# Set number of Bootstrap replications
n_boot <- 2000

# Set number of outer loop/Monte Carlo replications
n_outer_loop <- 250

for (i in 1:n_outer_loop) {
  # Sample original data from population parameters
  orig_data <- sample_logit_data(betas, n_data)
  
  # Calculate logit regression on sample
  orig_reg <-
    glm(y ~ x, data = orig_data, family = binomial(link = "logit"))
  
  # Store coefficient
  beta_hat <- orig_reg$coefficients[2]
  
  # Generate DP synthetic data
  dp_synth_data <- synthetic_data(orig_data, dp = T, rho = rho)
  
  # Run regression on DP synthetic data
  naive_dp_reg <-
    glm(y ~ x, data = dp_synth_data, family = binomial(link = "logit"))
  
  dp_beta_hat <- naive_dp_reg$coefficients[2]
  
  # Append beta_hat and dp_beta_hat to results vector
  betas_hat <- c(betas_hat, beta_hat)
  dp_betas_hat <- c(dp_betas_hat, dp_beta_hat)
  
  # Get bootstrap distribution based on the DP synthetic data
  # Including privacy noise
  bootstrap_distribution <- foreach(# Iterator i
    i = 1:n_boot,
    # Combine results to a vector
    .combine = 'c') %dopar% {
      # Use parallel backend
      logit_function(dp_synth_data, dp = T, rho = rho)
    }
  
  # Append bootstrap distribution to results list
  bootstrap_dists[[i]] <- bootstrap_distribution
  
  # Report what iteration was just finished
  cat("Iteration: ", i, "\n")
}


# Get percentile confidence intervals from bootstrapping distributions
dp_boot_cis <- t(sapply(bootstrap_dists, function(x) quantile(x, c(0.025, 0.975))))


# Generate a plot of the confidence intervals
par(mfrow = c(1, 1))
plot_cis(dp_boot_cis, true_pop_mean = 1.5, xlim = c(0, 3))

# Compare the standard deviation of p^hat, p^tilde and p^hattilde

# SD of p^hat
sd(betas_hat)

# SD of p^tilde
sd(dp_betas_hat, na.rm = T)

# Average SD of bootstrap distributions
mean(sapply(bootstrap_dists, sd))


# Ridge plot of bootstrap distributions


par(mfrow = c(1, 1))
plot(
  1,
  1,
  xlim = c(0, 3),
  ylim = c(0, 2*length(bootstrap_dists) + 3),
  type = "n",
  yaxt = "n",
  ylab = "",
  bty = "n",
  xlab = "Logit Coefficient"
)



for (i in 1:length(bootstrap_dists)) {
  den <- density(bootstrap_dists[[i]])
  den$y <- den$y + (2 * i - 1)
  
  lines(den, col = viridis::viridis(3, 0.7)[1])
}

# Include the true distribution as a reference
true_den <- density(dp_betas_hat)
true_den$y <- true_den$y + (2 * i)

lines(true_den, col = viridis::viridis(3, 0.7)[2])


# Question: Smaple size?
# Effect of sample size? (up for noisy, down for p_hat)
#
