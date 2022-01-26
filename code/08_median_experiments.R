source("06_cdf_functions.R")
library(isotone)
set.seed(100)


# Set number of data points
n_data <- 1000

# Number of bootstrap repetitions

B <- 1000

# Set hyperparameters for dp cdf and discretized cdf
upper_bound <- 4
lower_bound <- -4
granularity <- 0.01
cdp <- TRUE
epsilon <- 1








# Parameters for mixture distribution

weights <- c(0.5, 0.5)
comp_means <- c(-2, 2)
comp_sds <- c(0.5, 0.5)


# Quantity of interest

true_theta <- 0

# Corresponding quantile
q <- pmixnorm(true_theta, weights = weights, component_means = comp_means, component_sds = comp_sds)

# Set up empty objects to collect results
res_sampling_dist <- NULL
res_priv_sampling_dist <- NULL
boot_res <- list()

# Settings for loop
projection_step <- TRUE
bootstrap <- TRUE

# Experiment Loop

n_rep <- 10

for (i in 1:n_rep) {
  
  # Sample new P_hat
  samp <- rmixnorm(n = n_data, weights = weights, component_means = comp_means, component_sds = comp_sds)
  
  # Discretize
  discretized_cdf <-
    dpCDF(
      samp,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      epsilon = Inf,
      granularity = granularity,
      cdp = cdp,
      num_trials = 1
    )
  

  
  
  
  
  
  # return quantile
  
  res_sampling_dist <-
    c(res_sampling_dist, get_quantile(discretized_cdf[[1]], q) - true_theta)
  
  # Generate dp cdf
  private_cdf <-
    dpCDF(
      samp,
      lower_bound,
      upper_bound,
      epsilon = epsilon,
      granularity = granularity,
      cdp = cdp,
      num_trials = 1
    )
  
  
                        
  
  # plot(res[[1]][[2]], res[[1]][[1]], type = "n", ylab = "Cumulative Probability", xlab = "x", bty = "n", las = 1, main = paste0("Epsilon: ", epsilon))
  # lapply(res, function(x) lines(x[[2]], x[[1]], col = viridis::viridis(3, 1)[2]))
  
  if (projection_step) {
    pava <- isotone::gpava(private_cdf[[1]][[2]], private_cdf[[1]][[1]])
    private_cdf[[1]][[1]] <- pava$x
  }
  
  
  
  
  theta_hat <-  get_quantile(private_cdf[[1]], q)
  
  res_priv_sampling_dist <-
    c(res_priv_sampling_dist, theta_hat - true_theta)
  
  
  if (bootstrap) {
    boot_dist <-
      cdf_bootstrap(
        cdf = private_cdf[[1]],
        n_data = n_data,
        B = B,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        epsilon = epsilon,
        granularity = granularity,
        cdp = cdp,
        projection_step = projection_step
      )
    
    
    
    boot_res[[i]] <-
      sapply(boot_dist, function(x) get_quantile(q = q, cdf_x = private_cdf[[1]][[2]], cdf_y = x) - theta_hat)
  }
  
  
  cat(i, "\n")
}

# Look at results
par(mfrow = c(1, 2))


hist(
  res_sampling_dist,
  xlim = c(-4, 4),
  border = "white",
  col = viridis::viridis(2)[1],
  main = "Non-private Sampling Distribution"
)
hist(
  res_priv_sampling_dist,
  xlim = c(-4, 4),
  border = "white",
  col = viridis::viridis(2)[2],
  main = "Sampling + Privacy Distribution"
)

par(mfrow = c(3, 4))

lapply(boot_res, function(x)
  hist(x, xlim = c(-4, 4)))


