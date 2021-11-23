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

# Set true population parameters for the regression
betas <- c(0, 1.5)
p_x <- 0.8

n_data <- 2000

# Set rho parameter
rho <- 0.00337869
#rho <- 1
# Check achieved privacy in terms of delta and epsilon
delta <- 1 / (2 * n_data)

epsilon <- rho + 2 * sqrt(rho * log(1 / delta))
dp_res <- NULL
res <- NULL
boot_dists <- list()
dp_boot_dists <- list()
startval <- c(0, 0)

B <- 1000

logit_link <- function(x) {
  1 / (1 + exp(-x))
}

loglik <- function(theta, X, y, counts) {
  
  
  
  y_hat <- logit_link(cbind(1,X)%*%theta)
  
  y_hat0 <- logit_link(cbind(1, 0)%*%theta)
  y_hat1 <- logit_link(cbind(1, 1)%*%theta)
  
  # Negative Loglikelihood from histogram
  if(!is.null(counts)){
    ll <- as.numeric(-1/sum(counts)*(counts[3]*log(y_hat0) + counts[1]*log(1-y_hat0) + counts[4]*log(y_hat1) + counts[2]*log(1-y_hat1)))
  } else {
    ll <- -mean(log(y_hat * y + (1 - y_hat) * (1 - y)))
  }
  #ll <- -prod(((y_hat^ y) * ((1 - y_hat) ^ (1 - y))))
  
  #ll <- -log(sum(y_hat * y + (1 - y_hat) * (1 - y)))
  return(ll)
  
}

ll_fct <- function(counts, dp = F){
  
  if(dp){
    
    tmp_counts <- c(counts[3]+counts[4], counts[2]+counts[4])
    
    tmp_counts <-  tmp_counts + gaussian_mech_zCDP_vec(tmp_counts, 1, rho/length(tmp_counts))
    
    dp_data <- data.frame(y = rbinom(n_data, 1, tmp_counts[1]/n_data), x = rbinom(n_data, 1, tmp_counts[2]/n_data))
    
    tmp_counts <-  table(apply(dp_data, 1, function(x)
      paste(x, collapse = " ")))
    
    
    #tmp_counts <- counts + gaussian_mech_zCDP_vec(counts, 1, rho/length(counts))
    
    if (any(tmp_counts < 0)) {
      tmp_counts <- tmp_counts + abs(min(tmp_counts[tmp_counts < 0]))
    }
    
  } else {
    tmp_counts <- counts
  }
  
  samp_counts <- rowSums(rmultinom(sum(tmp_counts), 1, tmp_counts))
  
  optim(startval, loglik, X = orig_data$x, y = orig_data$y, counts = samp_counts)$par[2]
  
}

for(i in 1:100){
  
  x <- rbinom(n_data, 1, p_x)
  
  p <- 1 / (1 + exp(-cbind(1, x) %*% betas))
  
  y <- rbinom(n_data, 1, p)
  
  orig_data <- data.frame(y, x)
  
  counts <-
    table(apply(orig_data, 1, function(x)
      paste(x, collapse = " ")))
  
  b_hat <- optim(startval, loglik, X = orig_data$x, y = orig_data$y, counts = counts)$par[2]
  res <- c(res, b_hat)
  
  boot_dist <- foreach(
    i = 1:B,
    .combine = 'c'
  ) %dopar% {
    ll_fct(
      counts, dp = F)
  }
  
  boot_dists[[i]] <- boot_dist
  
  
  
  dp_counts <- colSums(orig_data) + gaussian_mech_zCDP_vec(colSums(orig_data), 1, rho/length(colSums(orig_data)))
  
  dp_data <- data.frame(y = rbinom(n_data, 1, dp_counts[1]/n_data), x = rbinom(n_data, 1, dp_counts[2]/n_data))
  
  dp_counts <-  table(apply(dp_data, 1, function(x)
    paste(x, collapse = " ")))
  
  dp_boot_dist <- foreach(
    i = 1:B,
    .combine = 'c'
  ) %dopar% {
    ll_fct(
      dp_counts, dp = T)
  }
  
  
  #dp_boot_dist <- apply(t(replicate(1000, rowSums(rmultinom(sum(tmp_counts), 1, tmp_counts)))), 1, function(x) optim(startval, loglik, X = orig_data$x, y = orig_data$y, counts = x)$par[2])
  dp_boot_dists[[i]] <- dp_boot_dist
  
  dp_res <- c(dp_res, optim(startval, loglik, X = orig_data$x, y = orig_data$y, counts = dp_counts)$par[2])
  
  cat("Iter: ", i, "\n")
}

pdf(paste0("~/Downloads/logit_cis_n_",n_data,"_epsilon_",round(epsilon,2),".pdf"),
    width = 11,
    height = 8.5)

par(mfrow = c(1, 2))
plot_cis(t(sapply(boot_dists, quantile, c(0.025, 0.975))), true_pop_mean = 1.5, xlim = c(-2, 5.5), label = "Bootstrap", xlab = "Logit Coefficient")
plot_cis(t(sapply(dp_boot_dists, quantile, c(0.025, 0.975))), true_pop_mean = 1.5, xlim = c(-2, 5.5), label = "DP Bootstrap", xlab = "Logit Coefficient")

dev.off()
