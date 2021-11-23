# Distribution of P and Private P

library(zeallot)

library(foreach)



gaussian_mech_zCDP_vec <- function(vec, sensitivity, rho) {
  sigma <- sqrt((sensitivity ^ 2) / (2 * rho))
  
  sapply(vec, function(v)
    rnorm(1, 0, sigma))
  
}

synthetic_data <- function(data, dp = F, boot = F, simplex = F) {
  counts <-
    table(apply(data, 1, function(x)
      paste(x, collapse = " ")))
  
  prob_m <- counts / nrow(data)
  
  if (dp) {
    
    if(simplex){
      
      noise <- rejection_sample(tol = 1e-2)
      
      counts <- counts + noise[1,]
      
    } else {
      noise <- gaussian_mech_zCDP_vec(counts, 1, rho / 4)
      
      counts <- counts + noise
      # for(i in 1:length(counts)) {
      #   counts[i] <- counts[i] + noise[i]
      #   counts[-i] <- counts[-i] - noise[i]/(length(counts)-1)
      # }
    }
    prob_m <- counts / sum(counts)
    
    if (any(prob_m < 0)) {
      tmp <- prob_m + abs(min(prob_m[prob_m < 0]))
      #tmp <- tmp - abs(max(prob_m[prob_m>1]))
      prob_m <- tmp / sum(tmp)
    }
    
    prob_m <- prob_m/sum(prob_m)
  }
  
  if (dp & boot) {
    synth <- t(rmultinom(n = nrow(data), 1, prob_m))
    #synth <- t(rmultinom(n = min(c(nrow(data), floor(sum(counts)))), 1, prob_m))
  } else {
    synth <- t(rmultinom(n = nrow(data), 1, prob_m))
    #synth <- t(rmultinom(n = min(c(nrow(data), floor(sum(counts)))), 1, prob_m))
  }
  synth_df <-
    data.frame(do.call(rbind, stringr::str_split(colnames(synth)[synth %*% 1:4], " ")))
  synth_df <- data.frame(lapply(synth_df, as.numeric))
  colnames(synth_df) <- colnames(data)
  return(list(synth_df, noise, counts))
}





logit_function <- function(data, dp = F, boot = F, simplex = F) {
  tmp <- synthetic_data(data, dp = dp, boot = boot, simplex = simplex)[[1]]
  
  naive_dp_reg <-
    glm(y ~ x, data = tmp, family = binomial(link = "logit"))
  
  dp_beta_hat <- naive_dp_reg$coefficients[2]
  
  return(dp_beta_hat)
}

betas <- c(0, 1.5)

rho <- 0.00337869

delta <- 1 / (2 * n_data)

epsilon <- rho + 2 * sqrt(rho * log(1 / delta))

n.cores <- parallel::detectCores() - 1

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)

#check cluster definition (optional)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()


betas_hat <- NULL
dp_betas_hat <- NULL

bootstrap_dists <- list()

n_data <- 6000
n_boot <- 2000

for(i in 1:100){
x <- rbinom(n_data, 1, 0.8)

p <- 1 / (1 + exp(-cbind(1, x) %*% betas))

y <- rbinom(n_data, 1, p)

orig_data <- data.frame(y, x)


orig_reg <-
  glm(y ~ x, data = orig_data, family = binomial(link = "logit"))

beta_hat <- orig_reg$coefficients[2]



c(dp_synth_data, noise_draw, noisy_histogram) %<-% synthetic_data(orig_data, dp = T)


naive_dp_reg <-
  glm(y ~ x, data = dp_synth_data, family = binomial(link = "logit"))

dp_beta_hat <- naive_dp_reg$coefficients[2]

betas_hat <- c(betas_hat, beta_hat)
dp_betas_hat <- c(dp_betas_hat, dp_beta_hat)


bootstrap_distribution <- foreach(
  i = 1:n_boot,
  .combine = 'c'
) %dopar% {
  logit_function(
    dp_synth_data, dp = T, boot = T, simplex = F
  )
}


bootstrap_dists[[i]] <- bootstrap_distribution

cat("Iteration: ", i, "\n")
}

par(mfrow = c(1, 2))
hist(betas_hat, breaks = 200, border = NA)
hist(dp_betas_hat, breaks = 200, border = NA)


mean(betas_hat)
mean(dp_betas_hat, na.rm = T)

sd(betas_hat)
sd(dp_betas_hat, na.rm = T)


mean(sapply(bootstrap_dists, sd))



plot(density(bootstrap_dists[[1]]), col = "grey")


par(mfrow = c(1, 1))
plot(1, 1, xlim = c(0.5, 2.1), ylim = c(0, 203), type = "n", yaxt = "n", ylab = "", bty = "n", xlab = "Logit Coefficient")



for(i in 1:100){
  
  den <- density(bootstrap_dists[[i]])
  den$y <- den$y+(2*i-1)
  
  
lines(den, col = viridis::viridis(3, 0.7)[1])

# Fill area
#polygon(den, col = "slateblue1")
# Private P
}

true_den <- density(dp_betas_hat)
true_den$y <- true_den$y+(2*i)

lines(true_den, col = viridis::viridis(3, 0.7)[2])

hist(bootstrap_dists[[1]])




replicate(2, logit_function(
  dp_synth_data, dp = T, boot = T, simplex = F
))



# Question: Smaple size?
# Effect of sample size? (up for noisy, down for p_hat)
#

orig_data
counts <-
  table(apply(orig_data, 1, function(x)
    paste(x, collapse = " ")))

logit_link <- function(x) {
  1 / (1 + exp(-x))
}

betas

theta <- c(0.1, 0.1)

X <- orig_data$x
y <- orig_data$y

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

n_data <- 1000
dp_res <- NULL
res <- NULL
boot_dists <- list()
dp_boot_dists <- list()


B <- 1000

ll_fct <- function(counts, dp = F){
  
  if(dp){
  tmp_counts <- counts + gaussian_mech_zCDP_vec(counts, 1, rho/length(counts))
  
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
  
  x <- rbinom(n_data, 1, 0.8)
  
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
  
  dp_counts <- counts + gaussian_mech_zCDP_vec(counts, 1, rho/length(counts))
  
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

t(sapply(boot_dists, quantile, c(0.025, 0.975)))

pdf(paste0("~/Downloads/logit_cis_n_",n_data,".pdf"),
    width = 11,
    height = 8.5)

par(mfrow = c(1, 2))
plot_cis(t(sapply(boot_dists, quantile, c(0.025, 0.975))), true_pop_mean = 1.5, xlim = c(0, 3.5), label = "Bootstrap", xlab = "Logit Coefficient")
plot_cis(t(sapply(dp_boot_dists, quantile, c(0.025, 0.975))), true_pop_mean = 1.5, xlim = c(0, 3.5), label = "DP Bootstrap", xlab = "Logit Coefficient")

dev.off()

optim(startval, loglik, X = orig_data$x, y = orig_data$y, counts = counts)$par[2]

mean(dp_res)
plot(density(dp_res))
lines(density(res), col = "red")
mean(res)
mean(dp_res)


library(zeallot)

synthetic_data <- function(data, dp = F, boot = F, simplex = F) {
  counts <-
    table(apply(data, 1, function(x)
      paste(x, collapse = " ")))
  
  prob_m <- counts / nrow(data)
  
  if (dp) {
    
    if(simplex){
      
      noise <- rejection_sample(tol = 1e-2)
      
      counts <- counts + noise[1,]
      
    } else {
      noise <- gaussian_mech_zCDP_vec(counts, 1, rho / 4)
      
      counts <- counts + noise
      # for(i in 1:length(counts)) {
      #   counts[i] <- counts[i] + noise[i]
      #   counts[-i] <- counts[-i] - noise[i]/(length(counts)-1)
      # }
    }
    prob_m <- counts / sum(counts)
    
    if (any(prob_m < 0)) {
      tmp <- prob_m + abs(min(prob_m[prob_m < 0]))
      #tmp <- tmp - abs(max(prob_m[prob_m>1]))
      prob_m <- tmp / sum(tmp)
    }
    
    prob_m <- prob_m/sum(prob_m)
  }
  
  if (dp & boot) {
    synth <- t(rmultinom(n = nrow(data), 1, prob_m))
    #synth <- t(rmultinom(n = min(c(nrow(data), floor(sum(counts)))), 1, prob_m))
  } else {
    synth <- t(rmultinom(n = nrow(data), 1, prob_m))
    #synth <- t(rmultinom(n = min(c(nrow(data), floor(sum(counts)))), 1, prob_m))
  }
  synth_df <-
    data.frame(do.call(rbind, stringr::str_split(colnames(synth)[synth %*% 1:4], " ")))
  synth_df <- data.frame(lapply(synth_df, as.numeric))
  colnames(synth_df) <- colnames(data)
  return(list(synth_df, noise, counts))
}



logit_function <- function(data, dp = F, boot = F, simplex = F) {
  tmp <- synthetic_data(data, dp = dp, boot = boot, simplex = simplex)[[1]]
  
  naive_dp_reg <-
    glm(y ~ x, data = tmp, family = binomial(link = "logit"))
  
  dp_beta_hat <- naive_dp_reg$coefficients[2]
  
  return(dp_beta_hat)
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


cis <- NULL
boot_cis <- NULL
naive_dp_cis <- NULL
boot_dp_cis <- NULL
recentered_dp_cis <- NULL
bc_dp_cis <- NULL

noise_df <- NULL

n_data <- 2000

for (i in 1:100) {
  x <- rbinom(n_data, 1, 0.8)
  
  p <- 1 / (1 + exp(-cbind(1, x) %*% betas))
  
  y <- rbinom(n_data, 1, p)
  
  orig_data <- data.frame(y, x)
  
  
  orig_reg <-
    glm(y ~ x, data = orig_data, family = binomial(link = "logit"))
  
  beta_hat <- orig_reg$coefficients[2]
  
  
  ci <-
    beta_hat + c(qnorm(0.025), qnorm(0.975)) * summary(orig_reg)$coefficients[2, 2]
  
  cis <- rbind(cis, ci)
  
  
  c(dp_synth_data, noise_draw, noisy_histogram) %<-% synthetic_data(orig_data, dp = T)
  
  
  naive_dp_reg <-
    glm(y ~ x, data = dp_synth_data, family = binomial(link = "logit"))
  
  dp_beta_hat <- naive_dp_reg$coefficients[2]
  
  
  naive_dp_ci <-
    dp_beta_hat + c(qnorm(0.025), qnorm(0.975)) * summary(naive_dp_reg)$coefficients[2, 2]
  
  naive_dp_cis <- rbind(naive_dp_cis, naive_dp_ci)
  
  boot_dp_replicates <- replicate(2000, logit_function(
    dp_synth_data, dp = T, boot = T, simplex = F
  ))
  
  boot_dp_ci <-
    quantile(boot_dp_replicates, c(0.025, 0.975))
  
  jackknife_distribution <- jackknife_logit(dp_synth_data)
  
  bc_dp_ci <- bc_ci(boot_dp_replicates, dp_beta_hat, jackknife_distribution)
  
  bc_dp_cis <- rbind(bc_dp_cis, bc_dp_ci)
  
  recentered_dp_ci <- boot_dp_ci + (dp_beta_hat - mean(boot_dp_replicates))
  
  boot_dp_cis <- rbind(boot_dp_cis, boot_dp_ci)
  recentered_dp_cis <- rbind(recentered_dp_cis, recentered_dp_ci)
  
  boot_ci <-
    quantile(replicate(100, logit_function(orig_data, dp = F)), c(0.025, 0.975))
  
  boot_cis <- rbind(boot_cis, boot_ci)
  noise_df <- rbind(noise_df, noise_draw)
  
  cat("Iteration: ", i, "\n")
}




jackknife_distribution <- jackknife_logit(dp_synth_data)


bootstrap_distribution <- boot_dp_replicates
initial_estimate <- dp_beta_hat

alpha <- 0.025





bc_ci(bootstrap_distribution, initial_estimate, jackknife_distribution)

pdf("~/Downloads/logit_cis_11.pdf",
    width = 11,
    height = 8.5)

par(mfrow = c(1, 5))
plot_cis(
  d = cis,
  true_pop_mean = 1.5,
  label = "Normal Approx.",
  xlim = c(-2, 5),
  xlab = "Logit Coefficient"
)

plot_cis(
  d = boot_cis,
  true_pop_mean = 1.5,
  label = "Bootstrapped",
  xlim = c(-2, 5),
  xlab = "Logit Coefficient"
)

plot_cis(
  d = naive_dp_cis,
  true_pop_mean = 1.5,
  label = "Naive DP",
  xlim = c(-2, 5),
  xlab = "Logit Coefficient"
)

plot_cis(
  d = boot_dp_cis,
  true_pop_mean = 1.5,
  label = "Bootstrapped DP",
  xlim = c(-2, 5),
  xlab = "Logit Coefficient"
)
plot_cis(
  d = bc_dp_cis,
  true_pop_mean = 1.5,
  label = "BC Bootstrapped DP",
  xlim = c(-2, 5),
  xlab = "Logit Coefficient"
)
dev.off()



plot_cis(
  d = boot_dp_cis,
  true_pop_mean = 1.5,
  label = "Bootstrapped DP",
  xlim = c(-2, 5),
  xlab = "Logit Coefficient"
)


apply(naive_dp_cis, 1, function(x) x[2] - x[1])

plot(apply(boot_dp_cis, 1, function(x) x[2] - x[1]), round(apply(noise_df, 1, function(x) abs(sum(x)))))
plot(apply(boot_dp_cis, 1, function(x) x[2] - x[1]), round(apply(noise_df, 1, var)))


cor(apply(boot_dp_cis, 1, function(x) x[2] - x[1]), round(apply(noise_df, 1, sum)))

text(rep(-1.5, 100), 1:100, label = round(apply(noise_df, 1, sum), 1), cex = 0.4)
text(rep(-1, 100), 1:100, label = round(apply(noise_df, 1, var), 1), cex = 0.4)
par(mfrow = c(1, 1))
plot_cis(
  d = cis,
  true_pop_mean = 1.5,
  label = "Normal Approx.",
  xlim = c(-2, 5)
)
plot_cis(
  d = naive_dp_cis,
  true_pop_mean = 1.5,
  label = "Naive DP",
  xlim = c(-2, 5)
)

plot_cis(
  d = boot_cis,
  true_pop_mean = 1.5,
  label = "Bootstrapped",
  xlim = c(-2, 5)
)
plot_cis(
  d = boot_dp_cis,
  true_pop_mean = 1.5,
  label = "Bootstrapped DP",
  xlim = c(-2, 5)
)

nrow(boot_cis)




def gaussian_mech_zCDP_vec(vec, sensitivity, rho):sigma = np.sqrt((sensitivity **
                                                                     2) / (2 * rho))
return [v + np.random.normal(loc = 0, scale = sigma) for v in vec]


gaussian_mech_zCDP_vec <- function(vec, sensitivity, rho) {
  sigma <- sqrt((sensitivity ^ 2) / (2 * rho))
  
  sapply(vec, function(v)
    rnorm(1, 0, sigma))
  
}


prob_m

quantile(replicate(1000, logit_function(dp_synth_data, dp = T)), c(0.025, 0.975))


counts <-
  table(apply(data, 1, function(x)
    paste(x, collapse = " ")))

# rho/m

rho <- 0.0337869

delta <- 1 / (2 * n_data)

epsilon <- rho + 2 * sqrt(rho * log(1 / delta))


a <- gaussian_mech_zCDP_vec(counts, 1, rho / 4)



counts + noise




sum(counts)

sapply(1:length(noise), function(x)
  - noise[x] / (length(noise) - 1))

counts[1] + noise[1]

counts[-1] - noise[1] / 3


synth_D <- cbind()


confint(glm(y ~ x, data = orig_data, family = binomial(link = "logit")))

replicate(500, synthetic_data(dp_synth_data), simplify = F)


rlaplace <- function(n, b = 1) {
  rexp(n, rate = 1 / b) - rexp(n, rate = 1 / b)
}

binom_lapl <- function(n, size, p, b) {
  p <- (p * n + rlaplace(1, b)) / n
  p <- ifelse(p > 1, 1, p)
  p <- ifelse(p < 0, 0, p)
  rbinom(n, size, p)
  
}

n_data <- 500
sensitivity <- 1
epsilon <- 1 / 3
epsilon <- 0.1
b <- sensitivity / epsilon

true_pop_mean <- 0.8

cis <- NULL
naive_dp_cis <- NULL

replicates <- list()

dp_replicates <- list()

dp_l_replicates <- list()


for (i in 1:100) {
  data_set <- rbinom(n_data, 1, true_pop_mean)
  
  p_hat <- mean(data_set)
  
  ci <-
    p_hat + c(-1.96, 1.96) * sqrt((p_hat * (1 - p_hat)) / length(data_set))
  
  cis <- rbind(cis, ci)
  
  dp_p_hat <- (sum(data_set) + rlaplace(1, b)) / length(data_set)
  
  dp_p_hat <- ifelse(dp_p_hat > 1, 1, dp_p_hat)
  dp_p_hat <- ifelse(dp_p_hat < 0, 0, dp_p_hat)
  
  naive_dp_ci <-
    dp_p_hat + c(-1.96, 1.96) * sqrt((dp_p_hat * (1 - dp_p_hat)) / length(data_set))
  
  naive_dp_cis <- rbind(naive_dp_cis, naive_dp_ci)
  
  replicates[[i]] <- replicate(1000, rbinom(n_data, 1, p_hat))
  
  
  dp_replicates[[i]] <- replicate(1000, rbinom(n_data, 1, dp_p_hat))
  
  dp_l_replicates[[i]] <-
    replicate(1000, binom_lapl(n_data, 1, dp_p_hat, b))
}


boot_cis <-
  t(sapply(replicates, function(x)
    quantile(colMeans(x), c(0.025, 0.975))))
dp_boot_cis <-
  t(sapply(dp_replicates, function(x)
    quantile(colMeans(x), c(0.025, 0.975))))
dp_l_boot_cis <-
  t(sapply(dp_l_replicates, function(x)
    quantile(colMeans(x), c(0.025, 0.975))))


sum(apply(cis, 1, function(x)
  true_pop_mean > x[1] & true_pop_mean < x[2]))
sum(apply(boot_cis, 1, function(x)
  true_pop_mean > x[1] & true_pop_mean < x[2]))

sum(apply(naive_dp_cis, 1, function(x)
  true_pop_mean > x[1] & true_pop_mean < x[2]))
sum(apply(dp_boot_cis, 1, function(x)
  true_pop_mean > x[1] & true_pop_mean < x[2]))
sum(apply(dp_l_boot_cis, 1, function(x)
  true_pop_mean > x[1] & true_pop_mean < x[2]))

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

pdf("~/Downloads/proportion_cis.pdf",
    width = 11,
    height = 8.5)
par(mfrow = c(1, 4))
plot_cis(
  d = cis,
  true_pop_mean = true_pop_mean,
  label = "Normal Approx.",
  xlim = c(0, 1)
)
plot_cis(
  d = boot_cis,
  true_pop_mean = true_pop_mean,
  label = "Bootstrapped",
  xlim = c(0, 1)
)
# plot_cis(
#   d = naive_dp_cis,
#   true_pop_mean = true_pop_mean,
#   label = "Naive DP, Normal Approx.",
#   xlim = c(0, 1)
# )
plot_cis(
  d = dp_boot_cis,
  true_pop_mean = true_pop_mean,
  label = "Naive DP",
  xlim = c(0, 1)
)
plot_cis(
  d = dp_l_boot_cis,
  true_pop_mean = true_pop_mean,
  label = "DP, Adj. Bootstrap",
  xlim = c(0, 1)
)
dev.off()




apply(dp_l_boot_cis, 1, function(x)
  true_pop_mean >= x[1] & true_pop_mean <= x[2])



## Multivariate

betas <- c(0, 1.5)

x <- rbinom(n_data, 1, 0.8)

p <- 1 / (1 + exp(-cbind(1, x) %*% betas))

y <- rbinom(n_data, 1, p)

orig_data <- data.frame(y, x)

paste(c(0, 1)[1], c(0, 1)[2])

counts <-
  table(apply(cbind(y, x), 1, function(x)
    paste(x[1], x[2])))

prob_m <- counts / n_data

prob_m_dp <- (counts + rlaplace(length(counts), b)) / n_data

synth <- t(rmultinom(n = n_data, 1, prob_m))
synth_df <-
  data.frame(do.call(rbind, stringr::str_split(colnames(synth)[synth %*% 1:4], " ")))
synth_df <- data.frame(lapply(synth_df, as.numeric))
colnames(synth_df) <- c("y", "x")

data <- orig_data
test <- rejection_sample(tol = 1e-2)

counts + test[1,]
synthetic_data <- function(data, dp = F, boot = F, simplex = F) {
  counts <-
    table(apply(data, 1, function(x)
      paste(x, collapse = " ")))
  
  prob_m <- counts / nrow(data)
  
  if (dp) {
    
    if(simplex){
      
      noise <- rejection_sample(tol = 1e-2)
      
      counts <- counts + noise[1,]
      
    } else {
      noise <- gaussian_mech_zCDP_vec(counts, 1, rho / 4)
      
      counts <- counts + noise
      # for(i in 1:length(counts)) {
      #   counts[i] <- counts[i] + noise[i]
      #   counts[-i] <- counts[-i] - noise[i]/(length(counts)-1)
      # }
    }
    prob_m <- counts / sum(counts)
    
    if (any(prob_m < 0)) {
      tmp <- prob_m + abs(min(prob_m[prob_m < 0]))
      #tmp <- tmp - abs(max(prob_m[prob_m>1]))
      prob_m <- tmp / sum(tmp)
    }
    
    prob_m <- prob_m/sum(prob_m)
  }
  
  if (dp & boot) {
    #synth <- t(rmultinom(n = nrow(data), 1, prob_m))
    synth <- t(rmultinom(n = min(c(nrow(data)-3, floor(sum(counts)))), 1, prob_m))
  } else {
    #synth <- t(rmultinom(n = nrow(data), 1, prob_m))
    synth <- t(rmultinom(n = min(c(nrow(data)-3, floor(sum(counts)))), 1, prob_m))
  }
  synth_df <-
    data.frame(do.call(rbind, stringr::str_split(colnames(synth)[synth %*% 1:4], " ")))
  synth_df <- data.frame(lapply(synth_df, as.numeric))
  colnames(synth_df) <- colnames(data)
  return(synth_df)
}

logit_function <- function(data, dp = F, boot = F, simplex = F) {
  tmp <- synthetic_data(data, dp = dp, boot = boot, simplex = simplex)
  
  naive_dp_reg <-
    glm(y ~ x, data = tmp, family = binomial(link = "logit"))
  
  dp_beta_hat <- naive_dp_reg$coefficients[2]
  
  return(dp_beta_hat)
}



cis <- NULL
boot_cis <- NULL
naive_dp_cis <- NULL
boot_dp_cis <- NULL
recentered_dp_cis <- NULL


for (i in 1:100) {
  x <- rbinom(n_data, 1, 0.8)
  
  p <- 1 / (1 + exp(-cbind(1, x) %*% betas))
  
  y <- rbinom(n_data, 1, p)
  
  orig_data <- data.frame(y, x)
  
  
  orig_reg <-
    glm(y ~ x, data = orig_data, family = binomial(link = "logit"))
  
  beta_hat <- orig_reg$coefficients[2]
  
  
  ci <-
    beta_hat + c(qnorm(0.025), qnorm(0.975)) * summary(orig_reg)$coefficients[2, 2]
  
  cis <- rbind(cis, ci)
  
  
  dp_synth_data <- synthetic_data(orig_data, dp = T)
  
  naive_dp_reg <-
    glm(y ~ x, data = dp_synth_data, family = binomial(link = "logit"))
  
  dp_beta_hat <- naive_dp_reg$coefficients[2]
  
  
  naive_dp_ci <-
    dp_beta_hat + c(qnorm(0.025), qnorm(0.975)) * summary(naive_dp_reg)$coefficients[2, 2]
  
  naive_dp_cis <- rbind(naive_dp_cis, naive_dp_ci)
  
  boot_dp_replicates <- replicate(1000, logit_function(
    dp_synth_data, dp = T, boot = T, simplex = F
  ))
  
  boot_dp_ci <-
    quantile(boot_dp_replicates, c(0.025, 0.975))
  
  recentered_dp_ci <- boot_dp_ci + (dp_beta_hat - mean(boot_dp_replicates))
  
  boot_dp_cis <- rbind(boot_dp_cis, boot_dp_ci)
  recentered_dp_cis <- rbind(recentered_dp_cis, recentered_dp_ci)
  
  boot_ci <-
    quantile(replicate(1000, logit_function(orig_data, dp = F)), c(0.025, 0.975))
  
  boot_cis <- rbind(boot_cis, boot_ci)
  
  cat("Iteration: ", i, "\n")
}


pdf("~/Downloads/logit_cis_5.pdf",
    width = 11,
    height = 8.5)
par(mfrow = c(1, 5))
plot_cis(
  d = cis,
  true_pop_mean = 1.5,
  label = "Normal Approx.",
  xlim = c(-2, 5),
  xlab = "Logit Coefficient"
)

plot_cis(
  d = boot_cis,
  true_pop_mean = 1.5,
  label = "Bootstrapped",
  xlim = c(-2, 5),
  xlab = "Logit Coefficient"
)

plot_cis(
  d = naive_dp_cis,
  true_pop_mean = 1.5,
  label = "Naive DP",
  xlim = c(-2, 5),
  xlab = "Logit Coefficient"
)

plot_cis(
  d = boot_dp_cis,
  true_pop_mean = 1.5,
  label = "Bootstrapped DP",
  xlim = c(-2, 5),
  xlab = "Logit Coefficient"
)
plot_cis(
  d = recentered_dp_cis,
  true_pop_mean = 1.5,
  label = "Recentered Bootstrapped DP",
  xlim = c(-2, 5),
  xlab = "Logit Coefficient"
)
dev.off()


par(mfrow = c(1, 1))
plot_cis(
  d = cis,
  true_pop_mean = 1.5,
  label = "Normal Approx.",
  xlim = c(-2, 5)
)
plot_cis(
  d = naive_dp_cis,
  true_pop_mean = 1.5,
  label = "Naive DP",
  xlim = c(-2, 5)
)

plot_cis(
  d = boot_cis,
  true_pop_mean = 1.5,
  label = "Bootstrapped",
  xlim = c(-2, 5)
)
plot_cis(
  d = boot_dp_cis,
  true_pop_mean = 1.5,
  label = "Bootstrapped DP",
  xlim = c(-2, 5)
)

nrow(boot_cis)




def gaussian_mech_zCDP_vec(vec, sensitivity, rho):sigma = np.sqrt((sensitivity **
                                                                     2) / (2 * rho))
return [v + np.random.normal(loc = 0, scale = sigma) for v in vec]


gaussian_mech_zCDP_vec <- function(vec, sensitivity, rho) {
  sigma <- sqrt((sensitivity ^ 2) / (2 * rho))
  
  sapply(vec, function(v)
    rnorm(1, 0, sigma))
  
}


prob_m

quantile(replicate(1000, logit_function(dp_synth_data, dp = T)), c(0.025, 0.975))


counts <-
  table(apply(data, 1, function(x)
    paste(x, collapse = " ")))

# rho/m

rho <- 0.0337869

delta <- 1 / (2 * n_data)

epsilon <- rho + 2 * sqrt(rho * log(1 / delta))


a <- gaussian_mech_zCDP_vec(counts, 1, rho / 4)



counts + noise




sum(counts)

sapply(1:length(noise), function(x)
  - noise[x] / (length(noise) - 1))

counts[1] + noise[1]

counts[-1] - noise[1] / 3


synth_D <- cbind()


confint(glm(y ~ x, data = orig_data, family = binomial(link = "logit")))

replicate(500, synthetic_data(dp_synth_data), simplify = F)

