source("06_cdf_functions.R")
set.seed(100)

epsilon_settings <- c(0.1, 0.5, 1, 10, Inf)
granularity_settings <- c(0.1, 0.01, 0.005)
bound_settings <- c(4, 5)
n_settings <- c(100, 500, 1000)
settings <- expand.grid(epsilon_settings, granularity_settings, bound_settings, n_settings)


for(setting in 1:nrow(settings)){
# Set number of data points
n_data <- settings[setting, 4]

# Set hyperparameters for dp cdf
upper_bound <- settings[setting, 3]
lower_bound <- -settings[setting, 3]
granularity <- settings[setting, 2]
cdp <- TRUE
epsilon <- settings[setting, 1]

# Number of bootstrap repetitions

B <- 1000

# Set number of repetitions of experiment
n_rep <- 100

#Set up empty list for results
res_list <- vector("list", n_rep)


for(rep in 1:n_rep){
x <- rnorm(n_data)

if(epsilon == Inf) {
  lower_bound <- min(x)
  upper_bound <- max(x)
}

res <- dpCDF(x, lower_bound, upper_bound, epsilon, granularity, cdp, num_trials = 1)

# True CDF for standard normal data
if(rep == 1) cu <- curve(pnorm, lower_bound, upper_bound, n = length(res[[1]][[1]]), add = T, col = "black")
#par(mfrow = c(2, 2))
# plot(res[[1]][[2]], res[[1]][[1]], type = "n", ylab = "Cumulative Probability", xlab = "x", bty = "n", las = 1, main = paste0("Epsilon: ", epsilon))
# lapply(res, function(x) lines(x[[2]], x[[1]], col = viridis::viridis(3, 1)[1]))




boot_dist <- vector("list", B)
for(i in 1:B){
#x_boot <- sample(x, length(x), replace = T)
u <- runif(length(x), min = min(res[[1]][[1]]), max = max(res[[1]][[1]]))
x_boot <- res[[1]][[2]][sapply(u, function(x) which.min(abs(x  - res[[1]][[1]])))]
boot_cdf <- dpCDF(x_boot, lower_bound, upper_bound, epsilon, granularity, cdp, num_trials = 1)
boot_dist[[i]] <- boot_cdf[[1]][[1]]
#lines(boot_cdf[[1]][[2]], boot_cdf[[1]][[1]], col = viridis::viridis(3, 0.3)[2])
}






#curve(pnorm, -4, 4, n = length(res[[1]][[1]]), add = T, col = viridis::viridis(2)[2])

boot_mat <- do.call(cbind, boot_dist)

boot_cis <- t(apply(boot_mat, 1, quantile, c(0.025, 0.975)))

# lines(res[[1]][[2]], boot_cis[,1], lty = "dashed", col = "red")
# lines(res[[1]][[2]], boot_cis[,2], lty = "dashed", col = "red")


boot_cis_sim <- t(apply(boot_mat, 1, quantile, c((1-0.95^(1/length(res[[1]][[1]])))/2, 1-(1-0.95^(1/length(res[[1]][[1]])))/2)))

# lines(res[[1]][[2]], boot_cis_sim[,1], lty = "dashed", col = "orange")
# lines(res[[1]][[2]], boot_cis_sim[,2], lty = "dashed", col = "orange")
# 
# lapply(res, function(x) lines(x[[2]], x[[1]], col = viridis::viridis(3, 1)[1]))


#cu$y


coverage_points <- sum(cu$y >= boot_cis[,1] & cu$y <= boot_cis[,2])/nrow(boot_cis)
coverage_sim <- all((cu$y >= boot_cis_sim[,1] & cu$y <= boot_cis_sim[,2]))

# First entry will always be 0 and last 1, so let's discard it for the coverage calculation
coverage_sim_adj <- all((cu$y >= boot_cis_sim[,1] & cu$y <= boot_cis_sim[,2])[2:(nrow(boot_cis_sim)-1)])

res_list[[rep]] <- list(coverage_points, coverage_sim, coverage_sim_adj)
cat("Repetition: ", rep, "\n")
}

saveRDS(res_list, paste0("../results/", n_data, "_", epsilon,"_",granularity, "_",upper_bound, ".RDS"))
}

n_data <- 1000
epsilon <- 0.5
granularity <- 0.005
upper_bound <- 5
res_list <- readRDS( paste0("../results/", n_data, "_", epsilon,"_",granularity, "_",upper_bound, ".RDS"))


mean(sapply(res_list, function(x) (x[[1]])), na.rm = T)
mean(sapply(res_list, function(x) (x[[2]])))
mean(sapply(res_list, function(x) (x[[3]])))

0.948537^257

lines(res[[1]][[2]], cbind(apply(boot_mat, 1, mean) + apply(sweep(boot_mat, 1, apply(boot_mat, 1, mean)), 1, quantile, c(0.025, 0.975))[1,], apply(boot_mat, 1, mean) + apply(sweep(boot_mat, 1, apply(boot_mat, 1, mean)), 1, quantile, c(0.025, 0.975))[2,])
[,1], lty = "dashed", col = "gold")
lines(res[[1]][[2]], cbind(apply(boot_mat, 1, mean) + apply(sweep(boot_mat, 1, apply(boot_mat, 1, mean)), 1, quantile, c(0.025, 0.975))[1,], apply(boot_mat, 1, mean) + apply(sweep(boot_mat, 1, apply(boot_mat, 1, mean)), 1, quantile, c(0.025, 0.975))[2,])
[,2], lty = "dashed", col = "gold")


lines(res[[1]][[2]], cbind(res[[1]][[1]] + apply(sweep(boot_mat, 1, res[[1]][[1]]), 1, quantile, c(0.025, 0.975))[1,], res[[1]][[1]] + apply(sweep(boot_mat, 1, res[[1]][[1]]), 1, quantile, c(0.025, 0.975))[2,])
      [,1], lty = "dashed", col = "blue")
lines(res[[1]][[2]], cbind(res[[1]][[1]] + apply(sweep(boot_mat, 1, res[[1]][[1]]), 1, quantile, c(0.025, 0.975))[1,], res[[1]][[1]] + apply(sweep(boot_mat, 1, res[[1]][[1]]), 1, quantile, c(0.025, 0.975))[2,])
      [,2], lty = "dashed", col = "blue")



x^513=0.95



0.9999^513
boot_dist



0.99995^1000




sum(rbinom(100000, 1, 0.95))

x2 <- x
x2[x2<lower] <- NA
x2[x2>upper] <- NA
lines(ecdf(x2), col = viridis::viridis(2)[2])



which.min(abs(u[1]  - res[[1]][[1]]))



dpCDF(x_tilde, lower_bound, upper_bound, epsilon, granularity, cdp, num_trials = 1)

par(mfrow=c(1, 2))
hist(x, breaks = 10)
hist(x_tilde, breaks = 20)

res[[1]][[2]][which.min(abs(u[1]  - res[[1]][[1]]))]

lines(seq(lower, upper, length.out = length(res[[1]])), res[[1]], col = viridis::viridis(3)[2])



points(-4, u[1], pch = 19)


res[[1]]

res[[1]]
i <- 2
tree[[1]][[1]]
child <- 
wBelows <- wBelow(tree)


counts <- countBelow(tree, wBelows)

child <- tree[[3]][[1]]




list(list(1000, c(-4, 4)), lapply(bins, function(b) dpHistogram(x, 1000, -4, 4, 0.1, b)))


breaks <- seq(lower_bound, upper_bound, length.out = bins+1)
true_hist <- .Call(graphics:::C_BinCount, x, breaks, T, T)




