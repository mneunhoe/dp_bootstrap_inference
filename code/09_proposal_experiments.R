source("06_cdf_functions.R")
library(isotone)



adult <-
  read.csv(
    "https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data",
    header = FALSE
  )
adult_age_full <- adult[, 1]

Ps <- list(
  "mixture" = list(
    P = function(n) {
      rmixnorm(
        n,
        weights = c(0.5, 0.5),
        component_means = c(-2, 2),
        component_sds = c(0.5, 0.5)
      )
    },
    settings = c(-4, 4, 0.01),
    true_value = 0
  ),
  "normal" = list(
    P = rnorm,
    settings = c(-4, 4, 0.01),
    true_value = 0
  ),
  "adult" = list(
    P = adult_age_full,
    settings = c(18, 99, 0.1),
    true_value = median(adult_age_full)
  )
)


n_rep <- 2
res_P <- vector("list", length(Ps))
names(res_P) <- names(Ps)

epsilons <- c(0.1, 1, 10, Inf)
epsilon <- 1
dataset <- "mixture"
for (dataset in names(Ps)) {
  res_list <- vector("list", length(epsilons))
  names(res_list) <- paste0(epsilons)
  
  for (epsilon in epsilons) {
    ci_list <- vector("list", n_rep)
    
    for (i in 1:n_rep) {
      P_hat <- draw_sample(Ps[[dataset]]$P, n = 1000)
      
      ci_list[[i]] <- dp_ci(
        P_hat,
        epsilon = epsilon,
        B = 100,
        lower_bound = Ps[[dataset]]$settings[1],
        upper_bound = Ps[[dataset]]$settings[2],
        granularity = Ps[[dataset]]$settings[3]
      )
      cat(i, "\n")
    }
    res_list[[paste0(epsilon)]] <- ci_list
  }
  res_P[[paste0(dataset)]] <- res_list
}

saveRDS(res_P, "../results/median_experiments.RDS")

par(mfrow = c(1, 2))
lapply(names(Ps), function(x) summarize_results(res_P[[x]], dataset = x, true_value = Ps[[x]]$true_value) )



### Get mean

cdf <- private_cdf[[1]]
probs <- c(cdf[[1]][1], diff(cdf[[1]]))  

sum(probs*cdf[[2]])
hist(samp)
mean(samp)
median(samp)
