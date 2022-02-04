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
  "lognormal" = list(
    P = rlnorm,
    settings = c(0, 15, 0.01),
    true_value = 1
  ),
  "adult" = list(
    P = adult_age_full,
    settings = c(18, 99, 0.05),
    true_value = median(adult_age_full)
  )
)

epsilons <- c(Inf)

n_rep <- 100
res_P <- vector("list", length(Ps))
names(res_P) <- names(Ps)

n_data <- 1000

alpha <- 0.05
projection_step <- TRUE
boot <- F

for (dataset in names(Ps)) {
  res_list <- vector("list", length(epsilons))
  names(res_list) <- paste0(epsilons)
  
  for (epsilon in epsilons) {
    exact_ci_list <- vector("list", n_rep)
    
    
    for (i in 1:n_rep) {
      P_hat <- draw_sample(Ps[[dataset]]$P, n = n_data)
      
      exact_ci_list[[i]] <- cimed(P_hat)
      
      cat(i, "\n")
    }
  
      res_list[[paste0(epsilon)]] <- list(exact_ci_list)
    
  }
  res_P[[paste0(dataset)]] <- res_list
}
