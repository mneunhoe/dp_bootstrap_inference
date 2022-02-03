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


n_rep <- 100
res_P <- vector("list", length(Ps))
names(res_P) <- names(Ps)

epsilons <- c(0.1, 1, 10, Inf)
B <- 1000

n_data <- 1000

alpha <- 0.05
projection_step <- TRUE
for (dataset in names(Ps)) {
  res_list <- vector("list", length(epsilons))
  names(res_list) <- paste0(epsilons)
  
  for (epsilon in epsilons) {
    ci_list <- vector("list", n_rep)
    pivot_ci_list <- vector("list", n_rep)
    percentile_ci_list <- vector("list", n_rep)
    
    for (i in 1:n_rep) {
      P_hat <- draw_sample(Ps[[dataset]]$P, n = n_data)
      
      # Generate dp cdf
      private_cdf <-
        dpCDF(
          P_hat,
          Ps[[dataset]]$settings[1],
          Ps[[dataset]]$settings[2],
          epsilon = epsilon,
          granularity = Ps[[dataset]]$settings[3],
          cdp = TRUE,
          num_trials = 1
        )
      
      if (projection_step) {
        pava <- isotone::gpava(private_cdf[[1]][[2]], private_cdf[[1]][[1]])
        private_cdf[[1]][[1]] <- pava$x
      }
      
      boot_dist <-
        cdf_bootstrap(
          cdf = private_cdf[[1]],
          n_data = n_data,
          B = B,
          lower_bound = Ps[[dataset]]$settings[1],
          upper_bound = Ps[[dataset]]$settings[2],
          epsilon = epsilon,
          granularity = Ps[[dataset]]$settings[3],
          cdp = TRUE,
          projection_step = projection_step
        )
      

      # setting up the sub directory
      sub_dir <- paste0("../results/", dataset)
      
      # check if sub directory exists 
      if (!file.exists(sub_dir)){
        dir.create(sub_dir)
      }
      
      saveRDS(
        list(private_cdf = private_cdf,
             boot_dist = boot_dist),
        paste(
          "../results/",
          dataset,
          "/",
          "eps_",
          epsilon,
          "_rep_",
          i,
          ".RDS",
          sep = ""
        )
      )
      
      
      
      percentile_ci <- quantile(sapply(boot_dist, function(x)
        get_quantile(q = 0.5, cdf_x = private_cdf[[1]][[2]], cdf_y = x)), c(alpha/2, 1-alpha/2))
        
      theta_hat <- get_quantile(private_cdf[[1]], 0.5)
      
      pivot_ci <- theta_hat - quantile(sapply(boot_dist, function(x)
        get_quantile(q = 0.5, cdf_x = private_cdf[[1]][[2]], cdf_y = x)) - theta_hat, c(1-alpha/2, alpha/2))
      
      pivot_ci_list[[i]] <- pivot_ci
      percentile_ci_list[[i]] <- percentile_ci
      
      
      # ci_list[[i]] <- dp_ci(
      #   P_hat,
      #   epsilon = epsilon,
      #   B = 1000,
      #   lower_bound = Ps[[dataset]]$settings[1],
      #   upper_bound = Ps[[dataset]]$settings[2],
      #   granularity = Ps[[dataset]]$settings[3]
      # )
      cat(i, "\n")
    }
    res_list[[paste0(epsilon)]] <- list(pivot_ci_list, percentile_ci_list)
  }
  res_P[[paste0(dataset)]] <- res_list
}

saveRDS(res_P, "../results/median_experiments.RDS")

#res_P <- readRDS("../results/median_experiments.RDS")

par(mfrow = c(1, 2))
lapply(names(Ps), function(x) summarize_results(res_P[[x]], dataset = x, true_value = Ps[[x]]$true_value) )



### Get mean

cdf <- private_cdf[[1]]
probs <- c(cdf[[1]][1], diff(cdf[[1]]))  

sum(probs*cdf[[2]])
hist(samp)
mean(samp)
median(samp)

