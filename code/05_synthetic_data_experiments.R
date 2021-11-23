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
eps <- "1"
# Sample original data from population parameters
for(i in 1:100){
orig_data <- sample_logit_data(betas, n_data)


dir.create(paste0("~/Desktop/github/private-pgm/data/epsilon_",eps,"/exp",i),recursive = T)
write.csv(orig_data, paste0("~/Desktop/github/private-pgm/data/epsilon_",eps,"/exp",i,"/logit.csv"), row.names = F)
}


for(i in 1:100){
 syn_df <-  read.csv(paste0("~/Desktop/github/private-pgm/data/epsilon_",eps,"/exp",i,"/synth_data.csv"))
 
 for(j in 1:100){
 samp <- sample(nrow(syn_df), nrow(syn_df), replace = T)
 
 write.csv(syn_df[samp,], paste0("~/Desktop/github/private-pgm/data/epsilon_",eps,"/exp",i,"/boot",j,".csv"), row.names = F)
 
 }
}

orig_cis <- NULL
naive_cis <- NULL
boot_dists <- vector(mode ="list")
eps <- "1"
for(j in 1:100){
  orig_data <- read.csv(paste0("~/Desktop/github/private-pgm/data/epsilon_",eps,"/exp",j,"/logit.csv"))
  orig_reg <-
    glm(y ~ x, data = orig_data, family = binomial(link = "logit"))

  orig_cis <- rbind(orig_cis, confint(orig_reg)[2,])
  
  synth_data <- read.csv(paste0("~/Desktop/github/private-pgm/data/epsilon_",eps,"/exp",j,"/synth_data.csv"))
  synth_reg <-
    glm(y ~ x, data = synth_data, family = binomial(link = "logit"))
  summary(synth_reg)
  #naive_cis <- rbind(naive_cis, confint(synth_reg)[2,])
  naive_cis <- rbind(naive_cis, tryCatch(confint(synth_reg)[2,], error= function(cond){ c(0, 0)}))
  
  
  res <- coef(synth_reg)[2]
  
  for(i in 1:100){

synth_data <- read.csv(paste0("~/Desktop/github/private-pgm/data/epsilon_",eps,"/exp",j,"/synth_boot",i,".csv"))
synth_reg <-
  glm(y ~ x, data = synth_data, family = binomial(link = "logit"))
res <- rbind(res, coef(synth_reg)[2])

  }
#   
 boot_dists[[j]] <- res

}


par(mfrow = c(1, 3))
plot_cis(d = orig_cis, true_pop_mean = 1.5, xlim = c(0, 3), label = "Original Data")
plot_cis(d = naive_cis, true_pop_mean = 1.5, xlim = c(0, 3), label = "Naive CI on Synthetic Data")
plot_cis(d = t(sapply(boot_dists, quantile, c(0.025, 0.975), na.rm = T)), true_pop_mean = 1.5, xlim = c(0, 3), label = "Bootstrap")


quantile(res[,2], c(0.025, 0.975))





orig_reg <-
  glm(y ~ x, data = orig_data, family = binomial(link = "logit"))

synth_reg <-
  glm(y ~ x, data = synth_data, family = binomial(link = "logit"))

