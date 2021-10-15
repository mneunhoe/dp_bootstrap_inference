
# Collect needed libraries in a vector
p_needed <-
  c("knitr", # Package to generate pdf/html output from Rmarkdown
    "MASS", # R port of Venables and Ripley, ``Modern Applied Statistics with S'' 
    "viridis", # Accessible color palettes
    "zeallot", #  Allows python style assignment of lists to multiple objects
    "foreach", # Package to parallelize for loops
    "parallel", # Parallelization backend
    "doParallel", # dopar
    "here" # Package for easy file paths
    )

# Get all installed packages
packages <- rownames(installed.packages())

# Check which packages of p_needed are not in packages
p_to_install <- p_needed[!(p_needed %in% packages)]

# If at least one package is not installed, install the new packages
if (length(p_to_install) > 0) {
  install.packages(p_to_install)
}

# Load packages to environment and return whether this was succesful
sapply(p_needed, require, character.only = TRUE)


