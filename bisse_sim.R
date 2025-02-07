# BiSSE Analysis Script
# Written by Miles Woodcock-Girard
# Igic Lab, University of Illinois at Chicago

library(diversitree)
library(dplyr)
library(ggplot2)
library(purrr)

set.seed(42)

# ==============================================================================
#                        DEFINE BiSSE SIMULATION PARAMETERS
# ==============================================================================
N <- 100                            # Number of simulations per parameter set
PARAM_FILE_PATH <- "params.csv"     # Path to parameters CSV file, formatted:
#
#    lambda0  lambda1   mu0   mu1   q01   q10
# 1      0.1      0.1  0.03  0.03  0.01  0.01
# 2      ...
# 3      ...
# ...

# Where:

#   lambda0 is birth/speciation rate for trait 0,
#   lambda1 is birth/speciation rate for trait 1,
#   mu0 is death/extinction rate for trait 0,
#   mu1 is death/extinction rate for trait 1,
#   q01 is rate of transition from 0 to 1,
#   q10 is rate of transition from 1 to 0

# Get parameters from CSV file
param_set <- read.csv(PARAM_FILE_PATH,
                      header = TRUE,
                      sep = ",")


# Define BiSSE simulation function
# Returns true and estimated likelihoods and parameters
simulate_bisse <- function(parameters, i) {
  # Repeat failed simulations
  repeat {
    # Run BiSSE simulation, obtain tree end states
    phy <- tree.bisse(parameters, max.t = 50, x0 = 0)
    end_states <- phy$tip.state
    
    #print(end_states)
    # Perform failed simulation checks
    
    # Check if tree too small
    if (length(end_states) <= 2) {
      message(sprintf("  Retrying simulation %d: tree has fewer than 2 tips.", i))
      next
    }
    
    # Check for NA values or invalid states
    if (any(is.na(end_states)) || !all(end_states %in% c(0, 1))) {
      message(sprintf("  Retrying simulation %d: Unexpected end states detected.", i))
      next
    }
    
    # Ensure both character states (0 and 1) are present in final tree
    if (length(unique(end_states)) < 2) {
      message(sprintf("  Retrying simulation %d: Only one state survived. Estimation impossible.", i))
      next
    }

    # Get BiSSE likelihood function  
    lik <- make.bisse(phy, end_states)
  
    # Log-likelihood of true parameters
    log_lik <- lik(parameters)
  
    # Estimate parameters via MLE
    pt <- starting.point.bisse(phy)
    
    fit <- tryCatch({
      find.mle(lik, pt, method = "subplex")
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(fit)) {
      warning(sprintf("  Retrying simulation %d: MLE failure.", i))
      next
    }
    
    message(sprintf("Simulations completed: %d / %d | Num. tips: %d | 0s: %d | 1s: %d",
                    i, N, length(phy$tip.state),
                    length(phy$tip.state[phy$tip.state == 0]),
                    length(phy$tip.state[phy$tip.state == 1])))
    
    # Log-likelihood of MLE
    log_lik_fit <- logLik(fit)
  
    # Store results as a named list
    return(list(tree = phy,
                parameters = parameters,
                log_lik_true = log_lik,
                estimated_params = coef(fit),
                log_lik_fit = log_lik_fit))
  }
}



# ==============================================================================
#                   MAIN BiSSE SIMULATION, ESTIMATION LOOP
# ==============================================================================
# Use `purrr::map()` to loop over parameter sets, returning list of results
# dataframes, one per parameter set.
all_results <- pmap(param_set, function(lambda0, lambda1, mu0, mu1, q01, q10) {
    # Prepare BiSSE simulation parameters
    parameters <- c(lambda0, lambda1, mu0, mu1, q01, q10)
    
    # Run N simulations for current parameter set, get results
    results <- map(1:N, ~ simulate_bisse(parameters, .))
    
    # Convert results to dataframe list
    results_df <- do.call(rbind.data.frame, lapply(results, function(res) {
      c(res$log_lik_true, res$log_lik_fit, res$parameters, res$estimated_params)
    }))

    # Rename columns in results dataframe list
    colnames(results_df) <- c("Log-likelihood (true)", "Log-likelihood (fit)",
                              "lambda0 (true)", "lambda1 (true)",
                              "mu0 (true)", "mu1 (true)",
                              "q01 (true)", "q10 (true)",
                              "lambda0 (est)", "lambda1 (est)",
                              "mu0 (est)", "mu1 (est)",
                              "q01 (est)", "q10 (est)")
    
    return(results_df)
  })

# NOT WORKING YET. PERFORM

# Convert dataframe list into singular dataframe
# all_results_df <- bind_rows(all_results)

# ==============================================================================
#               ANALYZE BiSSE SIMULATION, ESTIMATION RESULTS
# ==============================================================================


# # Calculate the mean of estimated parameters
# mean_estimates <- lapply(all_results, function(res) {
#   c(mean(res["lambda0 (est)"]), mean(res["lambda1 (est)"]),
#     mean(res["mu0 (est)"]), mean(res["mu1 (est)"]),
#     mean(res["q01 (est)"]), mean(res["q10 (est)"]))
# })
#   
#  all_results %>%
#  summarise(across(starts_with("lambda"), mean),
#            across(starts_with("mu"), mean),
#            across(starts_with("q"), mean))

# Calculate bias as the difference between true and estimated parameters
# bias <- all_results_df %>%
#   summarise(across(starts_with("lambda"), ~ mean(.x) - parameters[1]),
#             across(starts_with("mu"), ~ mean(.x) - parameters[3]),
#             across(starts_with("q"), ~ mean(.x) - parameters[5]))
# 
# # Calculate the standard deviation of estimates
# precision <- all_results_df %>%
#   summarise(across(starts_with("lambda"), sd),
#             across(starts_with("mu"), sd),
#             across(starts_with("q"), sd))