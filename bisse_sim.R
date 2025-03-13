# BiSSE Analysis Script
# Written by Miles Woodcock-Girard
# Igic Lab, University of Illinois at Chicago

library(diversitree)
library(dplyr)
library(ggplot2)
library(purrr)
library(R.utils)

args <- commandArgs(trailingOnly=TRUE)

# Set random seed for simulation
#set.seed(1)

# Set error rate to apply to tips
#error_rate <- 0.05
error_rate_str <- args[1]
error_rate <- as.double(error_rate_str)

# ==============================================================================
#                        DEFINE BiSSE SIMULATION PARAMETERS
# ==============================================================================
N <- 100                           # Number of simulations per parameter set
PARAM_FILE_PATH <- "params.csv"    # Path to parameters CSV file, formatted:
#
#    lambda0  lambda1   mu0   mu1   q01   q10
# 1      0.1      0.1  0.03  0.03  0.01  0.01
# 2      ...
# 3      ...
# 
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
simulate_bisse <- function(parameters, error_rate, i) {
  # Repeat failed simulations
  repeat {
    # Run BiSSE simulation, obtain tree end states
    phy <- tree.bisse(parameters, max.t = 5, x0 = 0)
    end_states <- phy$tip.state
        
    # Check for single-tip tree
    if (length(end_states) < 2) {
      cat(sprintf("  Retrying simulation %d: tree has fewer than 2 tips.\n", i))
      next
    }
    
    # Check for NA values or invalid states
    if (any(is.na(end_states)) || !all(end_states %in% c(0, 1))) {
      cat(sprintf("  Retrying simulation %d: Unexpected end states detected.\n", i))
      next
    }
    
    # Ensure both character states are in final tree
    if (length(unique(end_states)) < 2) {
      cat(sprintf("  Retrying simulation %d: Only one state survived. Estimation impossible.\n", i))
      next
    }

    cat(sprintf("  Tree built for sim: %d | Num. tips: %d (0s: %d | 1s: %d)\n",
		    i, length(phy$tip.state),
                    length(phy$tip.state[phy$tip.state == 0]),
                    length(phy$tip.state[phy$tip.state == 1])))

    # Perform MLE with true data
    cat("  Starting MLE ... ")
    lik_true <- make.bisse(phy, end_states)
    log_lik_true <- lik_true(parameters)
    fit_true <- tryCatch(withTimeout(find.mle(lik_true, starting.point.bisse(phy),
                                              method = "subplex"),
                                     timeout = 60, onTimeout = "error"),
	                 error = function(e) NULL)
    if (is.null(fit_true)) {
      cat(sprintf("  Retrying simulation %d: MLE failure to converge\n", i))
      #fit_true <- find.mle(lik_true, starting.point.bisse(phy))
      next
    }
    log_lik_true <- logLik(fit_true)
    cat("done.\n")

    # Apply random error to tip states
    # TODO: Explore bias in flipping values with error
    repeat {
      end_states_error <- end_states
      # Determine number of errors to add to tree
      num_errors <- round(length(end_states) * error_rate)
      
      # Ensure at least one error in end states
      if (num_errors == 0) {
        num_errors <- 1
      }
      
      # Ensure adding errors does not fix tree, preventing infinite loop
      if (num_errors >= length(end_states[end_states == 0]) && num_errors >= length(end_states[end_states == 1])) {
        num_errors <- num_errors - 1
      }
      error_indices <- sample(seq_along(end_states_error), num_errors)
      end_states_error[error_indices] <- ifelse(end_states_error[error_indices] == 0, 1, 0)
     
      # Ensure new end_states with error added contains both character states
      if (length(unique(end_states_error)) == 2) {
        cat(sprintf("    Added %d errors.\n", num_errors))
        break
      }
    }
    
    # Perform MLE with noisy data
    cat("  Starting MLE ... ")
    lik_error <- make.bisse(phy, end_states_error)
    log_lik_error <- lik_error(parameters)
    fit_error <- tryCatch(withTimeout(find.mle(lik_error, starting.point.bisse(phy),
	         	              method = "subplex"), timeout = 60, onTimeout = "error"),
		                      error = function(e) NULL)
    if (is.null(fit_error)) {
      cat(sprintf("  Retrying simulation %d: MLE failure to converge.\n", i))
      #fit_error <- find.mle(lik_error, starting.point.bisse(phy))
      next
    }
    log_like_error <- logLik(fit_error)
    cat("done.\n")
    
    cat(sprintf("Simulations completed: %d / %d | Num. tips: %d (0s: %d | 1s: %d) | Num. errors: %d\n",
                    i, N, length(phy$tip.state),
                    length(phy$tip.state[phy$tip.state == 0]),
                    length(phy$tip.state[phy$tip.state == 1]),
                    num_errors))
    
    h <- history.from.sim.discrete(phy, 0:1)
    plot(h, phy)
  
    # Store results as a named list
    return(list(tree = phy,
                parameters = parameters,
                log_lik_true = log_lik_true,
                estimated_params_true = coef(fit_true),
                log_lik_error = log_lik_error,
                estimated_params_error = coef(fit_error),
                num_errors = num_errors))
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

    cat("=====================================================\n")
    cat("=====================================================\n")
    cat("  NOW USING PARAMETER SET:\n")
    cat(sprintf("    λ0: %f\n", lambda0))
    cat(sprintf("    λ1: %f\n", lambda1))
    cat(sprintf("    μ0: %f\n", mu0))
    cat(sprintf("    μ1: %f\n", mu1))
    cat(sprintf("   q01: %f\n", q01))
    cat(sprintf("   q10: %f\n\n\n", q10))
    
    # Run N simulations for current parameter set, get results
    results <- map(1:N, ~ simulate_bisse(parameters, error_rate, .))
    
    # Convert results to dataframe list
    results_df <- do.call(rbind.data.frame, lapply(results, function(res) {
      c(res$log_lik_true, res$log_lik_error, res$parameters,
        res$estimated_params_true, res$estimated_params_error,
        res$num_errors, length(res$tree$tip.state))
    }))

    # Rename columns in results dataframe list
    colnames(results_df) <- c("Log-likelihood (true)", "Log-likelihood (error)",
                              "lambda0 (true)", "lambda1 (true)",
                              "mu0 (true)", "mu1 (true)",
                              "q01 (true)", "q10 (true)",
                              "lambda0 (est, true)", "lambda1 (est, true)",
                              "mu0 (est, true)", "mu1 (est, true)",
                              "q01 (est, true)", "q10 (est, true)",
                              "lambda0 (est, error)", "lambda1 (est, error)",
                              "mu0 (est, error)", "mu1 (est, error)",
                              "q01 (est, error)", "q10 (est, error)",
                              "Num errors", "Num Tips")
    
    return(results_df)
  })

# Convert dataframe list into singular dataframe
all_results_df <- bind_rows(all_results)

# ==============================================================================
#               ANALYZE BiSSE SIMULATION, ESTIMATION RESULTS
# ==============================================================================


# Group estimation results of simulations by parameter sets
df_grouped <- all_results_df %>% group_by(`lambda0 (true)`, `lambda1 (true)`,
                                          `mu0 (true)`, `mu1 (true)`,
                                          `q01 (true)`, `q10 (true)`)

# Get mean parameter estimates before and after error introduced
mean_estimates <- df_grouped %>% summarise(mean_lambda0_est = mean(`lambda0 (est, true)`),
                                           mean_lambda1_est = mean(`lambda1 (est, true)`),
                                           mean_mu0_est = mean(`mu0 (est, true)`),
                                           mean_mu1_est = mean(`mu1 (est, true)`),
                                           mean_q01_est = mean(`q01 (est, true)`),
                                           mean_q10_est = mean(`q10 (est, true)`),
                                           mean_lambda0_err = mean(`lambda0 (est, error)`),
                                           mean_lambda1_err = mean(`lambda1 (est, error)`),
                                           mean_mu0_err = mean(`mu0 (est, error)`),
                                           mean_mu1_err = mean(`mu1 (est, error)`),
                                           mean_q01_err = mean(`q01 (est, error)`),
                                           mean_q10_err = mean(`q10 (est, error)`),
                                           mean_num_tips = mean(`Num Tips`))

# Get average bias of parameter estimates before, after error
bias_estimates <- df_grouped %>% summarise(bias_lambda0_est = mean(`lambda0 (est, true)` - `lambda0 (true)`),
                                           bias_lambda1_est = mean(`lambda1 (est, true)` - `lambda1 (true)`),
                                           bias_mu0_est = mean(`mu0 (est, true)` - `mu0 (true)`),
                                           bias_mu1_est = mean(`mu1 (est, true)` - `mu1 (true)`),
                                           bias_q01_est = mean(`q01 (est, true)` - `q01 (true)`),
                                           bias_q10_est = mean(`q10 (est, true)` - `q10 (true)`),
                                           bias_lambda0_err = mean(`lambda0 (est, error)` - `lambda0 (true)`),
                                           bias_lambda1_err = mean(`lambda1 (est, error)` - `lambda1 (true)`),
                                           bias_mu0_err = mean(`mu0 (est, error)` - `mu0 (true)`),
                                           bias_mu1_err = mean(`mu1 (est, error)` - `mu1 (true)`),
                                           bias_q01_err = mean(`q01 (est, error)` - `q01 (true)`),
                                           bias_q10_err = mean(`q10 (est, error)` - `q10 (true)`))

# Observe correlation between noise and metrics (variance, biases, estimates)

# Save result dataframes to CSV files
all_results_file <- paste(error_rate_str, "all_results.csv", sep = "_")
write.csv(all_results_df, all_results_file)

mean_estimates_file <- paste(error_rate_str, "mean_estimates.csv", sep = "_")
write.csv(mean_estimates, mean_estimates_file)

bias_estimates_file <- paste(error_rate_str, "bias_estimates.csv", sep = "_")
write.csv(bias_estimates, bias_estimates_file)


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
