# Load necessary packages
library(stats4)          # For Maximum Likelihood Estimation (MLE)
library(MASS)            # For generalized statistical functions
library(fitdistrplus)    # For fitting distributions to data
library(statmod)         # For Gauss-Hermite quadrature

# ================================================================
# Function to generate synthetic data
# ---------------------------------------------------------------
# This function creates data where:
# - Residual variability can come from different distributions.
# - Includes optional outliers for testing robustness of methods.
# ================================================================
generate_data <- function(n) {
  # Generate random noise for the residuals (vi)
  vi <- rnorm(n, mean = 0, sd = 5)  # Gaussian noise
  
  # Uncomment the desired distribution for residual variability (ui):
  
  # (1) Half-normal distribution
  # ui <- abs(rnorm(n, mean = 0, sd = 15))
  
  # (2) Exponential distribution
  # ui <- rexp(n, rate = 1/15)
  
  # (3) Gamma distribution
  shape <- 2  # Shape parameter for Gamma distribution
  scale <- 2  # Scale parameter for Gamma distribution
  ui <- rgamma(n, shape = shape, scale = scale)
  
  # (4) Lindley distribution
  rlindley <- function(n, theta) {
    u <- runif(n)
    -1 / theta * log(1 - u * (1 + theta) / 
                       (1 + theta + sqrt((1 + theta)^2 - 4 * u * theta)))
  }
  theta <- 1  # Parameter for Lindley distribution
  # ui <- rlindley(n, theta)
  
  # Introduce 10% outliers by scaling randomly selected residuals
  vi[sample(1:n, size = n * 0.1)] <- vi[sample(1:n, size = n * 0.1)] * 3
  
  # Generate predictors (x) and response variable (y)
  x <- rnorm(n, mean = 5, sd = 2)
  y <- 15 + 20 * x + (vi - ui)  # Linear model with noise
  
  # Return the simulated dataset
  return(data.frame(x = x, y = y))
}

# ================================================================
# Integration methods for handling residuals
# ---------------------------------------------------------------
# Each method estimates model predictions accounting for residual variability.
# - G-HQ: Gauss-Hermite Quadrature
# - MCI: Monte Carlo Integration
# - SMLE: Simulation-based Maximum Likelihood Estimation
# ================================================================

# Gauss-Hermite Quadrature (G-HQ)
ghq_integration <- function(model, points = 10) {
  gh_weights <- gauss.quad(points, kind = "hermite")$weights
  gh_nodes <- gauss.quad(points, kind = "hermite")$nodes
  pred_ghq <- sapply(gh_nodes, function(node) {
    model$coefficients[1] + model$coefficients[2] * model$model$x + 
      node * sd(model$residuals)
  })
  weighted_preds <- pred_ghq %*% gh_weights
  final_pred_ghq <- weighted_preds / sum(gh_weights)
  return(final_pred_ghq)
}

# Monte Carlo Integration (MCI)
mci_integration <- function(model, n_samples = 1000) {
  random_samples <- rnorm(n_samples, mean = 0, sd = sd(model$residuals))
  pred_mci <- model$coefficients[1] + model$coefficients[2] * model$model$x + 
    mean(random_samples)
  return(pred_mci)
}

# Simulation-based MLE (SMLE)
smle_integration <- function(model, n_sim = 1000) {
  simulated_ui <- abs(rnorm(n_sim, mean = 0, sd = sd(model$residuals)))
  pred_smle <- model$coefficients[1] + model$coefficients[2] * model$model$x + 
    mean(simulated_ui)
  return(pred_smle)
}

# ================================================================
# Metrics calculation
# ---------------------------------------------------------------
# Functions to compute RMSE and MAE to evaluate prediction accuracy.
# ================================================================
calculate_metrics <- function(y, y_hat) {
  rmse <- sqrt(mean((y - y_hat)^2))
  mae <- mean(abs(y - y_hat))
  return(c(RMSE = rmse, MAE = mae))
}

# ================================================================
# Simulation function
# ---------------------------------------------------------------
# Simulates model fitting and integration methods over multiple iterations.
# ================================================================
run_simulation_for_methods <- function(n, iterations, model_type) {
  results <- data.frame(Numerical_Method = character(),
                        Sample_Size = integer(),
                        RMSE = numeric(),
                        MAE = numeric(),
                        stringsAsFactors = FALSE)
  
  for (i in 1:iterations) {
    data <- generate_data(n)
    model <- lm(y ~ x, data = data)
    
    # Apply each numerical method and compute metrics
    # G-HQ
    y_hat_ghq <- ghq_integration(model)
    metrics_ghq <- calculate_metrics(data$y, y_hat_ghq)
    results <- rbind(results, data.frame(Numerical_Method = "G-HQ", Sample_Size = n, 
                                         RMSE = metrics_ghq["RMSE"], MAE = metrics_ghq["MAE"]))
    
    # MCI
    y_hat_mci <- mci_integration(model)
    metrics_mci <- calculate_metrics(data$y, y_hat_mci)
    results <- rbind(results, data.frame(Numerical_Method = "MCI", Sample_Size = n, 
                                         RMSE = metrics_mci["RMSE"], MAE = metrics_mci["MAE"]))
    
    # SMLE
    y_hat_smle <- smle_integration(model)
    metrics_smle <- calculate_metrics(data$y, y_hat_smle)
    results <- rbind(results, data.frame(Numerical_Method = "SMLE", Sample_Size = n, 
                                         RMSE = metrics_smle["RMSE"], MAE = metrics_smle["MAE"]))
  }
  
  # Average the results across iterations
  avg_results <- aggregate(. ~ Numerical_Method + Sample_Size, data = results, FUN = mean)
  avg_results$Model <- model_type
  
  return(avg_results)
}

# ================================================================
# Main simulation loop
# ---------------------------------------------------------------
# Executes simulations across various models and sample sizes.
# ================================================================
iterations <- 1000
sample_sizes <- c(10, 250, 1000)
models <- c("C-SFA", "E-SFA", "G-SFA", "L-SFA")

# Store results
final_results <- data.frame(Model = character(),
                            Numerical_Method = character(),
                            Sample_Size = integer(),
                            RMSE = numeric(),
                            MAE = numeric(),
                            stringsAsFactors = FALSE)

# Run simulations
for (model in models) {
  for (sample_size in sample_sizes) {
    print(paste("Running simulation for", model, "with sample size", sample_size))
    model_results <- run_simulation_for_methods(sample_size, iterations, model)
    final_results <- rbind(final_results, model_results)
  }
}

# Display results
print(final_results)

# ================================================================
# ================================================================
# Load necessary packages
library(stats4)        # For maximum likelihood estimation (MLE)
library(MASS)          # For statistical functions, e.g., for fitting distributions
library(fitdistrplus)  # For fitting distributions
library(statmod)       # For Gauss-Hermite quadrature (G-HQ)

# Step 1: Define the bank data (during period)
bank_data_during <- data.frame(
  Bank = c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12"),
  # Variables like E1_during, E2_during, E3_during, etc.
)

# Step 2: Define y (dependent variable) and X (independent variables matrix)
y <- bank_data$E4_pre  # The dependent variable (E4_pre) from the dataset
X <- as.matrix(bank_data[, c("E1_pre", "E2_pre", "E3_pre", "E5_pre", "E6_pre")])  # Independent variables

# Add intercept to matrix X
X <- cbind(1, X)  # Add a column of 1s for the intercept term in the regression

# Step 3: Define a function to calculate error metrics
calculate_metrics <- function(y, y_hat) {
  # Calculate RMSE, MAE, MAPE, and Bias for model performance
  rmse <- sqrt(mean((y - y_hat)^2))
  mae <- mean(abs(y - y_hat))
  mape <- mean(abs((y - y_hat) / y)) * 100
  bias <- mean(y_hat - y)
  return(c(RMSE = rmse, MAE = mae, MAPE = mape, Bias = bias))
}

# Step 4: Define Gauss-Hermite Quadrature (G-HQ) for numerical integration
ghq_integration <- function(model, X, points = 10) {
  # Use Gauss-Hermite quadrature to integrate and calculate predictions
  gh_weights <- gauss.quad(points, kind = "hermite")$weights
  gh_nodes <- gauss.quad(points, kind = "hermite")$nodes
  
  preds <- matrix(0, nrow = nrow(X), ncol = length(gh_nodes))
  
  for (i in 1:length(gh_nodes)) {
    preds[, i] <- X %*% model$coefficients + gh_nodes[i] * sd(model$residuals)
  }
  
  weighted_preds <- preds %*% gh_weights  # Weighted average of predictions
  final_pred_ghq <- weighted_preds
  
  return(final_pred_ghq)
}

# Step 5: Define Monte Carlo Integration (MCI) method
mci_integration <- function(model, X, n_samples = 1000) {
  # Generate random samples and calculate predictions using Monte Carlo integration
  random_samples <- rnorm(n_samples, mean = 0, sd = sd(model$residuals))
  preds <- X %*% model$coefficients + mean(random_samples)
  return(preds)
}

# Step 6: Define Simulated Maximum Likelihood Estimation (SMLE)
smle_integration <- function(model, X, n_sim = 1000) {
  # Simulate residuals and calculate predictions
  simulated_ui <- abs(rnorm(n_sim, mean = 0, sd = sd(model$residuals)))
  preds <- X %*% model$coefficients + mean(simulated_ui)
  return(preds)
}

# Step 7: Define a function to run different SFA models
run_sfa_model <- function(model_type, X, y) {
  # Initialize a dataframe to store the results for each numerical method
  model_results <- data.frame(Numerical_Method = character(), RMSE = numeric(), MAE = numeric(),
                              MAPE = numeric(), Bias = numeric(), stringsAsFactors = FALSE)
  
  # Fit a linear model (adjust as per actual SFA models like C-SFA, E-SFA, etc.)
  model <- lm(y ~ X - 1)  # Fit linear model without intercept term
  
  # Apply numerical methods for Gauss-Hermite Quadrature (G-HQ)
  y_hat_ghq <- ghq_integration(model, X)
  metrics_ghq <- calculate_metrics(y, y_hat_ghq)
  model_results <- rbind(model_results, data.frame(Numerical_Method = "Gauss-Hermite Quadrature", 
                                                   RMSE = metrics_ghq["RMSE"], MAE = metrics_ghq["MAE"], 
                                                   MAPE = metrics_ghq["MAPE"], Bias = metrics_ghq["Bias"]))
  
  # Apply Monte Carlo Integration (MCI)
  y_hat_mci <- mci_integration(model, X)
  metrics_mci <- calculate_metrics(y, y_hat_mci)
  model_results <- rbind(model_results, data.frame(Numerical_Method = "Monte Carlo Integration", 
                                                   RMSE = metrics_mci["RMSE"], MAE = metrics_mci["MAE"], 
                                                   MAPE = metrics_mci["MAPE"], Bias = metrics_mci["Bias"]))
  
  # Apply Simulated Maximum Likelihood Estimation (SMLE)
  y_hat_smle <- smle_integration(model, X)
  metrics_smle <- calculate_metrics(y, y_hat_smle)
  model_results <- rbind(model_results, data.frame(Numerical_Method = "SMLE", 
                                                   RMSE = metrics_smle["RMSE"], MAE = metrics_smle["MAE"], 
                                                   MAPE = metrics_smle["MAPE"], Bias = metrics_smle["Bias"]))
  
  model_results$Model <- model_type  # Add the model type to the results
  return(model_results)
}

# Step 8: Run SFA models and collect results
models <- c("Classical SFA", "Exponential SFA", "Gamma SFA", "Lindley SFA")
results <- data.frame(Model = character(), Numerical_Method = character(), RMSE = numeric(), 
                      MAE = numeric(), MAPE = numeric(), Bias = numeric(), stringsAsFactors = FALSE)

# Loop through each model and store results
for (model in models) {
  model_results <- run_sfa_model(model, X, y)
  results <- rbind(results, model_results)
}

# Step 9: Display results
print(results)  # Output the results of the models and numerical methods


# ================================================================
# Load necessary libraries
library(stats4)   # For statistical modeling
library(MASS)     # For additional statistical methods
library(fitdistrplus)  # For fitting distributions
library(statmod)  # For Gauss-Hermite quadrature

# Step 1: Prepare the bank data
# Define the bank data with multiple attributes
bank_data_during <- data.frame(
  Bank = c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12"),
  # other data columns (E1_during, E2_during, etc.)...
)

# Set the response variable `y` and predictor matrix `X`
y <- bank_data$E4_pre  # Replace E4_pre with the actual response variable name
X <- as.matrix(bank_data[, c("E1_pre", "E2_pre", "E3_pre", "E5_pre", "E6_pre")])  # Predictor variables
X <- cbind(1, X)  # Add intercept column

# Step 2: Define a function to calculate error metrics
calculate_metrics <- function(y, y_hat) {
  # Calculate RMSE (Root Mean Squared Error), MAE (Mean Absolute Error), MAPE (Mean Absolute Percentage Error), and Bias
  rmse <- sqrt(mean((y - y_hat)^2))
  mae <- mean(abs(y - y_hat))
  mape <- mean(abs((y - y_hat) / y)) * 100
  bias <- mean(y_hat - y)
  return(c(RMSE = rmse, MAE = mae, MAPE = mape, Bias = bias))
}

# Step 3: Define Gauss-Hermite Quadrature (G-HQ) integration
ghq_integration <- function(model, X, points = 10) {
  # Perform Gauss-Hermite Quadrature (G-HQ) integration to compute weighted predictions
  # Involves using quadrature weights and nodes
  gh_weights <- gauss.quad(points, kind = "hermite")$weights
  gh_nodes <- gauss.quad(points, kind = "hermite")$nodes
  preds <- matrix(0, nrow = nrow(X), ncol = length(gh_nodes))
  # Loop through each node and apply the model to calculate predictions
  for (i in 1:length(gh_nodes)) {
    preds[, i] <- X %*% model$coefficients + gh_nodes[i] * sd(model$residuals)
  }
  weighted_preds <- preds %*% gh_weights  # Apply the weights
  return(weighted_preds)
}

# Step 4: Define Monte Carlo Integration (MCI)
mci_integration <- function(model, X, n_samples = 1000) {
  # Monte Carlo Integration (MCI) to generate predictions based on random sampling
  random_samples <- rnorm(n_samples, mean = 0, sd = sd(model$residuals))
  preds <- X %*% model$coefficients + mean(random_samples)
  return(preds)
}

# Step 5: Define Simulated Maximum Likelihood Estimation (SMLE)
smle_integration <- function(model, X, n_sim = 1000) {
  # Simulate errors (abs normal distribution) and apply to generate predictions
  simulated_ui <- abs(rnorm(n_sim, mean = 0, sd = sd(model$residuals)))
  preds <- X %*% model$coefficients + mean(simulated_ui)
  return(preds)
}

# Step 6: Define the main SFA model function
run_sfa_model <- function(model_type, X, y) {
  # Initialize a data frame to store the results
  model_results <- data.frame(Numerical_Method = character(), RMSE = numeric(), MAE = numeric(),
                              MAPE = numeric(), Bias = numeric(), stringsAsFactors = FALSE)
  
  # Step 6.1: Fit the model (using linear regression as a placeholder for SFA)
  model <- lm(y ~ X - 1)  # Model fitting (replace with actual SFA fitting method)
  
  # Step 6.2: Apply numerical methods and calculate metrics for each
  # Gauss-Hermite Quadrature
  y_hat_ghq <- ghq_integration(model, X)
  metrics_ghq <- calculate_metrics(y, y_hat_ghq)
  model_results <- rbind(model_results, data.frame(Numerical_Method = "Gauss-Hermite Quadrature", 
                                                   RMSE = metrics_ghq["RMSE"], MAE = metrics_ghq["MAE"], 
                                                   MAPE = metrics_ghq["MAPE"], Bias = metrics_ghq["Bias"]))
  
  # Monte Carlo Integration
  y_hat_mci <- mci_integration(model, X)
  metrics_mci <- calculate_metrics(y, y_hat_mci)
  model_results <- rbind(model_results, data.frame(Numerical_Method = "Monte Carlo Integration", 
                                                   RMSE = metrics_mci["RMSE"], MAE = metrics_mci["MAE"], 
                                                   MAPE = metrics_mci["MAPE"], Bias = metrics_mci["Bias"]))
  
  # Simulated Maximum Likelihood Estimation
  y_hat_smle <- smle_integration(model, X)
  metrics_smle <- calculate_metrics(y, y_hat_smle)
  model_results <- rbind(model_results, data.frame(Numerical_Method = "SMLE", 
                                                   RMSE = metrics_smle["RMSE"], MAE = metrics_smle["MAE"], 
                                                   MAPE = metrics_smle["MAPE"], Bias = metrics_smle["Bias"]))
  
  # Step 6.3: Add the model type to the results and return
  model_results$Model <- model_type
  return(model_results)
}

# Step 7: Define a list of models to run (e.g., Classical SFA, Exponential SFA, etc.)
models <- c("Classical SFA", "Exponential SFA", "Gamma SFA", "Lindley SFA")

# Initialize an empty data frame for storing the final results
results <- data.frame(Model = character(), Numerical_Method = character(), RMSE = numeric(), 
                      MAE = numeric(), MAPE = numeric(), Bias = numeric(), stringsAsFactors = FALSE)

# Step 8: Run the models and store the results
for (model in models) {
  model_results <- run_sfa_model(model, X, y)
  results <- rbind(results, model_results)
}

# Step 9: (Optional) Output results - This is where you'd print or save the results (not done here)
# print(results)   # Uncomment this if you want to view the results after execution


# ================================================================
# Load necessary libraries
library(frontier)      # For Stochastic Frontier Analysis (SFA)
library(fitdistrplus)  # For fitting distributions

# Step 1: Define Lindley distribution functions
# Define the probability density function (PDF) for the Lindley distribution
dlindley <- function(x, theta) {
  # Returns the density function value based on Lindley distribution
}

# Define the cumulative distribution function (CDF) for the Lindley distribution
plindley <- function(q, theta) {
  # Returns the CDF value based on Lindley distribution
}

# Define the random number generation function for Lindley distribution
rlindley <- function(n, theta) {
  # Generates random numbers following the Lindley distribution
}

# Step 2: Plot residuals for distribution check
plot_residuals <- function(residuals) {
  # Plots a histogram of the residuals and marks zero residuals with a red line
}

# Step 3: Fit Lindley distribution to residuals
fit_lindley <- function(residuals) {
  # Ensure all residuals are positive
  if (any(residuals <= 0)) {
    residuals <- residuals - min(residuals) + 0.01  # Shift residuals to positive
  }
  
  # Plot residuals
  plot_residuals(residuals)
  
  # Try fitting the Lindley distribution using fitdistrplus package
  fit <- tryCatch({
    fitdist(residuals, "lindley", start = list(theta = 0.5), 
            densfun = dlindley, 
            pfun = plindley, 
            rfun = rlindley)
  }, error = function(e) {
    warning("Lindley fitting failed: ", e$message)
    return(NULL)
  })
  
  return(fit)
}

# Step 4: Define the bank data (with attributes like E1_pre, E2_pre, etc.)
bank_data <- data.frame(
  Bank = c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12"),
  # Other variables such as E1_pre, E2_pre, etc.
)

# Step 5: Apply log transformation to the bank data (for Cobb-Douglas production function)
bank_data$log_E1 <- log(bank_data$E1_pre)
bank_data$log_E2 <- log(bank_data$E2_pre)
bank_data$log_E3 <- log(bank_data$E3_pre)
bank_data$log_E4 <- log(bank_data$E4_pre)
bank_data$log_E5 <- log(bank_data$E5_pre)
bank_data$log_E6 <- log(bank_data$E6_pre)

# Step 6: Define the Cobb-Douglas frontier model formula
formula <- log_E4 ~ log_E1 + log_E2 + log_E3 + log_E5 + log_E6

# Step 7: Classical SFA Model (C-SFA)
csfa_model <- sfa(formula, data = bank_data)
csfa_efficiency <- efficiencies(csfa_model)

# Step 8: Exponential SFA Model (E-SFA)
esfa_model <- sfa(formula, data = bank_data, truncNorm = TRUE)
esfa_efficiency <- efficiencies(esfa_model)

# Step 9: Apply Gamma SFA using shifted residuals
residuals_csfa <- as.numeric(residuals(csfa_model))
shift_value <- abs(min(residuals_csfa)) + 1  # Ensure positive residuals
shifted_residuals <- residuals_csfa + shift_value

# Fit Gamma distribution to shifted residuals
if (length(shifted_residuals) > 1 && is.numeric(shifted_residuals)) {
  gamma_fit <- fitdist(shifted_residuals, "gamma")
  
  # Gamma efficiency estimation (based on shifted residuals)
  gamma_efficiency <- 1 / (1 + shifted_residuals)  # Example approximation
} else {
  gamma_efficiency <- rep(NA, length(bank_data$Bank))
  warning("Invalid residuals for Gamma distribution fitting.")
}

# Step 10: Fit Lindley distribution to shifted residuals
lindley_fit <- fit_lindley(shifted_residuals)

# Lindley efficiency estimation
if (!is.null(lindley_fit)) {
  lindley_efficiency <- 1 / (1 + shifted_residuals)  # Example approximation
} else {
  lindley_efficiency <- rep(NA, length(bank_data$Bank))
}

# Step 11: Combine results into a data frame
efficiency_table <- data.frame(
  Bank = bank_data$Bank,
  CSFA = csfa_efficiency,
  ESFA = esfa_efficiency,
  GSFA = gamma_efficiency,  # Gamma efficiency estimation
  LSFA = lindley_efficiency  # Lindley efficiency estimation
)

# Step 12: (Optional) Output results
# print(efficiency_table)   # Uncomment to display results

# Load necessary packages
library(stats4)
library(MASS)
library(fitdistrplus)
library(statmod)  # For Gauss-Hermite quadrature

# Function to generate data
# Function to generate data with more variability in residuals
generate_data <- function(n) {
  vi <- rnorm(n, mean = 0, sd = 5)
  
  # Gamma distribution for ui
  shape <- 2  # Adjust the shape parameter as needed
  scale <- 2  # Adjust the scale parameter as needed
  ui <- rgamma(n, shape = shape, scale = scale)
  
  # Introduce more variability in residuals by adding some larger outliers
  vi[sample(1:n, size = n * 0.1)] <- vi[sample(1:n, size = n * 0.1)] * 3  # if necessary 10% outliers scaled up
  
  # Generate x and y values
  x <- rnorm(n, mean = 5, sd = 2)
  y <- 15 + 20 * x + (vi - ui)
  return(data.frame(x = x, y = y))
}

# Define the integration methods

# Gauss-Hermite Quadrature (G-HQ) integration
ghq_integration <- function(model, points = 10) {
  gh_weights <- gauss.quad(points, kind = "hermite")$weights
  gh_nodes <- gauss.quad(points, kind = "hermite")$nodes
  pred_ghq <- sapply(gh_nodes, function(node) {
    model$coefficients[1] + model$coefficients[2] * model$model$x + node * sd(model$residuals)
  })
  weighted_preds <- pred_ghq %*% gh_weights
  final_pred_ghq <- weighted_preds / sum(gh_weights)
  return(final_pred_ghq)
}

# Monte Carlo Integration (MCI)
mci_integration <- function(model, n_samples = 1000) {
  random_samples <- rnorm(n_samples, mean = 0, sd = sd(model$residuals))
  pred_mci <- model$coefficients[1] + model$coefficients[2] * model$model$x + mean(random_samples)
  return(pred_mci)
}

# Simulated Maximum Likelihood Estimation (SMLE)
smle_integration <- function(model, n_sim = 1000) {
  simulated_ui <- abs(rnorm(n_sim, mean = 0, sd = sd(model$residuals)))
  pred_smle <- model$coefficients[1] + model$coefficients[2] * model$model$x + mean(simulated_ui)
  return(pred_smle)
}

# Define function to calculate metrics (RMSE, MAE)
calculate_metrics <- function(y, y_hat) {
  rmse <- sqrt(mean((y - y_hat)^2))
  mae <- mean(abs(y - y_hat))
  return(c(RMSE = rmse, MAE = mae))
}

# Simulation function for methods and models
run_simulation_for_methods <- function(n, iterations, model_type) {
  results <- data.frame(Numerical_Method = character(),
                        Sample_Size = integer(),
                        RMSE = numeric(),
                        MAE = numeric(),
                        stringsAsFactors = FALSE)
  
  for (i in 1:iterations) {
    data <- generate_data(n)
    model <- lm(y ~ x, data = data)
    
    # G-HQ
    y_hat_ghq <- ghq_integration(model)
    metrics_ghq <- calculate_metrics(data$y, y_hat_ghq)
    results <- rbind(results, data.frame(Numerical_Method = "G-HQ", Sample_Size = n, RMSE = metrics_ghq["RMSE"], MAE = metrics_ghq["MAE"]))
    
    # MCI
    y_hat_mci <- mci_integration(model)
    metrics_mci <- calculate_metrics(data$y, y_hat_mci)
    results <- rbind(results, data.frame(Numerical_Method = "MCI", Sample_Size = n, RMSE = metrics_mci["RMSE"], MAE = metrics_mci["MAE"]))
    
    # SMLE
    y_hat_smle <- smle_integration(model)
    metrics_smle <- calculate_metrics(data$y, y_hat_smle)
    results <- rbind(results, data.frame(Numerical_Method = "SMLE", Sample_Size = n, RMSE = metrics_smle["RMSE"], MAE = metrics_smle["MAE"]))
  }
  
  # Average the results across iterations
  avg_results <- aggregate(. ~ Numerical_Method + Sample_Size, data = results, FUN = mean)
  avg_results$Model <- model_type
  
  return(avg_results)
}

# Set simulation parameters
iterations <- 1000
sample_sizes <- c(10, 250, 1000)
models <- c("C-SFA", "E-SFA", "G-SFA", "L-SFA")

# Initialize results storage
final_results <- data.frame(Model = character(),
                            Numerical_Method = character(),
                            Sample_Size = integer(),
                            RMSE = numeric(),
                            MAE = numeric(),
                            stringsAsFactors = FALSE)

# Run simulation for each model and sample size
for (model in models) {
  for (sample_size in sample_sizes) {
    print(paste("Running simulation for", model, "with sample size", sample_size))
    model_results <- run_simulation_for_methods(sample_size, iterations, model)
    final_results <- rbind(final_results, model_results)
  }
}

# Format and display final results
print(final_results)


# ================================================================
# Step 1: Load Necessary Libraries
# - `ggplot2` for plotting
# - `reshape2` for reshaping data

# Load required libraries
library(ggplot2)
library(reshape2)

# Step 2: Prepare Data for Visualization
# Create a data frame with bank names and efficiency scores
# Define the columns for Pre-COVID and Dur-COVID scores for each method (CSFA, ESFA, GSFA, LSFA)
efficiency_scores <- data.frame(
  Bank = c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12"),
  Pre_CSFA = c( ... ),  # Pre-COVID CSFA efficiency scores
  Pre_ESFA = c( ... ),  # Pre-COVID ESFA efficiency scores
  Pre_GSFA = c( ... ),  # Pre-COVID GSFA efficiency scores
  Pre_LSFA = c( ... ),  # Pre-COVID LSFA efficiency scores
  Dur_CSFA = c( ... ),  # Dur-COVID CSFA efficiency scores
  Dur_ESFA = c( ... ),  # Dur-COVID ESFA efficiency scores
  Dur_GSFA = c( ... ),  # Dur-COVID GSFA efficiency scores
  Dur_LSFA = c( ... )   # Dur-COVID LSFA efficiency scores
)

# Step 3: Calculate Percentage Changes
# Calculate the percentage change for each efficiency method (CSFA, ESFA, GSFA, LSFA)
efficiency_scores$Change_CSFA <- (efficiency_scores$Dur_CSFA - efficiency_scores$Pre_CSFA) / efficiency_scores$Pre_CSFA * 100
efficiency_scores$Change_ESFA <- (efficiency_scores$Dur_ESFA - efficiency_scores$Pre_ESFA) / efficiency_scores$Pre_ESFA * 100
efficiency_scores$Change_GSFA <- (efficiency_scores$Dur_GSFA - efficiency_scores$Pre_GSFA) / efficiency_scores$Pre_GSFA * 100
efficiency_scores$Change_LSFA <- (efficiency_scores$Dur_LSFA - efficiency_scores$Pre_LSFA) / efficiency_scores$Pre_LSFA * 100

# Step 4: Reshape Data for Plotting
# Use `melt()` to reshape the data from wide format to long format
# This makes it suitable for plotting with ggplot2
visual_data <- melt(efficiency_scores[, c("Bank", "Change_CSFA", "Change_ESFA", "Change_GSFA", "Change_LSFA")], id.vars = "Bank")

# Step 5: Create the Bar Plot
# Use `ggplot()` to create a bar plot to show the percentage changes in efficiency scores
# Different methods (CSFA, ESFA, GSFA, LSFA) will be represented by different colors
ggplot(visual_data, aes(x = Bank, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  # Bar plot with dodged bars
  theme_minimal() +  # Clean background theme
  labs(
    title = "Percentage Changes in Efficiency Scores Pre-COVID19 vs Dur-COVID19",
    x = "Banks",   # X-axis label
    y = "Percentage Change (%)",  # Y-axis label
    fill = "Efficiency Method"  # Legend title
  ) +
  scale_fill_manual(  # Custom color for each method
    values = c("Change_CSFA" = "#F8766D", "Change_ESFA" = "#00BFC4", 
               "Change_GSFA" = "#7CAE00", "Change_LSFA" = "#C77CFF"),
    labels = c("CSFA", "ESFA", "GSFA", "LSFA")
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
