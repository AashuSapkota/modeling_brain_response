# Install necessary packages if not already installed
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("GGally")) install.packages("GGally")
if (!require("viridis")) install.packages("viridis")
if (!require("ggExtra")) install.packages("ggExtra")



library(ggplot2)
library(dplyr)
library(tidyr)
library(GGally)
library(viridis)
library(ggExtra)



data <- read.csv("E:/MSC/Stat/data.csv")

# Select the relevant columns for signals x1, x2, x3, x4, x5
data_selected <- data %>% select(x1, x2, x3, x4, x5)

# Display the structure of the dataset
str(data_selected)

# Display the first few rows of the dataset
head(data_selected)

# Summary statistics for each column
summary(data_selected)

# Check for NA values in the dataset
colSums(is.na(data_selected))

# Assuming there is a time or trial column in your dataset, if not, we can assume a sequence for time
# If the dataset has a 'time' or 'trial' column, replace 'time' below with the correct column name
data_selected$time <- 1:nrow(data_selected) # Simulating time from 1 to number of rows

# Reshape the data from wide format to long format for ggplot
data_long <- data_selected %>%
  gather(key = "Signal", value = "BloodOxygenation", x1, x2, x3, x4, x5)


# Task1.1 Time series plots (of input and output signal) 
# Create the time series plot
ggplot(data_long, aes(x = time, y = BloodOxygenation, color = Signal)) +
  geom_line() +
  labs(title = "Time Series of Brain Activity (Blood Oxygenation)",
       x = "Time (in seconds or trials)",
       y = "Blood Oxygenation Level") +
  theme_minimal() +
  theme(legend.title = element_blank())

# Create the time series plot with faceting
ggplot(data_long, aes(x = time, y = BloodOxygenation, color = Signal)) +
  geom_line() +
  labs(title = "Time Series of Brain Activity (Blood Oxygenation)",
       x = "Time (in seconds or trials)",
       y = "Blood Oxygenation Level") +
  theme_minimal() +
  theme(legend.position = "none") +  # Remove legend since each facet is labeled
  facet_wrap(~ Signal, ncol = 1)  # Create separate facets for each Signal

# Task1.2 Distribution for each signal (time-series) 
# Density plot with faceting for each signal
ggplot(data_long, aes(x = BloodOxygenation, color = Signal)) +
  geom_density(alpha = 0.3) +  # Semi-transparent density curves
  labs(title = "Distribution of Blood Oxygenation Levels for Each Signal",
       x = "Blood Oxygenation Level",
       y = "Density") +
  theme_minimal()
  theme(legend.title = element_blank())
  
  
# Individual Bins  
ggplot(data_long, aes(x = BloodOxygenation, fill = Signal)) +
    geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.5, color = "black") +
    geom_density(alpha = 0.7) +
    facet_wrap(~ Signal, scales = "free", ncol=2) +  # Separate panel for each signal
    labs(title = "Individual Distribution of Blood Oxygenation Levels",
         x = "Blood Oxygenation Level",
         y = "Density") +
    theme_minimal() +
    theme(legend.position = "none")
  
# Task 1.3 Correlation and scatter plots (between different combination of input and output signals) to
# examine their dependencies

# Select relevant columns for correlation
input_output <- data_selected %>% select(x2, x1, x3, x4, x5)

# Compute correlation matrix
cor_matrix <- cor(input_output, use = "complete.obs")

# Print the correlation matrix
print(cor_matrix)

# Reshape data for faceted plotting
data_long <- data_selected %>%
  select(x2, x1, x3, x4, x5) %>%
  pivot_longer(cols = c(x1, x3, x4, x5), names_to = "Input", values_to = "Value")

# Create faceted scatter plots
ggplot(data_long, aes(x = Value, y = x2)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(~Input, scales = "free_x") +
  labs(title = "Scatter Plots of Inputs vs Output", x = "Inputs", y = "Output (x2)") +
  theme_minimal()


# Pairwise scatter plots and correlations
ggpairs(input_output,
        lower = list(continuous = "smooth"),
        diag = list(continuous = "densityDiag"),
        upper = list(continuous = "cor"),
        title = "Pairwise Relationships between Inputs and Output")

# Task 2
# Selecting input and output signals
y <- data_selected$x2 # output variable x2
x1 <- data_selected$x1 # input variable x1
x3 <- data_selected$x3 # input variable x2
x4 <- data_selected$x4 # input variable x3
x5 <- data_selected$x5 # input variable x4

# Define the bias term
bias <- rep(1, length(y))

# Add Gaussian noise to the response variable
set.seed(42)  # For reproducibility
noise_sd <- 0.3  # Standard deviation of the noise
noise <- rnorm(length(y), mean = 0, sd = noise_sd)  # Generate Gaussian noise

# Noisy y variable (the true value + noise)
y_noisy <- y + noise

# List of models and corresponding design matrices
models <- list(
  list(X = cbind(x4, x3^2, bias), name = "Model 1"),
  list(X = cbind(x4, x3^2, x5, bias), name = "Model 2"),
  list(X = cbind(x3, x4, x5^3), name = "Model 3"),
  list(X = cbind(x4, x3^2, x5^3, bias), name = "Model 4"),
  list(X = cbind(x4, x1^2, x3^2, bias), name = "Model 5")
)

# Task 2.1
# Function to compute theta and RSS for a given design matrix

compute_model_with_lse <- function(X, y) {
  theta <- solve(t(X) %*% X) %*% t(X) %*% y  # Compute theta
  return(list(theta = theta))
}


# Compute theta and LSE for each model
results_lse <- lapply(models, function(model) {
  compute_model_with_lse(model$X, y_noisy)
})

# Print results LSE
for (i in 1:length(results_lse)) {
  cat(models[[i]]$name, "\n")  # Print model name
  cat("Thetahat:\n", results_lse[[i]]$theta, "\n")
}


# Task 2.2
# Function to compute theta and RSS for a given design matrix
compute_model_with_rss <- function(X, y) {
  theta <- solve(t(X) %*% X) %*% t(X) %*% y  # Compute theta
  y_hat <- X %*% theta  # Predicted values
  residuals <- y - y_hat  # Residuals
  rss <- sum(residuals^2)  # Compute residual sum of squares (RSS)
  return(list(theta = theta, rss = rss, y_hat = y_hat, residuals = residuals))
}

# Compute theta and RSS for each model
results_rss <- lapply(models, function(model) {
  compute_model_with_rss(model$X, y_noisy)
})


# Print results RSS
for (i in 1:length(results_rss)) {
  cat(models[[i]]$name, "\n")  # Print model name
  cat("Theta:\n", results_rss[[i]]$theta, "\n")
  cat("RSS:", results_rss[[i]]$rss, "\n\n")
}


# Task 2.3: Function to compute Log-Likelihood
compute_log_likelihood <- function(X, y) {
  # Compute theta (least squares estimate)
  theta <- solve(t(X) %*% X) %*% t(X) %*% y

  # Compute residuals and RSS
  y_hat <- X %*% theta  # Predicted values
  residuals <- y - y_hat  # Residuals
  rss <- sum(residuals^2)  # Residual sum of squares (RSS)

  # Estimate the variance (sigma^2) from the residuals using RSS / (n - 1)
  n <- length(y)
  sigma_squared <- rss / (n - 1)

  # Compute the log-likelihood
  log_likelihood <- -n / 2 * log(2 * pi * sigma_squared) - rss / (2 * sigma_squared)

  return(list(theta = theta, rss = rss, log_likelihood = log_likelihood, sigma_squared = sigma_squared))
}

# Compute Log-Likelihood for each model (Task 2.3)
results_log_likelihood <- lapply(models, function(model) {
  compute_log_likelihood(model$X, y_noisy)
})

# Print results (Log-Likelihood for each model)
cat("Task 2.3: Log-Likelihood for each model\n")
for (i in 1:length(results_log_likelihood)) {
  cat(models[[i]]$name, "\n")
  cat("Theta Hat:\n", results_log_likelihood[[i]]$theta, "\n")
  cat("Log-Likelihood:", results_log_likelihood[[i]]$log_likelihood, "\n")
  cat("Estimated Sigma^2:", results_log_likelihood[[i]]$sigma_squared, "\n\n")
}

# Task 2.4:

# Function to compute AIC and BIC
compute_aic_bic <- function(log_likelihood, k, n) {
  # AIC = 2k - 2 * log-likelihood
  aic <- 2 * k - 2 * log_likelihood

  # BIC = k * log(n) - 2 * log-likelihood
  bic <- k * log(n) - 2 * log_likelihood

  return(list(aic = aic, bic = bic))
}

# Number of observations (n)
n <- length(y)

# Define the number of parameters (k) for each model
k_values <- c(3, 4, 3, 4, 4)  # Number of parameters for each model

# Compute AIC and BIC for each model (Task 2.4)
results_aic_bic <- mapply(function(log_likelihood, k) {
  compute_aic_bic(log_likelihood, k, n)
}, log_likelihood = sapply(results_log_likelihood, function(x) x$log_likelihood),
k = k_values, SIMPLIFY = FALSE)

# Print results for Task 2.4 (AIC and BIC)
cat("Task 2.4: AIC and BIC for each model\n")
for (i in 1:length(results_aic_bic)) {
  cat(models[[i]]$name, "\n")
  cat("AIC:", results_aic_bic[[i]]$aic, "\n")
  cat("BIC:", results_aic_bic[[i]]$bic, "\n\n")
}



# # Task 2.5: Plot residuals distribution for each model
# Function to plot residuals distribution (Histogram and Q-Q plot)
plot_residuals_distribution <- function(residuals, model_name) {
  # Plot Histogram of Residuals
  hist(residuals, breaks = 20, main = paste("Residuals Distribution for", model_name),
       xlab = "Residuals", col = "lightblue", border = "black", probability = TRUE)
  lines(density(residuals), col = "red", lwd = 2)  # Add density curve

  # Q-Q plot
  qqnorm(residuals, main = paste("Q-Q Plot for", model_name))
  qqline(residuals, col = "red", lwd = 2)  # Add reference line
}

# Task 2.5: Plot residuals distribution for each model
cat("Task 2.5: Checking Distribution of Model Prediction Errors (Residuals)\n")

# Set up a 5-row by 2-column layout for the residuals distribution (Histogram and Q-Q Plot)
par(mfrow = c(5, 2), mar = c(4, 4, 2, 1))  # Adjust margins for better spacing

# Loop through each model to plot residuals distribution
for (i in 1:length(results_rss)) {
  # Get residuals from results_rss for each model
  residuals <- results_rss[[i]]$residuals
  model_name <- models[[i]]$name

  # Plot the residuals distribution (Histogram and Q-Q plot)
  plot_residuals_distribution(residuals, model_name)

  # Print some basic residual statistics
  cat(model_name, "Residuals Summary:\n")
  cat("Mean of residuals: ", mean(residuals), "\n")
  cat("Standard deviation of residuals: ", sd(residuals), "\n\n")
}

# Reset the plotting layout
par(mfrow = c(1, 1))  # Reset to default layout


# Task 2.7: Train-Test Split and Model Validation for Best Model (Model 2)

# Create indices for training and testing data sets (70% train, 30% test)
train_indices <- sample(1:nrow(data_selected), size = 0.7 * nrow(data_selected))
train_data <- data_selected[train_indices, ]
test_data <- data_selected[-train_indices, ]

# Extract the output variable (y) and inputs (x3, x4, x5, bias) for Model 2
train_y <- train_data$x2
test_y <- test_data$x2

train_X <- cbind(train_data$x4, train_data$x3^2, train_data$x5, rep(1, nrow(train_data))) # Design matrix
test_X <- cbind(test_data$x4, test_data$x3^2, test_data$x5, rep(1, nrow(test_data)))

# Function to compute model parameters (theta) using training data
compute_theta <- function(X, y) {
  solve(t(X) %*% X) %*% t(X) %*% y
}

# Train the model on training data
theta_hat <- compute_theta(train_X, train_y)

# Predict output on test data
test_predictions <- test_X %*% theta_hat

# Compute residuals on the training data
train_y_hat <- train_X %*% theta_hat
train_residuals <- train_y - train_y_hat
sigma_squared <- sum(train_residuals^2) / (nrow(train_data) - ncol(train_X))

# Compute the standard errors of predictions
pred_var <- rowSums((test_X %*% solve(t(train_X) %*% train_X)) * test_X)
pred_std_err <- sqrt(sigma_squared * (1 + pred_var))


# Compute 95% confidence intervals for predictions
alpha <- 0.05
z_value <- qnorm(1 - alpha / 2)
conf_lower <- test_predictions - z_value * pred_std_err
conf_upper <- test_predictions + z_value * pred_std_err


# Combine test data, predictions, and confidence intervals for plotting
test_results <- data.frame(
  Actual = test_y,
  Predicted = test_predictions,
  Lower_CI = conf_lower,
  Upper_CI = conf_upper
)

# Plot predictions, confidence intervals, and actual test data
library(ggplot2)

ggplot(test_results, aes(x = 1:nrow(test_results))) +
  geom_point(aes(y = Actual), color = "blue", size = 2, alpha = 0.7, shape = 1) +
  geom_line(aes(y = Predicted), color = "red", size = 1) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.2, fill = "orange") +
  labs(title = "Model 2: Predictions vs Actual Test Data",
       x = "Test Data Points",
       y = "Output (x2)") +
  theme_minimal()

# Task3
# Task 3.1
# Least squares estimates from Task 2.1 (as per your code)
theta_hat_model2 <- c(0.7471398, 0.03920352, 0.1803357, -0.0641826)
rss_model2 <- 151.4747

# Calculate the absolute values of the parameters
abs_theta <- abs(theta_hat_model2)

# Identify the indices of the two parameters with the largest absolute values
top_two_indices <- order(abs_theta, decreasing = TRUE)[1:2]

# Extract the parameter names for the largest two absolute values
top_two_params <- theta_hat_model2[top_two_indices]

# Task 3.2
# Prior distributions for theta1 and theta3 (assuming uniform priors)
prior_theta1 <- function() runif(1, min = -abs(top_two_params[1]), max = abs(top_two_params[1]))
prior_theta3 <- function() runif(1, min = -abs(top_two_params[2]), max = abs(top_two_params[2]))

# Task 3.3
# Simulate the model with fixed parameters (fix theta2 and theta4, sample theta1 and theta3)
simulate_model <- function(theta1, theta3, X) {
  # Fixed parameters from Task 2.1
  theta2 <- 0.03920352  # Fixed value of theta2 from Task 2.1
  theta4 <- -0.0641826  # Fixed value of theta4 from Task 2.1
  
  # Recompute the fitted values based on the new parameter set
  # y_sim <- theta1 * X[,1] + theta2 * X[,2]^2 + theta3 * X[,3] + theta4 * X[,4]
  y_sim <- theta1 * X[,1] + theta2 * X[,2] + theta3 * X[,3] + theta4 * X[,4]
  return(y_sim)
}

# ABC Rejection Sampling
set.seed(42)  # For reproducibility
n_simulations <- 10000  # Number of simulations

# tolerance <- 600  # Distance tolerance for accepting a sample
tolerance <- rss_model2 * 4
print(tolerance)
accepted_samples <- list()  # Store accepted parameter samples

# Loop through rejection sampling
for (i in 1:n_simulations) {
  # Sample theta1 and theta3 from their priors
  theta1_sample <- prior_theta1()
  theta3_sample <- prior_theta3()
  
  # Simulate the model output with the sampled parameters
  y_sim <- simulate_model(theta1_sample, theta3_sample, models[[2]]$X)  # Using Model 2 design matrix
  
  # Compute the distance (sum of squared residuals) between simulated and observed data
  residuals <- y_noisy - y_sim  # Using the noisy data from Task 2
  distance <- sum(residuals^2)
  
  # Debugging output to check distance and residuals
  if (i %% 1000 == 0) {  # Print every 1000th iteration for clarity
    cat("Iteration:", i, "Distance:", distance, "\n")
  }
  
  # If the distance is within tolerance, accept the sample
  if (distance < tolerance) {
    accepted_samples[[length(accepted_samples) + 1]] <- c(theta1_sample, theta3_sample)
  }
}

# Convert the accepted samples into a data frame for easier analysis
accepted_samples_df <- do.call(rbind, accepted_samples)
accepted_samples_df <- as.data.frame(accepted_samples_df)


# Task 3.4
# Plot the accepted samples for theta1 and theta3
library(ggplot2)
p <- ggplot(accepted_samples_df, aes(x = V1, y = V2)) +
  geom_point(alpha = 0.5, color = "blue") +
  labs(title = "Accepted Samples from ABC Rejection Sampling",
       x = "Theta 1", y = "Theta 3") +
  theme_minimal()

# Add marginal histograms
ggExtra::ggMarginal(p, type = "histogram", fill = "blue", alpha = 0.5)

# Summary of accepted samples
cat("Number of accepted samples: ", nrow(accepted_samples_df), "\n")
summary(accepted_samples_df)
