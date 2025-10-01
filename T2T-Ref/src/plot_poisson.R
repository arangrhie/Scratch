#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)

# Generate Poisson distributions
# For Poisson distribution, lambda is both the mean and the parameter
# Lambda = 10 will have peak around 10
# Lambda = 8 will have peak around 8

dat10=100
dat8=80

# Create data for both distributions
x_values <- 0:200  # Range of x values to plot

# Generate Poisson probabilities
poisson_data <- data.frame(
  x = rep(x_values, 2),
  probability = c(
    dpois(x_values, lambda = dat10),  # Poisson with lambda = 100
    dpois(x_values, lambda = dat8)    # Poisson with lambda = 80
  ),
  distribution = rep(c("Copy = 10", "Copy = 8"), each = length(x_values))
)

# Create the plot
p1 <- ggplot(poisson_data, aes(x = x, y = probability, fill = distribution)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  labs(
    title = "Poisson Distributions with Different Lambda Values",
    x = "Value (k)",
    y = "Probability P(X = k)",
    fill = "Distribution"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "top"
  ) +
  scale_fill_manual(values = c("Copy = 10" = "steelblue", "Copy = 8" = "coral")) +
  scale_x_continuous(breaks = seq(0, 200, by = 20))

# Print the plot
print(p1)

# Save the plot
ggsave("../output/poisson_distribution_plot.pdf", plot = p1, width = 10, height = 6, dpi = 300)

# Alternative: Line plot instead of bar plot
p2 <- ggplot(poisson_data, aes(x = x, y = probability, color = distribution)) +
  geom_line(linewidth = 1.2) +
  # geom_point(size = 2) +
  labs(
    title = "Poisson Distributions - Line Plot",
    x = "Value (k)",
    y = "Probability P(X = k)",
    color = "Distribution"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "top"
  ) +
  scale_color_manual(values = c("Copy = 10" = "steelblue", "Copy = 8" = "coral")) +
  scale_x_continuous(breaks = seq(0, 200, by = 20))


print(p2)
ggsave("../output/poisson_distribution_line_plot.pdf", plot = p2, width = 10, height = 6, dpi = 300)

# Generate random samples and create histogram
set.seed(123)
n_samples <- 1000

sample_data <- data.frame(
  values = c(
    rpois(n_samples, lambda = dat10),
    rpois(n_samples, lambda = dat8)
  ),
  distribution = rep(c("Lambda = 10", "Lambda = 8"), each = n_samples)
)

# Histogram of random samples
p3 <- ggplot(sample_data, aes(x = values, fill = distribution)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 20) +
  facet_wrap(~distribution, ncol = 1) +
  labs(
    title = "Random Samples from Poisson Distributions",
    x = "Value",
    y = "Frequency",
    fill = "Distribution"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("Copy = 10" = "steelblue", "Copy = 8" = "coral")) +
  scale_x_continuous(breaks = seq(0, 200, by = 20))

print(p3)
ggsave("../output/poisson_samples_histogram.pdf", plot = p3, width = 10, height = 8, dpi = 300)

# Print summary statistics
cat("\nSummary Statistics:\n")
cat("Lambda = 8 distribution:\n")
cat("  Mean:", mean(dpois(0:200, dat8) * 0:200), "\n")
cat("  Variance:", dat8, "\n\n")

cat("Lambda = 10 distribution:\n")
cat("  Mean:", mean(dpois(0:200, dat10) * 0:200), "\n")
cat("  Variance:", dat10, "\n")
