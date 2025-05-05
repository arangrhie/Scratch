# Load necessary libraries
library(ggplot2)
library(reshape2)

# Read the data
data <- read.table("input/inova.copynum_dip.4187.txt", header = TRUE)

# Melt the data
data_melt <- melt(data, id.vars = "Sample")

# Plot distributions
ggplot(data_melt, aes(x = variable, y = value)) +
	geom_violin() +
	geom_boxplot(width = 0.1, outlier.shape = NA) +
	geom_whisker() +
	geom_jitter(alpha = 0.3, size = 0.2)