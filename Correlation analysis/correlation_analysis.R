library(tidyverse)
library(RColorBrewer)
library(pals)
library(ggthemes)
library(ggrepel)
library(gtable)
library(grid)
library(gridExtra)
library(dplyr)

setwd('/path/to/input/')

if (!dir.exists("correlation_plots")) {
  dir.create("correlation_plots")
}

filename = 'file name'
picname = 'plot name'
savepath = "/path/to/output/"

file_1 = 'counts_repaired.txt'
data_1 = read.table(file_1, header=TRUE, sep='\t')
data_1$counts = rowSums(data_1[2:5])

file_2 = 'counts_repaired_2.txt'
data_2 = read.table(file_2, header=TRUE, sep='\t')
data_2$counts = rowSums(data_2[2:5])

data_1_filtered <- subset(data_1, bp>=336 & bp<=8714)
data_2_filtered <- subset(data_2, bp>=336 & bp<=8714)

cross_correlation <- function(a, b, s) {
  n <- length(a)
  mu_a <- mean(a)
  mu_b <- mean(b)
  sigma_a <- sd(a)
  sigma_b <- sd(b)
  
  if (s >= 0) {
    sum((a[1:(n-s)] - mu_a) * (b[(1+s):n] - mu_b)) / (sigma_a * sigma_b * (n - s))
  } else {
    sum((a[(1-s):n] - mu_a) * (b[1:(n+s)] - mu_b)) / (sigma_a * sigma_b * (n + s))
  }
}

separations <- -400:400
n_correlations <- sapply(separations, function(s) cross_correlation(data_1_filtered$counts, data_2_filtered$counts, s))

df <- data.frame(separation = separations, correlation = n_correlations)


#  bootstrap and calculate correlations
bootstrap_correlation <- function(a, b, n_bootstrap) {
  boot_correlations <- matrix(NA, nrow = n_bootstrap, ncol = length(separations))
  n <- length(a)
  
  for (i in 1:n_bootstrap) {
    sample_indices <- sample(1:n, n, replace = TRUE)
    boot_correlations[i, ] <- sapply(separations, function(s) cross_correlation(a[sample_indices], b[sample_indices], s))
  }
  
  boot_correlations
}

set.seed(50)
n_bootstrap = 200
bootstrap_results <- bootstrap_correlation(data_1_filtered$counts, data_2_filtered$counts, n_bootstrap)

probs <- c("low_001" = 0.0005, "high_001" = 0.9995,
           "low_01" = 0.005, "high_01" = 0.995,
           "low_05" = 0.025, "high_05" = 0.975)
ci_limits <- apply(bootstrap_results, 2, quantile, probs = probs)

ci_data <- data.frame(
  separation = rep(separations, times = 3),
  ci_lower = c(ci_limits[1,], ci_limits[3,], ci_limits[5,]),
  ci_upper = c(ci_limits[2,], ci_limits[4,], ci_limits[6,]),
  level = factor(rep(alpha, each = length(separations)))
)

max_correlations <- apply(bootstrap_results, 2, max)
min_correlations <- apply(bootstrap_results, 2, min)

ci_data$max_correlation <- rep(max_correlations, times = 3)
ci_data$min_correlation <- rep(min_correlations, times = 3)


filtered_rows <- ci_data %>% filter(separation == 0 | separation == 1 | separation == -1)
print(filtered_rows)

breaks <- seq(-0.5, 1, by = 0.1)
labels <- ifelse(breaks %in% c(-0.5, 0, 0.5, 1), as.character(breaks), "")

g <- ggplot(df, aes(x = separation, y = correlation)) +
  geom_line() +
  ggtitle("Correlation Analysis") +
  xlab("Separation (nucleotides)") +
  ylab("Correlation Function") +
  theme_classic() +
  scale_y_continuous(breaks=breaks, labels=labels, limits = c(-0.5, 1.0)) +
  geom_hline(yintercept = c(-0.5, 0, 0.5, 1.0), linetype="dashed", color="#bebebe82", alpha = 0.5) +
  theme(
    legend.position = 'none',
    axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(size=13),
    axis.text.y = element_text(size=13)
  )


g <- g + 
  geom_area(data = subset(ci_data, level == 0.001),
            aes(x = separation, y = ci_lower), fill = "grey", alpha = 0.3) +
  geom_area(data = subset(ci_data, level == 0.001),
            aes(x = separation, y = ci_upper), fill = "grey", alpha = 0.3) +
  geom_area(data = subset(ci_data, level == 0.01),
            aes(x = separation, y = ci_lower), fill = "grey", alpha = 0.3) +
  geom_area(data = subset(ci_data, level == 0.01),
            aes(x = separation, y = ci_upper), fill = "grey", alpha = 0.3) +
  geom_area(data = subset(ci_data, level == 0.05),
            aes(x = separation, y = ci_lower), fill = "grey", alpha = 0.3) +
  geom_area(data = subset(ci_data, level == 0.05),
            aes(x = separation, y = ci_upper), fill = "grey", alpha = 0.3)

g <- g + scale_fill_manual(values = c("0.001" = "grey46", "0.01" = "grey46", "0.05" = "grey46"),
                           name = "Significance Level",
                           breaks = levels(ci_data$level),
                           labels = paste("p <", levels(ci_data$level)))

g <- g + annotate("text", x = 300, y = 0.9, label = "p > 0.05", col = "grey46", size = 4)

g <- g + annotate("text", x = -200, y = 0.9, label = picname, col = "black", size = 3)

ggsave(paste(savepath,filename,'.png',sep=''), plot = g, width = 4, height = 4, units='in', dpi=300)