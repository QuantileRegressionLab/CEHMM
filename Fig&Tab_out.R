##### Preliminary Analysis of Real Data #####
# Load required libraries
library(PerformanceAnalytics)
library(tseries)
library(stats)

# Load the returns data
load("data/returns.RData")

# Remove any rows with NA values from the returns data
ret_all <- na.omit(ret.df) 

# Check for any remaining NA values
sum(as.numeric(is.na(ret_all)))

# Calculate and select specific statistics for the returns data
table_stats <- table.Stats(ret_all[,])
table_stats <- table_stats[c(3, 6, 9, 14, 15, 16),]

# Calculate the Jarque-Bera and Augmented Dickey-Fuller test statistics for each return series
m <- dim(ret_all)[2]
jbt <- sapply(1:m, function(i) jarque.bera.test(ret_all[, i])$statistic)
jbt_pv <- sapply(1:m, function(i) jarque.bera.test(ret_all[, i])$p.value)
jbt <- as.numeric(jbt)

adf_stat <- sapply(1:m, function(i) as.numeric(adf.test(ret_all[, i])$statistic))

# Add test statistics to the table
table_stats <- rbind(table_stats, jbt, adf_stat)
rownames(table_stats) <- c("Minimum", "Mean", "Maximum", "Stdev", "Skewness", "Kurtosis", "Jarque-Bera test", "ADF test")

# Adjust row names and format the table
rownames(table_stats)[7] <- "SPUSBT"
rownames(table_stats)[9] <- "CL"
rownames(table_stats)[10] <- "GC"
table_stats <- as.data.frame(table_stats)
table_stats <- round(table_stats, digits = 2)
table_stats[6:10, 7] <- "-"

##### Correlation Matrix #######
m2 <- cor(ret_all[,1:5])
m2[lower.tri(m2, diag = T)] <- NA


##### Analysis of Non-Gaussian Residuals #####
# Linear models for each cryptocurrency
lm_BTC <- lm(ret.df[, 1] ~ ret.df[, 6] + ret.df[, 7] + ret.df[, 8] + ret.df[, 9] + ret.df[, 10])
BTC_stdres <- rstandard(lm_BTC)

lm_ETH <- lm(ret.df[, 2] ~ ret.df[, 6] + ret.df[, 7] + ret.df[, 8] + ret.df[, 9] + ret.df[, 10])
ETH_stdres <- rstandard(lm_ETH)

lm_LTC <- lm(ret.df[, 3] ~ ret.df[, 6] + ret.df[, 7] + ret.df[, 8] + ret.df[, 9] + ret.df[, 10])
LTC_stdres <- rstandard(lm_LTC)

lm_XRP <- lm(ret.df[, 4] ~ ret.df[, 6] + ret.df[, 7] + ret.df[, 8] + ret.df[, 9] + ret.df[, 10])
XRP_stdres <- rstandard(lm_XRP)

lm_BCH <- lm(ret.df[, 5] ~ ret.df[, 6] + ret.df[, 7] + ret.df[, 8] + ret.df[, 9] + ret.df[, 10])
BCH_stdres <- rstandard(lm_BCH)

# Prepare data for plotting
d <- data.frame(
  Group = c(rep("Bitcoin", length(BTC_stdres)), rep("Ethereum", length(ETH_stdres)), 
            rep("Litecoin", length(LTC_stdres)), rep("Ripple", length(XRP_stdres)), 
            rep("Bitcoin Cash", length(BCH_stdres))), 
  Sample = c(BTC_stdres, ETH_stdres, LTC_stdres, XRP_stdres, BCH_stdres)
)

# Plot the QQ plots of standardized residuals
ggplot(d, aes(sample = Sample, colour = as.factor(Group))) +
  stat_qq() + 
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2) + 
  labs(color = "Cryptocurrencies") +
  xlab("Normal Scores") + 
  ylab("Standardized Residuals") +
  theme(
    legend.text = element_text(size = 15), 
    legend.title = element_text(size = 15),
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28), 
    axis.text.x = element_text(size = 20), 
    axis.text.y = element_text(size = 20)
  )


###### Prices and Returns Plot ######
# Clear the workspace
rm(list = ls())

# Load required libraries
library(ggplot2)
library(reshape2)
library(scales)
library(ggpubr)

# Load returns data
load("data/returns.RData")

# Remove any rows with NA values from the returns data
ret_all <- na.omit(ret.df)

# Prepare the data for plotting returns
date <- as.Date(rownames(ret_all))
ret_plots_df <- data.frame(date, 
                           ret_all[, 1], ret_all[, 2], ret_all[, 3], ret_all[, 4], ret_all[, 5])
colnames(ret_plots_df) <- c("Date", paste(colnames(ret_all)[1:5]))

# Melt the data for ggplot
gg <- melt(ret_plots_df, id = "Date", variable.name = "Cryptocurrency", value.name = "Log.Return")


# Create the returns plot
ret_plot <- ggplot(gg, aes(x = Date, y = Log.Return, color = Cryptocurrency)) + ylab("Log Returns") +
  theme(legend.text = element_text(size = 24), legend.title = element_text(size = 24),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28), axis.text.x = element_text(size = 24), 
        axis.text.y = element_text(size = 24)) +
  geom_point(alpha = 0) + geom_line(alpha = 0.8) +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") + 
  annotate("rect", xmin = as.Date("2017-11-01"), xmax = as.Date("2018-02-01"), ymin = -Inf, ymax = Inf,
           alpha = 0.2, color = "gray") +
  annotate("rect", xmin = as.Date("2018-09-01"), xmax = as.Date("2019-01-01"), ymin = -Inf, ymax = Inf,
           alpha = 0.2, color = "gray") +
  annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2020-03-30"), ymin = -Inf, ymax = Inf,
           alpha = 0.2, color = "gray") +
  annotate("rect", xmin = as.Date("2020-11-01"), xmax = as.Date("2021-07-01"), ymin = -Inf, ymax = Inf,
           alpha = 0.2, color = "gray") 
# geom_vline(xintercept = as.numeric(date[c(106, 654, 820, 1146)]), linetype = 4, colour = "black")

# Load data prices
load("data/prices_df.RData")
# Function for normalizing time series data
normalTS <- function(x) (x - min(x)) / (max(x) - min(x))

# Normalize the prices data
prices_norm <- apply(pf.df, 2, normalTS)
date <- as.Date(rownames(pf.df))

# Prepare the data for plotting normalized prices
prices_plots_df <- data.frame(date, 
                              prices_norm[, 1], prices_norm[, 2], prices_norm[, 3], prices_norm[, 4], prices_norm[, 5])
colnames(prices_plots_df) <- c("Date", "BTC", "ETH", "LTC", "XRP", "BCH")

# Melt the data for ggplot
gg <- melt(prices_plots_df, id = "Date", variable.name = "Cryptocurrency", value.name = "Price")

# Create the prices plot
prices_plot <- ggplot(gg, aes(x = Date, y = Price, color = Cryptocurrency)) +
  theme(legend.text = element_text(size = 24), legend.title = element_text(size = 24),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28), axis.text.x = element_text(size = 24), 
        axis.text.y = element_text(size = 24)) +
  geom_point(alpha = 0) + geom_line() +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year")
# scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) #+
# geom_vline(xintercept = as.numeric(date[c(106, 654, 820, 1146)]), linetype = 4, colour = "black")

# Arrange the plots in a single figure
ggarrange(prices_plot, ret_plot, nrow = 2, ncol = 1, common.legend = TRUE, legend = "top")

png("prices+ret.png", width = 1414, height = 549)
#labels = c("A", "B"), font.label = list(size = 15, color = "red", face = "bold"))
#w:1414, h:549
print(prices_plot)
dev.off()
