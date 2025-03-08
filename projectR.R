library(rstan);
library("bayesplot")
library("rstanarm")
library("ggplot2")
library("mcmcse")
library(dplyr)
library(plyr)
library(mcmcse)
library(corrplot)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library("PerformanceAnalytics")
library(GGally)
library(viridis)

# Preparation

# 1 response (obsy - observed height of wave)
# y = Alpha * sin ( omega * t )+ error
# parameter to be estimated 

data <- read.csv("./CaiMin_544_data.csv",header = TRUE, sep = ",")#read.table("./CaiMin_544_data.csv",header=T)


# prepare data for Stan
nobs  <- nrow(data)                 # number of observation
time  <- data$t             # artificial time
obsv  <- data$observation           #observation
  
X11()

plot(time, obsv, type="b", col="purple", lwd=1, xlab="time", ylab = "Observation", main = "Sine Wave Data")
grid()
stan_data  <- list(nobs = nobs,
                   t = time,
                   y = obsv)

# run MCMC
# Read the Stan model code from a file
stan_model <- readLines("project_544-00.stan", warn = TRUE)


# Run the Bayesian analysis
fit2 <- stan(model_code = stan_model,
            data = stan_data,
            seed = 42,
            chains = 3,
            iter =15000, #15000
            cores= 3,
            warmup=5000 #7500
)

saveRDS(fit2, "project.rds")
fit <- readRDS("project.rds")


# Output credible interval
print(fit,pars=c("alpha", "gamma","sigma","omega","yt"),intervals=c(0.025, 0.975), digits=3)


#visualise
x11()
color_scheme_set("pink")
pairs(fit,pars=c("alpha","omega", "gamma","sigma"))

x11()
plot(fit,pars=c("alpha", "gamma","sigma","omega"),intervals=c(0.025, 0.975), digits=3)



# all chains combined
sampler_params <- get_sampler_params(fit, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

# each chain separately
lapply(sampler_params, summary, digits = 2)
#x11()
#lapply(extract(r, c("alpha", "gamma","sigma","omega","mu")))

#Traceplot to see convergence
x11()
traceplot(fit, pars = c("alpha", "gamma","sigma","omega"), inc_warmup = TRUE, nrow = 5)

x11()
traceplot(fit, pars = c("alpha", "gamma","sigma","omega"), inc_warmup = FALSE, nrow = 5)

x11() #mu got 100 pointssss
traceplot(fit, pars = c("mu[1]","mu[20]","mu[60]","mu[80]","mu[100]"), inc_warmup = TRUE, nrow = 5)

x11() #mu got 100 pointssss
traceplot(fit, pars = c("mu[1]"), inc_warmup = FALSE)

#correlation plot how



###Generate posterior
posterior <- as.matrix(fit)
#posterior <- as.array(fit) #separately wrt chain

# Assuming 'fit' is the result of your Stan sampling
posterior_samples <- extract(fit, pars = c("alpha", "gamma","sigma","omega","yt","y_rep"), permuted = TRUE)


plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(posterior,
           pars = c("alpha", "gamma","sigma","omega"),
           prob = 0.95) + plot_title


# Extract relevant parameters from the posterior samples
#alpha_samples <- posterior_samples$alphaplot_title <- ggtitle("Posterior distributions",
#                      "with medians and 95% intervals")
mcmc_areas(posterior,
           pars = c("alpha", "gamma","sigma","omega"),
           ) + plot_title
omega_samples <- posterior_samples$omega
gamma_samples <- posterior_samples$gamma
sigma_samples <- posterior_samples$sigma
yt_samples <- posterior_samples$yt
alpha_samples <- posterior_samples$alpha
yrep_samples <- posterior_samples$y_rep

# Number of posterior samples
n_samples <- nrow(posterior)
#mu_mean <- colMeans(mu_samples)
yt_mean <- colMeans(yt_samples)
yrep_mean <- colMeans(yrep_samples)

# Create a data frame for plotting
plot_data_yt <- data.frame(
  time = rep(data$t, each = n_samples),  # Assuming 'data$t' represents the time points
  yt_mean
)

plot_data_yrep <- data.frame(
  time = rep(data$t, each = n_samples),  # Assuming 'data$t' represents the time points
  yrep_mean
)

# Convert 'time' to a factor
plot_data_yt$time <- factor(plot_data_yt$time)

# Plotting
x11()


# Combine the data into a data frame
# Create a data frame with posterior samples
posterior_df <- data.frame(
  alpha = alpha_samples,
  omega = omega_samples,
  gamma = gamma_samples,
  sigma = sigma_samples,
  #mu = mu_samples,
  yt = yt_samples,
  time = time,
  yrep = yrep_samples
)

df <- data.frame(time = time, yt_mean, obsv = obsv,yrep=yrep_mean)

# Plot the data
ggplot(df, aes(x = time)) +
  geom_line(aes(y = yrep_mean, color = "Inference mu"), size = 1.5) +
  geom_point(aes(y = obsv, color = "Observation"), size = 3) +
  labs(x = "Time", y = "Value", title = "Sine Wave Data") +
  scale_color_manual(values = c("Inference mu" = "cyan", "Observation" = "purple")) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.box.background = element_rect(color = "black"),
        legend.box.margin = margin(0, 0, -1, 0),
        legend.margin = margin(0, 0, 10, 0)) +
  guides(color = guide_legend(override.aes = list(size = c(1.5, 3))))

# Add grid
theme_set(theme_classic())

# Quantiles of alpha
round(as.data.frame(mcmcse::mcse.q(posterior_samples$alpha, .975)),3)
round(as.data.frame(mcmcse::mcse.q(posterior_samples$alpha, .025)),3)

# Quantiles of gamma
round(as.data.frame(mcmcse::mcse.q(posterior_samples$gamma, .975)),3)
round(as.data.frame(mcmcse::mcse.q(posterior_samples$gamma, .025)),3)

# Quantiles of omega
round(as.data.frame(mcmcse::mcse.q(posterior_samples$omega, .975)),3)
round(as.data.frame(mcmcse::mcse.q(posterior_samples$omega, .025)),3)


# Quantiles of omega
round(as.data.frame(mcmcse::mcse.q(posterior_samples$sigma, .975)),3)
round(as.data.frame(mcmcse::mcse.q(posterior_samples$sigma, .025)),3)

# Quantiles of omega
round(as.data.frame(mcmcse::mcse.q(posterior_samples$mu, .975)),3)
round(as.data.frame(mcmcse::mcse.q(posterior_samples$mu, .025)),3)
x11()
plot(cumsum(posterior_samples$mu)/1:length(posterior_samples$mu),type="l",
     ylab="Posterior mean of mu", xlab="Iteration(t)")
x11()
plot(cumsum(posterior_samples$alpha)/1:length(posterior_samples$alpha),type="l",
     ylab="Posterior mean of alpha", xlab="Iteration(t)")
x11()
plot(cumsum(posterior_samples$gamma)/1:length(posterior_samples$gamma),type="l",
     ylab="Posterior mean of gamma", xlab="Iteration(t)")
x11()
plot(cumsum(posterior_samples$omega)/1:length(posterior_samples$omega),type="l",
     ylab="Posterior mean of omega", xlab="Iteration(t)")
x11()
plot(cumsum(posterior_samples$sigma)/1:length(posterior_samples$sigma),type="l",
     ylab="Posterior mean of sigma", xlab="Iteration(t)")

x11()
# Plot the data
ggplot(df, aes(x = time)) +
  geom_line(aes(y = yrep_mean, color = "Inference yrep"), size = 1.5) +
  geom_point(aes(y = obsv, color = "Observation"), size = 3) +
  labs(x = "Time", y = "Value", title = "Sine Wave Data") +
  scale_color_manual(values = c("Inference yrep" = "cyan", "Observation" = "purple")) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.box.background = element_rect(color = "black"),
        legend.box.margin = margin(0, 0, -1, 0),
        legend.margin = margin(0, 0, 10, 0)) +
  guides(color = guide_legend(override.aes = list(size = c(1.5, 3))))

#Print certain summary
summary_stats <- summary(fit, pars = c("alpha", "gamma", "sigma", "omega"))$summary[, c("mean", "sd", "2.5%", "97.5%", "n_eff","Rhat")]
summary_stats <- round(summary_stats, 2) # Round to two decimal points
print(summary_stats)

#Print certain summary
summary_stats <- summary(fit, pars = c("mu[]", "gamma", "sigma", "omega"))$summary[, c("mean", "sd", "2.5%", "97.5%", "n_eff","Rhat")]
summary_stats <- round(summary_stats, 5) # Round to two decimal points
print(summary_stats)

#print y
x11()
color_scheme_set("purple")
ppc_intervals(y = obsv, yrep = posterior_samples$y_rep, x = time,
           prob = 0.025, prob_outer = 0.975, size = 2, fatten = 2) +panel_bg(fill="gray90", color = NA)+
  grid_lines(color = "white")+ ggplot2::xlab("time(s)")+ ggplot2::ylab("y(t)")

x11()
# Create PPC ribbon plot
# Create PPC ribbon plot
ppc_ribbon(y = obsv, yrep = posterior_samples$y_rep, x = time,
           prob = 0.025, prob_outer = 0.975) +
  geom_point(aes(y = obsv)) +
  grid_lines(color = "white")+ ggplot2::xlab("time(s)")+ ggplot2::ylab("y(t)")

+scale_color_manual(values = "blue")+
  labs(x = "time(s)", y = "y(t)", title = "Posterior Distribution of y") +
  scale_color_manual(values = c("Observation" = "red")) 

x11()
# histograms of some parameters
color_scheme_set("pink")
mcmc_hist(fit, pars = c("alpha", "gamma", "sigma", "omega"))


corr_df <- data.frame(
  alpha = alpha_samples,
  omega = omega_samples,
  gamma = gamma_samples,
  sigma = sigma_samples
)

corre <- cor(corr_df)
round(corre, 2)
ggcorrplot(correlation)
x11()
corrplot(corre, type = "upper", order = "hclust", 
         tl.col = "blue", tl.srt = 45)
x11()
chart.Correlation(corr_df, histogram=TRUE, pch=19,bg= "blue", col= "blue", hist.col = "yellowgreen")
corrplot(corre, method="number")

x11()
ggpairs(corr_df)+  scale_color_viridis(discrete = TRUE)

x11()
pairs(corr_df,   # data frame used
             pch = 21,           # point symbol
             col= "blue",     # point border color
             bg= "blue",      # point background color
             digits= 2,          # number of significant digits to include
             cor = TRUE,         # TRUE reports correlations
             hist.col ="blue",# histogram color
             show.points = TRUE, # shows data pointss
             stars = TRUE,)      # statistical significant indication
