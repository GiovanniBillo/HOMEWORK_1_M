## Homework 1
#### TO USE GITHUB ####
install.packages("usethis")
library(usethis)

# FSDS: Chapter 2, exercises 2.8, 2.16, 2.21, 2.26, 2.52, 2.53, 2.70

# 2.8
#a
# The probabilities of catching a virus at  each single visit are 
# independent from the previous visits, they don;t sum up
#b
# assuming the question means the probability of catching a virus AT LEAST ONCE
dbinom(1, 100, 0.01) # with "success" = catching the virus

# 2.16
obs = c(10, 8, 14, 7, 21, 44, 60)
mean(obs)
dpois(10, mean(obs))
#a
# The Poisson distribution might not model the phenomenon adequately
# because there is a non-negligeable probability that multiple
# events occur simoultaneously ( e.g car accident with multiple people 
# being taken to the emergency room)
# Also, variance is high, and there seems to be a hike in the weekends
#b
# The number of weekly admissions to a hospital for a rare disease could
# be better modelled by a poisson distribution, as the events are independent
# and (practically) cannot occur simoultaneously

#2.21
# Parameters for the gamma distribution
shape <- 3    # Shape parameter
scale = c(0.5, 1, 2, 3, 4, 5)

# Generate a color for each curve
colors <- rainbow(length(scale))  # Generate a unique color for each curve

# Initial plot with the first gamma curve, setting up empty plot area
y <- dgamma(x, shape = scale[1], rate = rate)
plot(x, y, type = "n", ylim = c(0, 0.5),
     main = "Gamma Distributions with Varying Shape Parameters",
     xlab = "x", ylab = "Density")

# Loop through each shape parameter, adding a curve with a unique color
for (i in 1:length(scale)) {
  y <- dgamma(x, shape = shape, rate = 1/scale[i])
  lines(x, y, col = colors[i], lwd = 2)  # Plot each line with a unique color
}

# Add a legend with corresponding colors
legend("topright", legend = paste("shape =", shape_values),
       col = colors, lwd = 2, cex = 0.8)

# As we increase the scale parameter, the distribution becomes much more sparse


# 2.26
#a
# i)
probabilities <- c(0.080, 0.043, 0.017, 0.198, 0.254, 0.105, 0.079, 0.143, 0.081)
x <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
y <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)
install.packages("wCorr")
library(wCorr)
# after entering install.packages("wCorr")
corr_1 = weightedCorr(x, y, weights=probabilities)

#INTERPRETATION:  by assigning equal weights 1, 2, 3 to family income and happiness
# there seems to be a weak positive correlation between the variables

# ii)
probabilities <- c(0.080, 0.043, 0.017, 0.198, 0.254, 0.105, 0.079, 0.143, 0.081)
x <- c(1, 1, 1, 4, 4, 4, 5, 5, 5)
y <- c(1, 4, 5, 1, 4, 5, 1, 4, 5)
library(wCorr)
corr_2 = weightedCorr(x, y, weights=probabilities)

#INTERPRETATION:  results are like the previous example, only now the correlation 
# is slightly larger by assigning larger scores to "pretty happy" and "very happy"

#b
#??? idk
  
# 2.52
  
  
# 2.53
R = 1000000
Y = rnorm(R)
X = pnorm(R)

hist(X) #uniform distribution
hist(Y)

# 2.70
# b & c missing
###################################################
# Define parameters for alpha = beta cases
alpha = beta = c(0.5, 1, 10, 100)

# Generate a unique color for each curve
colors <- rainbow(length(alpha))

# Generate x-values
x <- seq(0, 1, length.out = 100)

# Set up the plot with no lines initially
plot(x, dbeta(x, alpha[1], beta[1]), type = "n",
     main = "Beta Distributions with same alpha-beta Parameters",
     xlab = "x", ylab = "Density")

# Loop through each shape parameter, add each curve
for (i in 1:length(alpha)) {
  # Calculate the beta PDF for each alpha = beta value
  y <- dbeta(x, shape1 = alpha[i], shape2 = beta[i])
  
  # Plot each line with a unique color
  lines(x, y, col = colors[i], lwd = 2)
}

# Add a legend
legend("topright", legend = paste("alpha = beta =", alpha),
       col = colors, lwd = 2, cex = 0.8)

###############################################
# with alpha > beta
# Define parameters for alpha = beta cases
alpha = c(0.5, 1, 10, 100)
beta = c(0.3, 0.4, 5, 50)

# Generate a unique color for each curve
colors <- rainbow(length(alpha))

# Generate x-values
x <- seq(0, 1, length.out = 100)

# Set up the plot with no lines initially
plot(x, dbeta(x, alpha[1], beta[1]), type = "n",
     main = "Beta Distributions with alpha > beta Parameters",
     xlab = "x", ylab = "Density")

# Loop through each shape parameter, add each curve
for (i in 1:length(alpha)) {
  # Calculate the beta PDF for each alpha = beta value
  y <- dbeta(x, shape1 = alpha[i], shape2 = beta[i])
  
  # Plot each line with a unique color
  lines(x, y, col = colors[i], lwd = 2)
}

# Add a legend
legend("topright", legend = paste("alpha = ", alpha, "beta = ", beta),
       col = colors, lwd = 2, cex = 0.8)


#########################################################
# with alpha < beta
# Define parameters for alpha = beta cases
alpha = c(0.5, 1, 10, 100)
beta = c(0.7, 2, 20, 150)


# Generate a unique color for each curve
colors <- rainbow(length(alpha))

# Generate x-values
x <- seq(0, 1, length.out = 100)

# Set up the plot with no lines initially
plot(x, dbeta(x, alpha[1], beta[1]), type = "n",
     main = "Beta Distributions with alpha < beta Parameters",
     xlab = "x", ylab = "Density")

# Loop through each shape parameter, add each curve
for (i in 1:length(alpha)) {
  # Calculate the beta PDF for each alpha = beta value
  y <- dbeta(x, shape1 = alpha[i], shape2 = beta[i])
  
  # Plot each line with a unique color
  lines(x, y, col = colors[i], lwd = 2)
}

# Add a legend
legend("topright", legend = paste("alpha = ", alpha, "beta = ", beta),
       col = colors, lwd = 2, cex = 0.8)

# when alpha = beta, the distributions are centered and 
# the tails get fatter as the parameters decrease

# when alpha > beta the distributions are skewed to the right, with the corresponding tail being fatter

# when alpha < beta the distributions are skewed to the left, with the corresponding tail being fatter

######################################################################################################
# FSDS: Chapter 3, exercises 3.18, 3.28, 3.24 (use R)
##### 3.18
n = 90000
mean_age = 72
sd = 12

n_sample = 100
sample_mean_age = 70 
sample_sd = 11

# having a large n for both population and sample, we can assume convergence to the normal distribution
#a
X_population = rnorm(n, mean = mean_age, sd = sd)
mean(X_population)
sd(X_population)
hist(X_population)


# the sample distribution will probably be more skewed, due to the smaller number of individuals
# and the high standard deviation.

#b
X_sample = sample(X_population, size = n_sample)
mean(X_sample)
sd(X_sample)
hist(X_sample)

# the center is actually more right skewed, and describes the larger number of "younger" senior citizens that lives 
# in Sunshine City. There are not as many people older than 80 as there are younger than 70 ( for biological reasons).

#c
# IT would not be unusual to sample a person whose age is 60 because it is within 1 standard deviation from the population mean (72), 
# and in a normally distributed RV the interval within 1 sd from the mean contains approximately 68% of the variance.
# However getting a sample with mean = 60 would be highly unusual because it would mean that we sampled a group of much younger individuals
# than the actual population, which is very unlikely if the sample is truly random.

#d
hist(sample(X_population, size = 1))
# not surprisingly, the sampling distribution uniform, consisting of just one value
hist(sample(X_population, n))
# now instead the sample is equal to the size of the population and is symmetric
#