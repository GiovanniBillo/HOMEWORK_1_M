## Homework 1
#### TO USE GITHUB ####
install.packages("usethis")
library(usethis)

#### other libraries ####
library(dplyr)

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



#########################################################################
# FSDS: Chapter 4, exercises 4.14, 4.16, 4.48
#### 4.14 ####
data = read.table("https://stat4ds.rwth-aachen.de/data/Students.dat", header = TRUE)
summary(data)
#a
x_bar = mean(data$tv)
s = sd(data$tv)
# H0: the mean is 7.2
# H1: the mean is different from 7.2
z = qnorm(0.975) # having n > 60, we can assume normality
SE = s/sqrt(length(data$tv))
CI = x_bar + c(-1, 1)*z*SE
# we can say that 95% of the students on average spend between 5.30 and almost 9 hours watching TV per week

#b
only_male = data %>% filter(gender == 0)
only_female = data %>% filter(gender == 1)

boxplot(only_male$tv, only_female$tv, names = c("Males", "Females"), ylab = "weekly hours watching tv")

# assuming
# -independent populations
# -equal variances
# - normality -> weakest assumption, barely 30 observations
# --> use t student's distribution

x_bar_male = mean(only_male$tv)
x_bar_female = mean(only_female$tv)
s_male = sd(only_male$tv)
s_female = sd(only_female$tv)
n_male = length(only_male$tv)
n_female = length(only_female$tv)

# assuming equal variance
s_p = ((n_male -1)*(s_male^2) + (n_female - 1)*(s_female^2))/(n_male + n_female -1)
# t = (x_bar_male - x_bar_female)/(s_p*(1/n_male + 1/n_female))
SE_diff = (x_bar_male - x_bar_female)/((s_p^2)*(1/n_male + 1/n_female))
mean_diff = x_bar_male - x_bar_female
z = qt(0.975, n_male + n_female - 2)
CI_diff = mean_diff + c(-1, 1)*z*SE_diff

# Test for equality of means
t.test(only_male$tv, only_female$tv, conf.level = 0.05)
# we can not state that the variances are significantly different betweent the 2 groups
# However the confidence interval seems to suggest that females watch slightly more tv than males

# assuming unequal variance
SE_diff = sqrt((s_male^2)/n_male + (s_female^2)/n_female)
# t = (x_bar_male - x_bar_female)/(s_p*(1/n_male + 1/n_female))
mean_diff = x_bar_male - x_bar_female
z = qt(0.975, n_male + n_female - 2)
CI_diff_uvars = mean_diff + c(-1, 1)*z*SE_diff

# now we the interval crosses 0?? how is such a different result possible?

#### 4.16 ####
data = read.table("https://stat4ds.rwth-aachen.de/data/Substance.dat", header = TRUE)

# compare alcohol users and non-users
alpha = 0.05
n = sum(data$count)
# find the total number of students that have or haven't used alcohol
N_alcohol_total = sum(data[data$alcohol == "yes", 4])
N_NOalcohol_total = sum(data[data$alcohol == "no", 4])
N_marijuana = sum(data[data$alcohol == "yes" & data$marijuana == "yes", 4])
N_NOalcohol_marijuana = sum(data[data$alcohol == "no" & data$marijuana == "yes", 4])

pi_1_hat = N_marijuana/N_alcohol_total
pi_2_hat = N_NOalcohol_marijuana/N_NOalcohol_total

# by hand
z = qnorm(1-alpha/2)
SE = (pi_1_hat*(1 - pi_1_hat)/N_alcohol_total + pi_2_hat*(1 - pi_2_hat)/N_NOalcohol_total)
CI_prop = (pi_1_hat - pi_2_hat) + c(-1, 1)*z*sqrt(SE)

#with software
prop.test(c(955, 5), c(1949, 327), conf.level = alpha, correct = F)

## interpretation: there seems to be a significant mean difference in the use of marijuana 
## between students that used alcohol before and students who didn't.


###################
#                 #
# cs  es3.5; 3.6  #
#                 #
#                 #
###################

#Chapter 3 es 3.5

#Ax = y
set.seed(0)
n <- 1000
A <- matrix(runif(n * n), n, n)
x.true <- runif(n)
y <- A %*% x.true

#b)

start_time <- Sys.time()
A_inv <- solve(A)
x1 <- A_inv %*% y
end_time <- Sys.time()
time_taken <- end_time - start_time

# Calcolare la differenza assoluta media
mean_absolute_error <- mean(abs(x1 - x.true))

print(time_taken)
print(mean_absolute_error)

#c)
start_time_direct <- Sys.time()
x2 <- solve(A, y)
end_time_direct <- Sys.time()
time_taken_direct <- end_time_direct - start_time_direct

# Calcolare la differenza assoluta media
mean_absolute_error_direct <- mean(abs(x2 - x.true))

print(time_taken_direct)
print(mean_absolute_error_direct)

# d) Conclusione:

# Conclusion:

# Calculate time to explicitly form the inverse matrix A^-1 and solve Ax = y
start_time <- Sys.time()
A_inv <- solve(A)
x1 <- A_inv %*% y
end_time <- Sys.time()
time_taken <- end_time - start_time

# Calculate the mean absolute error between x1 and x.true
mean_absolute_error <- mean(abs(x1 - x.true))

print(time_taken)
print(mean_absolute_error)

# Calculate time to directly solve Ax = y without forming A^-1
start_time_direct <- Sys.time()
x2 <- solve(A, y)
end_time_direct <- Sys.time()
time_taken_direct <- end_time_direct - start_time_direct

# Calculate the mean absolute error between x2 and x.true
mean_absolute_error_direct <- mean(abs(x2 - x.true))

print(time_taken_direct)
print(mean_absolute_error_direct)

# Conclusion:
# Using 'solve' to directly solve the equation Ax = y is significantly faster (approximately 0.39 seconds) 
# compared to explicitly forming A^-1 and then multiplying it by y (approximately 2.29 seconds).

# While both solutions are very precise, directly solving with 'solve' yields a slightly smaller mean absolute error 
# (approximately 1.36e-12) compared to forming the inverse (approximately 2.96e-11).

# Therefore, solving the linear system directly without calculating the explicit inverse is preferable 
# in terms of both efficiency and accuracy.


#es 3.6

#a)
# Funzione per calcolare la ECDF
ecdf_values <- function(x) {
  # Numero totale di osservazioni
  n <- length(x)
  
  # Ordinare i valori di x con sort.int mantenendo gli indici originali
  sorted_x <- sort.int(x, index.return = TRUE)
  
  # Calcolare la ECDF per ogni valore in x
  ecdf <- (1:n) / n
  
  # Riordinare i valori della ECDF nell'ordine originale
  ecdf_original_order <- ecdf[order(sorted_x$ix)]
  
  return(ecdf_original_order)
}

# Test della funzione con un esempio
set.seed(123)
x <- rnorm(10)
ecdf_values(x)


#b)
# Funzione per calcolare la ECDF e opzionalmente tracciarla
ecdf_values <- function(x, plot.cdf = FALSE) {
  # Numero totale di osservazioni
  n <- length(x)
  
  # Ordinare i valori di x
  sorted_x <- sort(x)
  
  # Calcolare la ECDF per ogni valore in x
  ecdf <- sapply(x, function(xi) sum(sorted_x < xi) / n)
  
  # Se plot.cdf è TRUE, tracciare la ECDF
  if (plot.cdf) {
    plot(sort(x), ecdf[order(x)], type = "s", main = "Empirical Cumulative Distribution Function",
         xlab = "x", ylab = "ECDF")
  }
  
  return(ecdf)
}

# Test della funzione con un esempio
set.seed(123)
x <- rnorm(10)
ecdf_values(x, plot.cdf = TRUE)

# Generare campioni casuali
set.seed(123)
x_norm <- rnorm(100)
x_unif <- runif(100)

# Testare la funzione con distribuzione normale
ecdf_norm <- ecdf_values(x_norm, plot.cdf = TRUE)

# Testare la funzione con distribuzione uniforme
ecdf_unif <- ecdf_values(x_unif, plot.cdf = TRUE)







