# Exercise 2.8 Solution

## Problem Statement

Each time a person shops at a grocery store, the probability of catching a cold or virus is constant at \(0.01\), independent from visit to visit.

Let:
  
  - \( p = 0.01 \), the probability of catching a cold on a single visit.

We aim to calculate the probability of catching a cold at least once over \( n \) grocery visits.

## Solution

### Step 1: Probability of Not Catching a Cold on a Single Visit

Since the probability of catching a cold on a single visit is \( p = 0.01 \), the probability of not catching a cold on a single visit is:
  
  
  \[
    1 - p = 1 - 0.01 = 0.99
    \]

### Step 2: Probability of Not Catching a Cold Over \( n \) Visits

The probability of not catching a cold over \( n \) independent visits is given by:
  
  \[
    (1 - p)^n
    \]

### Step 3: Probability of Catching a Cold at Least Once Over \( n \) Visits

Using the complement rule, the probability of catching a cold at least once over \( n \) visits is:
  
  \[
    1 - (1 - p)^n
    \]

## Example Calculation

Let’s calculate the probability of catching a cold at least once if a person visits the grocery store \( n = 10 \), \( n = 50 \), and \( n = 100 \) times.

```{r}
# Define probability and number of visits
p <- 0.01
n_values <- c(10, 50, 100)

# Calculate probability of catching a cold at least once
prob_cold <- 1 - (1 - p)^n_values
names(prob_cold) <- paste("n =", n_values)
prob_cold
```

The results provide the probabilities for each value of \( n \).

## Conclusion

As the number of visits increases, the probability of catching a cold at least once also increases, approaching certainty.

## Exercise 2.16 Solution

### Problem Statement

A hospital records the daily number of people who come to the emergency room. We analyze two parts:
  
  - **(a)** Daily admissions from Sunday to Saturday are 10, 8, 14, 7, 21, 44, and 60. Assess whether a Poisson distribution can adequately model this data.
- **(b)** Discuss whether a Poisson model could better describe weekly admissions for a rare disease.

### Solution

#### Part (a): Daily Emergency Room Visits Analysis

The Poisson distribution describes events occurring independently in a fixed time with a constant mean rate. For adequacy, the mean and variance should be roughly equal.

Given data:
  
  \[
    \text{Observations} = 10, 8, 14, 7, 21, 44, 60
    \]

Calculating mean and variance:
  
  \[
    \text{Mean} = \frac{10 + 8 + 14 + 7 + 21 + 44 + 60}{7} = \frac{164}{7} \approx 23.43
    \]

\[
  \text{Variance} = \frac{(10 - 23.43)^2 + (8 - 23.43)^2 + (14 - 23.43)^2 + (7 - 23.43)^2 + (21 - 23.43)^2 + (44 - 23.43)^2 + (60 - 23.43)^2}{6} \approx 366.57
  \]

The variance (\( \approx 366.57 \)) is much greater than the mean (\( \approx 23.43 \)), showing high variability. Since Poisson distributions expect mean and variance to be close, this data likely does not follow a Poisson distribution well.

#### Part (b): Poisson Model for Weekly Rare Disease Admissions

For rare disease weekly admissions, the Poisson distribution may be suitable because:
  
  - Rare admissions occur infrequently and independently.
- With low expected rates, the Poisson distribution effectively captures sparse data, estimating the probability of few admissions in any given week.

Thus, the Poisson distribution may better describe rare disease admissions than daily emergency room visits in (a).

## Exercise 2.21 Solution

### Problem Statement

Plot the gamma distribution by fixing the shape parameter \( k = 3 \) and setting the scale parameter \( \theta = 0.5, 1, 2, 3, 4, 5 \). What is the effect of increasing the scale parameter?
  
  ### Solution
  
  The gamma distribution with shape parameter \( k \) and scale parameter \( \theta \) has the probability density function (PDF) given by:
  
  \[
    f(x; k, \theta) = \frac{x^{k-1} e^{-\frac{x}{\theta}}}{\theta^k \Gamma(k)}
    \]

where \( x \geq 0 \), \( k > 0 \), \( \theta > 0 \), and \( \Gamma(k) \) is the gamma function.

In this solution, we set \( k = 3 \) and vary the scale parameter \( \theta \) over the values \( \theta = 0.5, 1, 2, 3, 4, 5 \). Below is the plot of the PDF for each scale value.

```{r}
# Load necessary library
library(ggplot2)

# Define parameters
k <- 3
theta_values <- c(0.5, 1, 2, 3, 4, 5)
x <- seq(0, 20, length.out = 500)

# Create a data frame for plotting
gamma_data <- data.frame(
  x = rep(x, times = length(theta_values)),
  density = unlist(lapply(theta_values, function(theta) dgamma(x, shape = k, scale = theta))),
  theta = factor(rep(theta_values, each = length(x)))
)

# Plotting the gamma distributions
ggplot(gamma_data, aes(x = x, y = density, color = theta)) +
  geom_line() +
  labs(title = "Gamma Distribution with Shape Parameter k = 3 and Various Scale Parameters θ",
       x = "x", y = "Density", color = "Scale (θ)") +
  theme_minimal()
```
## Exercise 2.26 Solution

### Problem Statement

Refer to Table 2.4 cross-classifying happiness with family income.

| **Relative Family Income** | **Not Too Happy** | **Pretty Happy** | **Very Happy** | **Total** |
  |----------------------------|-------------------|------------------|----------------|-----------|
  | Below Average              | 0.080            | 0.198           | 0.079         | 0.357     |
  | Average                    | 0.043            | 0.254           | 0.143         | 0.440     |
  | Above Average              | 0.017            | 0.105           | 0.081         | 0.203     |
  | **Total**                  | 0.140            | 0.557           | 0.303         | 1.000     |
  
  - **(a)** Find and interpret the correlation using scores (i) (1,2,3) for each variable, and (ii) (1,2,3) for family income and (1,4,5) for happiness.
- **(b)** Construct the joint distribution that has these marginal distributions and exhibits independence of \( X \) and \( Y \).

### Solution

#### (a) Finding and Interpreting the Correlation

To compute the correlation between **Happiness** (Y) and **Relative Family Income** (X), we need to:
  
  1. Assign scores to each variable as given:
  - For case (i): \( X = (1, 2, 3) \) and \( Y = (1, 2, 3) \)
- For case (ii): \( X = (1, 2, 3) \) and \( Y = (1, 4, 5) \)
2. Use the formula for correlation:
  
  \[
    \text{Correlation} = \frac{\sum (X_i - \bar{X})(Y_j - \bar{Y})P(X_i, Y_j)}{\sigma_X \sigma_Y}
    \]

where:
  - \( \bar{X} \) and \( \bar{Y} \) are the means of \( X \) and \( Y \),
- \( \sigma_X \) and \( \sigma_Y \) are the standard deviations of \( X \) and \( Y \),
- \( P(X_i, Y_j) \) represents the joint probabilities from the table.

Using this formula, we calculate the following results:
  
  - **With scores (1, 2, 3) for both variables**:
  
  \[
    \text{Correlation} = 0.191
    \]

This indicates a mild positive correlation, meaning that as family income increases, happiness tends to increase slightly.

- **With scores (1, 2, 3) for family income and (1, 4, 5) for happiness**:
  
  \[
    \text{Correlation} = 0.190
    \]

This also indicates a mild positive correlation, even with different scores for happiness. The interpretation remains similar: higher family income is associated with greater happiness, though the relationship is not very strong.

#### (b) Constructing the Joint Distribution under Independence

To construct a joint distribution assuming independence of \( X \) and \( Y \), we calculate \( P(X_i, Y_j) \) as the product of the marginal probabilities:
  
  \[
    P(X_i, Y_j) = P(X_i) \cdot P(Y_j)
    \]

Using the marginal distributions from the table, we compute each entry:
  
  ```{r}
# Define marginal probabilities
marginal_income <- c(0.357, 0.440, 0.203)
marginal_happiness <- c(0.140, 0.557, 0.303)

# Compute joint probabilities under independence
joint_distribution <- outer(marginal_income, marginal_happiness)
joint_distribution
```
## Exercise 2.52 Solution

### Problem Statement

The probability density function (pdf) \( f \) of a \( N(\mu, \sigma^2) \) distribution can be derived from the standard normal pdf \( \phi \).

- **(a)** Show that the normal cumulative distribution function (cdf) \( F \) relates to the standard normal cdf \( \Phi \) by \( F(y) = \Phi \left( \frac{y - \mu}{\sigma} \right) \).
- **(b)** From (a), show that \( f(y) = \frac{1}{\sigma} \phi \left( \frac{y - \mu}{\sigma} \right) \).

### Solution

#### (a) Showing that \( F(y) = \Phi \left( \frac{y - \mu}{\sigma} \right) \)

The cumulative distribution function (cdf) \( F(y) \) of a normal random variable \( Y \sim N(\mu, \sigma^2) \) is given by:
  
  \[
    F(y) = P(Y \leq y).
    \]

We can standardize the variable \( Y \) by rewriting it in terms of a standard normal variable \( Z \), where \( Z = \frac{Y - \mu}{\sigma} \). Then we have:
  
  \[
    F(y) = P(Y \leq y) = P \left( \frac{Y - \mu}{\sigma} \leq \frac{y - \mu}{\sigma} \right).
    \]

Since \( Z \sim N(0, 1) \), we can write the probability in terms of the cdf \( \Phi \) of the standard normal distribution:
  
  \[
    F(y) = \Phi \left( \frac{y - \mu}{\sigma} \right).
    \]

This completes the solution for part (a), showing that:
  
  \[
    F(y) = \Phi \left( \frac{y - \mu}{\sigma} \right).
    \]

#### (b) Showing that \( f(y) = \frac{1}{\sigma} \phi \left( \frac{y - \mu}{\sigma} \right) \)

The probability density function \( f(y) \) of a normal random variable \( Y \sim N(\mu, \sigma^2) \) can be derived from the pdf \( \phi(z) \) of the standard normal distribution.

The standard normal pdf \( \phi(z) \) is given by:
  
  \[
    \phi(z) = \frac{1}{\sqrt{2 \pi}} e^{-z^2 / 2}.
    \]

To find \( f(y) \), we use the transformation \( Z = \frac{Y - \mu}{\sigma} \), which implies \( Y = \mu + \sigma Z \). The pdf \( f(y) \) can be found using the change of variables method:
  
  \[
    f(y) = \frac{d}{dy} \left( P(Y \leq y) \right) = \frac{d}{dy} \left( \Phi \left( \frac{y - \mu}{\sigma} \right) \right).
    \]

Applying the chain rule, we get:
  
  \[
    f(y) = \Phi' \left( \frac{y - \mu}{\sigma} \right) \cdot \frac{d}{dy} \left( \frac{y - \mu}{\sigma} \right).
\]

Since \( \Phi'(z) = \phi(z) \) (the derivative of the standard normal cdf is the standard normal pdf), we have:
  
  \[
    f(y) = \phi \left( \frac{y - \mu}{\sigma} \right) \cdot \frac{1}{\sigma}.
    \]

Substituting \( \phi \left( \frac{y - \mu}{\sigma} \right) = \frac{1}{\sqrt{2 \pi}} e^{-(y - \mu)^2 / (2 \sigma^2)} \), we get:
  
  \[
    f(y) = \frac{1}{\sigma} \cdot \frac{1}{\sqrt{2 \pi}} e^{-(y - \mu)^2 / (2 \sigma^2)}.
    \]

Therefore,

\[
  f(y) = \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(y - \mu)^2 / (2 \sigma^2)},
  \]

which matches the form given in equation (2.8):
  
  \[
    f(y; \mu, \sigma) = \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(y - \mu)^2 / (2 \sigma^2)}.
    
    ## Exercise 2.53 Solution
    
    ### Problem Statement
    
    If \( Y \) is a standard normal random variable, with cdf \( \Phi \), what is the probability distribution of \( X = \Phi(Y) \)? Illustrate by randomly generating a million standard normal random variables, applying the \( \text{cdf} \) function \( \Phi() \) to each, and plotting histograms of the (a) \( y \) values and (b) \( x \) values.
    
    ### Solution
    
    #### Understanding the Distribution of \( X = \Phi(Y) \)
    
    Let \( Y \sim N(0,1) \) be a standard normal random variable, where \( \Phi \) is the cumulative distribution function (CDF) of \( Y \). The transformation \( X = \Phi(Y) \) applies the standard normal CDF to \( Y \), which maps the values of \( Y \) to the interval \([0, 1]\).
    
    Since the CDF of a random variable transforms it into a uniform distribution over \([0, 1]\), \( X = \Phi(Y) \) follows a **uniform distribution** on the interval \([0, 1]\). This is a well-known result: any continuous random variable, when transformed by its own CDF, results in a uniform distribution over \([0,1]\).
    
    #### Illustration with a Simulation
    
    To illustrate this, we can:
      
      1. Generate a large number of samples (e.g., 1 million) from a standard normal distribution for \( Y \).
    2. Compute \( X = \Phi(Y) \) for each sample, where \( \Phi \) is the CDF of the standard normal distribution.
    3. Plot histograms of the generated values of \( Y \) and \( X \) to verify their distributions.
    
    #### R Code and Histograms
    
    The following R code generates 1 million samples from a standard normal distribution, applies the CDF to each sample, and plots histograms of both \( Y \) and \( X \).
    
    ```{r}
    # Load necessary library
    library(ggplot2)
    
    # Number of samples
    num_samples <- 1e6
    
    # Step 1: Generate standard normal random variables for Y
    Y <- rnorm(num_samples)
    
    # Step 2: Apply the CDF function to each Y to get X = Φ(Y)
    X <- pnorm(Y)
    
    # Step 3: Plot histograms for Y and X
    par(mfrow = c(1, 2))
    
    # Histogram of Y (Standard Normal Distribution)
    hist(Y, breaks = 50, col = "skyblue", main = "Histogram of Y (Standard Normal Distribution)", xlab = "Y values", freq = FALSE)
    
    # Histogram of X (Uniform Distribution on [0,1])
    hist(X, breaks = 50, col = "lightcoral", main = "Histogram of X = Φ(Y) (Uniform Distribution on [0,1])", xlab = "X values", freq = FALSE)
    ```
    ## Exercise 2.70 Solution
    
    ### Problem Statement
    
    The **beta distribution** is a probability distribution over \( (0, 1) \) that is often used in applications for which the random variable is a proportion. The beta pdf is
    
    \[
      f(y; \alpha, \beta) = \frac{\Gamma(\alpha + \beta)}{\Gamma(\alpha)\Gamma(\beta)} y^{\alpha - 1}(1 - y)^{\beta - 1}, \quad 0 \leq y \leq 1,
      \]
    
    for parameters \( \alpha \) and \( \beta \), where \( \Gamma(\cdot) \) denotes the gamma function.
    
    ### Solution
    
    #### (a) Showing that the Uniform Distribution is the Special Case \( \alpha = \beta = 1 \)
    
    When \( \alpha = \beta = 1 \), the beta pdf becomes
    
    \[
      f(y; 1, 1) = \frac{\Gamma(2)}{\Gamma(1)\Gamma(1)} y^{1 - 1}(1 - y)^{1 - 1} = \frac{\Gamma(2)}{\Gamma(1)\Gamma(1)} = 1, \quad 0 \leq y \leq 1.
      \]
    
    Since \( f(y; 1, 1) = 1 \) over \( [0, 1] \), this is the pdf of the uniform distribution over \( [0, 1] \). Thus, the beta distribution reduces to the uniform distribution when \( \alpha = \beta = 1 \).
    
    #### (b) Showing that \( \mu = E(Y) = \frac{\alpha}{\alpha + \beta} \)
    
    For the beta distribution, the mean \( E(Y) \) is known to be
    
    \[
      E(Y) = \frac{\alpha}{\alpha + \beta}.
      \]
    
    This can be derived by using the properties of the beta distribution, specifically that the mean of a Beta(\( \alpha, \beta \)) distribution is given by \( \frac{\alpha}{\alpha + \beta} \). Therefore, we have shown that \( \mu = \frac{\alpha}{\alpha + \beta} \).
    
    #### (c) Finding \( E(Y^2) \) and Showing that \( \text{var}(Y) = \frac{\alpha \beta}{(\alpha + \beta)^2(\alpha + \beta + 1)} = \frac{\mu (1 - \mu)}{\alpha + \beta + 1} \)
    
    The second moment \( E(Y^2) \) for the beta distribution is given by
    
    \[
      E(Y^2) = \frac{\alpha (\alpha + 1)}{(\alpha + \beta)(\alpha + \beta + 1)}.
      \]
    
    Now, the variance \( \text{var}(Y) \) is computed as
    
    \[
      \text{var}(Y) = E(Y^2) - (E(Y))^2.
      \]
    
    Substituting \( E(Y) = \frac{\alpha}{\alpha + \beta} \) and \( E(Y^2) = \frac{\alpha (\alpha + 1)}{(\alpha + \beta)(\alpha + \beta + 1)} \), we get
    
    \[
      \text{var}(Y) = \frac{\alpha (\alpha + 1)}{(\alpha + \beta)(\alpha + \beta + 1)} - \left( \frac{\alpha}{\alpha + \beta} \right)^2.
      \]
    
    Simplifying this expression yields
    
    \[
      \text{var}(Y) = \frac{\alpha \beta}{(\alpha + \beta)^2(\alpha + \beta + 1)}.
      \]
    
    Alternatively, using \( \mu = \frac{\alpha}{\alpha + \beta} \), we can rewrite this as
    
    \[
      \text{var}(Y) = \frac{\mu (1 - \mu)}{\alpha + \beta + 1}.
      \]
    
    #### (d) Plotting the Beta PDF in R
    
    In R, we can use the following code to plot the beta pdf for the specified values of \( \alpha \) and \( \beta \):
      
      ```{r}
    # Define a sequence of y values
    y <- seq(0, 1, length.out = 100)
    
    # Plot for alpha = beta = 0.5, 1, 10, 100
    plot(y, dbeta(y, 0.5, 0.5), type = "l", col = "blue", ylim = c(0, 5),
         main = "Beta PDF for Various Values of alpha and beta",
         ylab = "Density", xlab = "y")
    lines(y, dbeta(y, 1, 1), col = "red")
    lines(y, dbeta(y, 10, 10), col = "green")
    lines(y, dbeta(y, 100, 100), col = "purple")
    legend("topright", legend = c("alpha = beta = 0.5", "alpha = beta = 1",
                                  "alpha = beta = 10", "alpha = beta = 100"),
           col = c("blue", "red", "green", "purple"), lty = 1)
    
    # Plot for some values of alpha > beta and alpha < beta
    plot(y, dbeta(y, 2, 5), type = "l", col = "blue", ylim = c(0, 3),
         main = "Beta PDF for alpha > beta and alpha < beta",
         ylab = "Density", xlab = "y")
    lines(y, dbeta(y, 5, 2), col = "red")
    legend("topright", legend = c("alpha = 2, beta = 5", "alpha = 5, beta = 2"),
           col = c("blue", "red"), lty = 1)
    ```