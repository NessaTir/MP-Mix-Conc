# Generate random data with a linear relationship
set.seed(123)
x <- seq(1, 100, by = 1)
y <- seq(1, 100, by = 1)
y <- x^2


# Plot the points
plot(log(x), log(y), main = "Scatter Plot with Linear Model", 
     xlab = "X values", ylab = "Y values")

plot(x, y, main = "Scatter Plot with Linear Model", 
     xlab = "X values", ylab = "Y values")


# Add a linear model
lm_model <- lm(y ~ x)
abline(lm_model, col = "red")

