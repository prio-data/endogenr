# Set random seed for reproducibility
set.seed(123)

# Sample size
n <- 1000

# Generate predictors
x1 <- rnorm(n, mean = 2, sd = 1)
x2 <- rnorm(n, mean = 0, sd = 2)

# True parameters
beta0 <- 1.5  # intercept
beta1 <- 2.0  # coefficient for x1
beta2 <- -0.5 # coefficient for x2

# Generate error term
epsilon <- rnorm(n, mean = 0, sd = 1)

# Generate dependent variable
y <- beta0 + beta1*x1 + beta2*x2 + epsilon

# Create data frame
data <- data.frame(y = y, x1 = x1, x2 = x2)
train <- data[1:500,]
test <- data[501:1000,]

# Fit linear model
model <- lm(y ~ x1 + x2, data = train)

# View results
summary(model)

lmpred <- predict(model, newdata = test, se.fit = T)

test$pred <- getpi(lmpred)
test$ae <- test$y - test$pred


predict(model, newdata = test, interval = "prediction") |> head()
# These should approximate the values above
pi <- getpi(lmpred, 1000)
pi[1,] |> quantile(probs = c(0.05, 0.975))
pi[2,] |> quantile(probs = c(0.05, 0.975))
pi[3,] |> quantile(probs = c(0.05, 0.975))
pi[4,] |> quantile(probs = c(0.05, 0.975))
pi[5,] |> quantile(probs = c(0.05, 0.975))
