n <- 40
x <- rnorm(n, mean = 0, sd = 4)
gamma <- 0.3
c0 <- log(gamma / (1 - gamma))
c1 <- log((1-gamma) / gamma)

expit <- function(y) {
  return(exp(y) / (1 + exp(y)))
}
logit <- function(z) {
  return(log(z / (1-z)))
}

f1 <- function(x, cx) {
  n <- length(x)
  factor1 <- 1 / (n * mean(x) * (1-mean(x)))
  factor2 <- (expit(logit(mean(x)) - cx) / (1 + exp(logit(mean(x)) - cx)))^2
  return(factor1 * factor2)
}
f2 <- function(x, cx) {
  n <- length(x)
  num <- mean(x) * (1 - mean(x)) * exp(2*cx)
  den <- n * ((1-mean(x)) * exp(cx) + mean(x))^4
  return(num / den)
}

# Try it out
f1(x, c0)
f2(x, c0)
f1(x, c1)
f2(x, c1)
