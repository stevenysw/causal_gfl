# load packages
library(grf)
library(rlearner)

# auxiliary funcions
d25 <- function(x){
  sqrt(sum((x-0.25)^2))
}

d75 <- function(x){
  sqrt(sum((x-0.75)^2))
}

d_tau <- function(x){
  ans = 0
  if (min(x) >= 0.25) ans = 1
  if (max(x) < 0.25) ans = -1
  ans
}

d05 <- function(x){
  ans = 0.2
  if (sum((x-0.5)^2) < 0.8) ans = 0.8
  ans
}

d10 <- function(x){
  q1 = c(0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5)
  q2 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25, 0.25)
  q3 = c(0.75, 0.75, 0.75, 0.75, 0.75, 0.5, 0.5, 0.5, 0.5, 0.5)
  q4 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75, 0.75, 0.75)
  
  ans = -1
  
  if (sqrt(sum((x-q1)^2)) < min(sqrt(sum((x-q2)^2)), sqrt(sum((x-q3)^2)), sqrt(sum((x-q4)^2)))) ans = 2
  if (sqrt(sum((x-q2)^2)) < min(sqrt(sum((x-q1)^2)), sqrt(sum((x-q3)^2)), sqrt(sum((x-q4)^2)))) ans = 1
  if (sqrt(sum((x-q3)^2)) < min(sqrt(sum((x-q1)^2)), sqrt(sum((x-q2)^2)), sqrt(sum((x-q4)^2)))) ans = 0
  ans
}

# Scenario 1
set.seed(123)

N = 500
n = 50
x1 = rep((1:n) / n, each = n)
x2 = rep((1:n) / n, n)
X = cbind(x1, x2)
e = 0.5

m = (1/25 * x1 + 0.75 * x2) < 0.5
tau0 = x1 < 0.6 & x2 < 0.6
sigma = 1

mse_r = rep(0, N)
mse_s = rep(0, N)
mse_t = rep(0, N)
mse_x = rep(0, N)
mse_f = rep(0, N)

for (i in 1:N){
  z = rbinom(n^2, 1, e)
  
  y = m + (z - 0.5) * tau0 + sigma * rnorm(n^2)
    
  rlasso_fit = rlasso(X, z, y)
  rlasso_est = predict(rlasso_fit, X)
  mse_r[i] = mean((rlasso_est - tau0)^2)
  
  slasso_fit = slasso(X, z, y)
  slasso_est = predict(slasso_fit, X)
  mse_s[i] = mean((slasso_est - tau0)^2)
  
  tlasso_fit = tlasso(X, z, y)
  tlasso_est = predict(tlasso_fit, X)
  mse_t[i] = mean((tlasso_est - tau0)^2)
  
  xlasso_fit = xlasso(X, z, y)
  xlasso_est = predict(xlasso_fit, X)
  mse_x[i] = mean((xlasso_est - tau0)^2)
  
  
  f_fit = causal_forest(X, y, z, num.trees = 2000, sample.fraction = 0.1)
  pred = predict(f_fit)
  f_est  = pred$predictions
  mse_f[i] = mean((f_est - tau0)^2)
}

mean(mse_r)
mean(mse_s)
mean(mse_t)
mean(mse_x)
mean(mse_f)

# Scenario 2
set.seed(123)
N = 500
n = 50
x1 = rep((1:n) / n, each = n)
x2 = rep((1:n) / n, n)
X = cbind(x1, x2)

e = as.numeric(apply(X, 1, d25) < apply(X, 1, d75))
e[e == 0] = 0.1
e[e == 1] = 0.9

m = (1/25 * x1 + 0.75 * x2) < 0.5
tau0 = x1 < 0.6 & x2 < 0.6
sigma = 1

mse_r = rep(0, N)
mse_s = rep(0, N)
mse_t = rep(0, N)
mse_x = rep(0, N)
mse_f = rep(0, N)

for (i in 1:N){
  z = rbinom(n^2, 1, e)
  
  y = m + (z - 0.5) * tau0 + sigma * rnorm(n^2)
    
  rlasso_fit = rlasso(X, z, y)
  rlasso_est = predict(rlasso_fit, X)
  mse_r[i] = mean((rlasso_est - tau0)^2)
  
  slasso_fit = slasso(X, z, y)
  slasso_est = predict(slasso_fit, X)
  mse_s[i] = mean((slasso_est - tau0)^2)
  
  tlasso_fit = tlasso(X, z, y)
  tlasso_est = predict(tlasso_fit, X)
  mse_t[i] = mean((tlasso_est - tau0)^2)
  
  xlasso_fit = xlasso(X, z, y)
  xlasso_est = predict(xlasso_fit, X)
  mse_x[i] = mean((xlasso_est - tau0)^2)
  
  f_fit = causal_forest(X, y, z, num.trees = 2000, sample.fraction = 0.1)
  pred = predict(f_fit)
  f_est  = pred$predictions
  mse_f[i] = mean((f_est - tau0)^2)
}

mean(mse_r)
mean(mse_s)
mean(mse_t)
mean(mse_x)
mean(mse_f)

# Scenario 3
set.seed(123)

N = 500
n = 10000
d = 5

e = 0.5

sigma = 0.75

mse_r = rep(0, N)
mse_s = rep(0, N)
mse_t = rep(0, N)
mse_x = rep(0, N)
mse_f = rep(0, N)

for (i in 1:N){
  X = matrix(runif(d*n), ncol = d)
  m = sin(pi*X[,1]*X[,2]) + 2*(X[,3]-0.5)^2 + X[,4] + 0.5 * X[,5]

  tau0 = apply(X, 1, d_tau)
  
  z = rbinom(n, 1, e)
  
  y = m + (z - 0.5) * tau0 + sigma * rnorm(n)
  
  rlasso_fit = rlasso(X, z, y)
  rlasso_est = predict(rlasso_fit, X)
  mse_r[i] = mean((rlasso_est - tau0)^2)
  
  slasso_fit = slasso(X, z, y)
  slasso_est = predict(slasso_fit, X)
  mse_s[i] = mean((slasso_est - tau0)^2)
  
  tlasso_fit = tlasso(X, z, y)
  tlasso_est = predict(tlasso_fit, X)
  mse_t[i] = mean((tlasso_est - tau0)^2)
  
  xlasso_fit = xlasso(X, z, y)
  xlasso_est = predict(xlasso_fit, X)
  mse_x[i] = mean((xlasso_est - tau0)^2)
  
  f_fit = causal_forest(X, y, z, num.trees = 2000, sample.fraction = 0.1)
  pred = predict(f_fit)
  f_est  = pred$predictions
  mse_f[i] = mean((f_est - tau0)^2)
}

mean(mse_r)
mean(mse_s)
mean(mse_t)
mean(mse_x)
mean(mse_f)

# Scenario 4
set.seed(123)

N = 500
d = 10
n = 10000

sigma = 0.5

mse_r = rep(0, N)
mse_s = rep(0, N)
mse_t = rep(0, N)
mse_x = rep(0, N)
mse_f = rep(0, N)

for (i in 1:N){
  X = matrix(runif(d*n), ncol = d)
  beta = runif(d)
  
  m = 2 * as.numeric(drop(X %*% beta) > round(d/2) - 1)
  m[m == 0] = -2

  e = as.numeric(drop(X^2 %*% beta) > floor(sqrt(d)) - 1)
  e[e == 0] = 0.25
  e[e == 1] = 0.75
  
  tau0 = 1
  
  z = rbinom(n, 1, e)
  
  y = m + (z - 0.5) * tau0 + sigma * rnorm(n)
  
  rlasso_fit = rlasso(X, z, y)
  rlasso_est = predict(rlasso_fit, X)
  mse_r[i] = mean((rlasso_est - tau0)^2)
  
  slasso_fit = slasso(X, z, y)
  slasso_est = predict(slasso_fit, X)
  mse_s[i] = mean((slasso_est - tau0)^2)
  
  tlasso_fit = tlasso(X, z, y)
  tlasso_est = predict(tlasso_fit, X)
  mse_t[i] = mean((tlasso_est - tau0)^2)
  
  xlasso_fit = xlasso(X, z, y)
  xlasso_est = predict(xlasso_fit, X)
  mse_x[i] = mean((xlasso_est - tau0)^2)
  
  ulasso_fit = ulasso(X, z, y)
  ulasso_est = predict(ulasso_fit, X)
  mse_u[i] = mean((ulasso_est - tau0)^2)
  
  f_fit = causal_forest(X, y, z, num.trees = 2000, sample.fraction = 0.1)
  pred = predict(f_fit)
  f_est  = pred$predictions
  mse_f[i] = mean((f_est - tau0)^2)
}

mean(mse_r)
mean(mse_s)
mean(mse_t)
mean(mse_x)
mean(mse_f)

# Scenario 5
set.seed(123)

N = 500
d = 10
n = 10000

sigma = 0.5

mse_r = rep(0, N)
mse_s = rep(0, N)
mse_t = rep(0, N)
mse_x = rep(0, N)
mse_f = rep(0, N)

for (i in 1:N){
  X = matrix(runif(d*n), ncol = d)
  beta = runif(d/2)
  
  xt = t(matrix(c(X[,1] * X[,2], X[,3] * X[,4], X[,5] * X[,6], X[,7] * X[,8], X[,9] * X[,10]), nrow = 5))
  
  m = drop(xt %*% beta)

  e = apply(X, 1, d05)
  
  tau0 = apply(X, 1, d10)
  
  z = rbinom(n, 1, e)
  
  y = m + (z) * tau0 + sigma * rnorm(n)
  
  rlasso_fit = rlasso(X, z, y)
  rlasso_est = predict(rlasso_fit, X)
  mse_r[i] = mean((rlasso_est - tau0)^2)
  
  slasso_fit = slasso(X, z, y)
  slasso_est = predict(slasso_fit, X)
  mse_s[i] = mean((slasso_est - tau0)^2)
  
  tlasso_fit = tlasso(X, z, y)
  tlasso_est = predict(tlasso_fit, X)
  mse_t[i] = mean((tlasso_est - tau0)^2)
  
  xlasso_fit = xlasso(X, z, y)
  xlasso_est = predict(xlasso_fit, X)
  mse_x[i] = mean((xlasso_est - tau0)^2)
  
  ulasso_fit = ulasso(X, z, y)
  ulasso_est = predict(ulasso_fit, X)
  mse_u[i] = mean((ulasso_est - tau0)^2)
  
  f_fit = causal_forest(X, y, z)
  pred = predict(f_fit)
  f_est  = pred$predictions
  mse_f[i] = mean((f_est - tau0)^2)
}

mean(mse_r)
mean(mse_s)
mean(mse_t)
mean(mse_x)
mean(mse_f)
