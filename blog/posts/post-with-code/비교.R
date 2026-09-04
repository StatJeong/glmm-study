####################################################
## Logistic random-intercept GLMM
## Unknown sigma_b^2 (Inv-Gamma prior)
## Stan (NUTS) vs PG-Gibbs vs MH (fixed iter, progress)
## Slightly more complex setting than before
####################################################

library(rstan)
library(BayesLogit)
library(MASS)
library(coda)
library(ggplot2)
library(bayesplot)
library(dplyr)
library(tidyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(2025)

####################################################
## 1. Simulate data (조금 더 복잡한 셋팅)
####################################################

J      <- 120   # groups (80 -> 120)
n_per  <- 4     # obs per group (5 -> 4)
N      <- J * n_per

beta0_true   <- 0
beta1_true   <- 1
sigma_b_true <- 1.5

group <- rep(1:J, each = n_per)
b_j   <- rnorm(J, 0, sigma_b_true)
b_i   <- b_j[group]

x   <- rnorm(N)
eta <- beta0_true + beta1_true * x + b_i
p   <- 1 / (1 + exp(-eta))
y   <- rbinom(N, 1, p)

dat <- data.frame(y = y, x = x, group = group)

## prior for sigma_b^2: Inv-Gamma(a0, b0)
a0 <- 2
b0 <- 2

####################################################
## 2. Stan model (unknown sigma_b^2)
####################################################

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1,upper=J> group[N];
  vector[N] x;
  int<lower=0,upper=1> y[N];
  real<lower=0> a0;
  real<lower=0> b0;
}
parameters {
  real beta0;
  real beta1;
  real<lower=0> sigma2_b;
  vector[J] b;
}
transformed parameters {
  real<lower=0> sigma_b;
  sigma_b = sqrt(sigma2_b);
}
model {
  beta0   ~ normal(0, 10);
  beta1   ~ normal(0, 10);
  sigma2_b ~ inv_gamma(a0, b0);
  b       ~ normal(0, sigma_b);
  y       ~ bernoulli_logit(beta0 + beta1 * x + b[group]);
}
"

stan_dat <- list(
  N     = N,
  J     = J,
  group = group,
  x     = x,
  y     = y,
  a0    = a0,
  b0    = b0
)

sm <- stan_model(model_code = stan_code)

time_stan <- system.time({
  fit_stan <- sampling(
    sm,
    data   = stan_dat,
    iter   = 3000,
    warmup = 500,
    chains = 1,
    refresh = 0
  )
})

post_stan <- as.matrix(fit_stan,
                       pars = c("beta0","beta1","sigma_b", paste0("b[",1:J,"]")))
n_stan    <- nrow(post_stan)

####################################################
## 3. PG-Gibbs sampler (beta,b | sigma2_b via PG,
##    sigma2_b | b via Inv-Gamma)
####################################################

X_beta <- cbind(1, x)
Z      <- matrix(0, nrow = N, ncol = J)
Z[cbind(1:N, group)] <- 1
X_full <- cbind(X_beta, Z)        # N x (2+J)
K      <- ncol(X_full)

prior_var_beta <- 10^2

n_iter_pg <- 3000
burn_pg   <- 500

theta_pg   <- matrix(NA, nrow = n_iter_pg, ncol = K)
sigma2_pg  <- numeric(n_iter_pg)
colnames(theta_pg) <- c("beta0","beta1",paste0("b[",1:J,"]"))

theta_curr  <- c(0,0,rep(0,J))
sigma2_curr <- 1.0   # init

time_pg <- system.time({
  for (it in 1:n_iter_pg) {
    eta_pg <- as.numeric(X_full %*% theta_curr)
    omega  <- BayesLogit::rpg(N, h = 1, z = eta_pg)
    WX     <- X_full * omega
    XtOmegaX <- crossprod(X_full, WX)
    
    prior_var_vec <- c(rep(prior_var_beta,2), rep(sigma2_curr,J))
    V0_inv        <- diag(1/prior_var_vec)
    
    Sigma_post <- solve(XtOmegaX + V0_inv)
    mu_post    <- Sigma_post %*% crossprod(X_full, y - 0.5)
    
    theta_curr <- as.numeric(MASS::mvrnorm(1, mu = mu_post, Sigma = Sigma_post))
    
    b_curr <- theta_curr[3:(2+J)]
    shape_post <- a0 + J/2
    rate_post  <- b0 + sum(b_curr^2)/2
    sigma2_curr <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
    
    theta_pg[it, ]  <- theta_curr
    sigma2_pg[it]   <- sigma2_curr
  }
})

theta_pg_keep  <- theta_pg[(burn_pg+1):n_iter_pg,,drop=FALSE]
sigma2_pg_keep <- sigma2_pg[(burn_pg+1):n_iter_pg]
post_pg <- cbind(
  beta0   = theta_pg_keep[,"beta0"],
  beta1   = theta_pg_keep[,"beta1"],
  sigma_b = sqrt(sigma2_pg_keep),
  theta_pg_keep[,paste0("b[",1:J,"]"),drop=FALSE]
)
n_pg <- nrow(post_pg)

####################################################
## 4. MH sampler (fixed iter = 3000, show acc rate)
##    - sigma2_b는 Gibbs 업데이트
####################################################

log_post_beta_b <- function(theta, sigma2_b, y, x, group,
                            prior_sd_beta = 10) {
  beta0 <- theta[1]
  beta1 <- theta[2]
  b     <- theta[3:length(theta)]
  
  eta <- beta0 + beta1 * x + b[group]
  ll  <- sum(y * eta - log1p(exp(eta)))
  
  lp_beta0 <- dnorm(beta0, 0, prior_sd_beta, log=TRUE)
  lp_beta1 <- dnorm(beta1, 0, prior_sd_beta, log=TRUE)
  lp_b     <- sum(dnorm(b, 0, sqrt(sigma2_b), log=TRUE))
  
  return(ll + lp_beta0 + lp_beta1 + lp_b)
}

n_iter_mh <- 3000
burn_mh   <- 500
theta_dim <- 2 + J

theta_mh  <- matrix(NA, nrow = n_iter_mh, ncol = theta_dim)
sigma2_mh <- numeric(n_iter_mh)
colnames(theta_mh) <- c("beta0","beta1",paste0("b[",1:J,"]"))

theta_curr  <- c(0,0,rep(0,J))
sigma2_curr <- 1.0
lp_curr     <- log_post_beta_b(theta_curr, sigma2_curr, y, x, group)

## 제안분포 표준편차 (조금 키움)
prop_sd_beta <- 0.25
prop_sd_b    <- 0.40
prop_sd_vec  <- c(prop_sd_beta, prop_sd_beta, rep(prop_sd_b,J))

acc_cnt <- 0

time_mh <- system.time({
  for (it in 1:n_iter_mh) {
    
    theta_prop <- theta_curr + rnorm(theta_dim, 0, prop_sd_vec)
    lp_prop    <- log_post_beta_b(theta_prop, sigma2_curr, y, x, group)
    
    log_alpha <- lp_prop - lp_curr
    if (log(runif(1)) < log_alpha) {
      theta_curr <- theta_prop
      lp_curr    <- lp_prop
      acc_cnt    <- acc_cnt + 1
    }
    
    ## sigma2_b Gibbs 업데이트
    b_curr     <- theta_curr[3:(2+J)]
    shape_post <- a0 + J/2
    rate_post  <- b0 + sum(b_curr^2)/2
    sigma2_curr <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
    
    theta_mh[it, ]  <- theta_curr
    sigma2_mh[it]   <- sigma2_curr
  }
})

accept_rate_mh <- acc_cnt / n_iter_mh
cat("MH acceptance rate:", accept_rate_mh, "\n")

theta_mh_keep  <- theta_mh[(burn_mh+1):n_iter_mh,,drop=FALSE]
sigma2_mh_keep <- sigma2_mh[(burn_mh+1):n_iter_mh]

post_mh <- cbind(
  beta0   = theta_mh_keep[,"beta0"],
  beta1   = theta_mh_keep[,"beta1"],
  sigma_b = sqrt(sigma2_mh_keep),
  theta_mh_keep[,paste0("b[",1:J,"]"),drop=FALSE]
)
n_mh <- nrow(post_mh)


####################################################
## 5. ESS, ESR1(ESS/n), ESR2(ESS/time) + runtime
####################################################

mcmc_stan <- mcmc(post_stan[,c("beta0","beta1","sigma_b")])
mcmc_pg   <- mcmc(post_pg[,c("beta0","beta1","sigma_b")])
mcmc_mh   <- mcmc(post_mh[,c("beta0","beta1","sigma_b")])

ess_stan <- effectiveSize(mcmc_stan)
ess_pg   <- effectiveSize(mcmc_pg)
ess_mh   <- effectiveSize(mcmc_mh)

esr1_stan <- ess_stan / n_stan
esr1_pg   <- ess_pg   / n_pg
esr1_mh   <- ess_mh   / n_mh

esr2_stan <- ess_stan / time_stan["elapsed"]
esr2_pg   <- ess_pg   / time_pg["elapsed"]
esr2_mh   <- ess_mh   / time_mh["elapsed"]

cat("ESS (beta0, beta1, sigma_b)\n")
print(rbind(Stan = ess_stan, PG = ess_pg, MH = ess_mh))

cat("ESR1 = ESS / #samples\n")
print(rbind(Stan = esr1_stan, PG = esr1_pg, MH = esr1_mh))

cat("ESR2 = ESS / elapsed_time\n")
print(rbind(Stan = esr2_stan, PG = esr2_pg, MH = esr2_mh))

cat("Run times (seconds)\n")
runtime <- data.frame(
  method  = c("Stan","PG Gibbs","MH"),
  elapsed = c(time_stan["elapsed"],
              time_pg["elapsed"],
              time_mh["elapsed"])
)
print(runtime)

####################################################
## 6. True vs posterior (mean & 95% CI) summary
####################################################

true_par <- c(beta0 = beta0_true,
              beta1 = beta1_true,
              sigma_b = sigma_b_true)

summ_one <- function(mat, method) {
  pars <- c("beta0","beta1","sigma_b")
  m    <- apply(mat[,pars,drop=FALSE], 2, mean)
  l95  <- apply(mat[,pars,drop=FALSE], 2, quantile, probs = 0.025)
  u95  <- apply(mat[,pars,drop=FALSE], 2, quantile, probs = 0.975)
  
  data.frame(
    method = method,
    param  = pars,
    true   = true_par[pars],
    mean   = m,
    low95  = l95,
    high95 = u95,
    row.names = NULL
  )
}

summ_stan_par <- summ_one(post_stan, "Stan (NUTS)")
summ_pg_par   <- summ_one(post_pg,   "PG Gibbs")
summ_mh_par   <- summ_one(post_mh,   "MH")

summ_all_par <- rbind(summ_stan_par, summ_pg_par, summ_mh_par)
print(summ_all_par)

## Plot for true vs posterior intervals
ggplot(summ_all_par,
       aes(x = param, y = mean,
           ymin = low95, ymax = high95,
           colour = method)) +
  geom_pointrange(position = position_dodge(width = 0.4)) +
  geom_hline(aes(yintercept = true),
             linetype = 2, colour = "black") +
  facet_wrap(~ method, nrow = 1) +
  labs(x = "Parameter", y = "Posterior mean & 95% CI") +
  theme_bw()

####################################################
## 7. Trace plots (beta0, beta1, sigma_b)
####################################################

bayesplot::mcmc_trace(as.matrix(post_stan[,c("beta0","beta1","sigma_b")]),
                      pars = c("beta0","beta1","sigma_b"),
                      facet_args = list(ncol = 1)) +
  ggtitle("Stan (NUTS) trace")

bayesplot::mcmc_trace(as.matrix(post_pg[,c("beta0","beta1","sigma_b")]),
                      pars = c("beta0","beta1","sigma_b"),
                      facet_args = list(ncol = 1)) +
  ggtitle("PG Gibbs trace")

bayesplot::mcmc_trace(as.matrix(post_mh[,c("beta0","beta1","sigma_b")]),
                      pars = c("beta0","beta1","sigma_b"),
                      facet_args = list(ncol = 1)) +
  ggtitle("MH trace")

####################################################
## 8. Random-effect plots (true vs posterior)
####################################################

samples_stan_b <- post_stan[,paste0("b[",1:J,"]")]
samples_pg_b   <- post_pg[,paste0("b[",1:J,"]")]
samples_mh_b   <- post_mh[,paste0("b[",1:J,"]")]

summ_re <- function(samples_b_mat, method_name, true_b) {
  J <- ncol(samples_b_mat)
  means  <- apply(samples_b_mat, 2, mean)
  q05    <- apply(samples_b_mat, 2, quantile, probs = 0.05)
  q25    <- apply(samples_b_mat, 2, quantile, probs = 0.25)
  q75    <- apply(samples_b_mat, 2, quantile, probs = 0.75)
  q95    <- apply(samples_b_mat, 2, quantile, probs = 0.95)
  
  data.frame(
    method  = method_name,
    group_id = 1:J,
    group_label = paste0("G",1:J),
    mean   = means,
    low90  = q05,
    high90 = q95,
    low50  = q25,
    high50 = q75,
    true   = true_b
  )
}

summ_stan <- summ_re(samples_stan_b, "Stan (NUTS)", b_j)
summ_pg   <- summ_re(samples_pg_b,   "PG Gibbs",   b_j)
summ_mh   <- summ_re(samples_mh_b,   "MH",         b_j)

summ_all <- rbind(summ_mh, summ_pg, summ_stan)

ord <- order(b_j)
group_levels <- paste0("G",1:J)[ord]
summ_all$group_factor <- factor(summ_all$group_label, levels = group_levels)

ggplot(summ_all, aes(x = group_factor)) +
  geom_linerange(aes(ymin = low90, ymax = high90),
                 colour = "grey80", size = 0.7) +
  geom_linerange(aes(ymin = low50, ymax = high50),
                 colour = "grey40", size = 1.4) +
  geom_point(aes(y = mean),
             colour = "black", size = 1.3) +
  geom_point(aes(y = true),
             colour = "red", size = 1.3, shape = 16) +
  facet_wrap(~ method, ncol = 1) +
  labs(x = "Group", y = "Random Effect") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor = element_blank()
  )


########################################

####################################################
## 0. 패키지 & 실제 데이터 불러오기
####################################################

library(HSAUR2)      # toenail 데이터
data("toenail")
head(toenail)
library(dplyr)
library(rstan)
library(BayesLogit)
library(MASS)
library(coda)
library(ggplot2)
library(bayesplot)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

####################################################
## 1. 데이터 전처리 + y, X, Z 정의
####################################################

toenail_bin <- toenail %>%
  mutate(
    y = ifelse(outcome == "moderate or severe", 1L, 0L),
    patientID = factor(patientID),
    visit = as.numeric(visit),
    treatment = factor(treatment)
  ) %>%
  filter(!is.na(y))

y  <- toenail_bin$y
X  <- model.matrix(~ treatment * visit, toenail_bin)      # 고정효과 디자인행렬
Z  <- model.matrix(~ 0 + patientID, toenail_bin)          # 랜덤효과 (환자별 절편)

n <- length(y)
p <- ncol(X)
J <- ncol(Z)

patient_index <- as.integer(toenail_bin$patientID)        # Stan용 인덱스 (1..J)

####################################################
## 공통 prior 파라미터
####################################################

a0 <- 2
b0 <- 2
prior_var_beta <- 10^2

####################################################
## 2. Stan (NUTS)
####################################################

stan_code_toenail <- "
data {
  int<lower=1> N;
  int<lower=1> p;
  int<lower=1> J;
  int<lower=1,upper=J> patient[N];
  matrix[N,p] X;
  int<lower=0,upper=1> y[N];
  real<lower=0> a0;
  real<lower=0> b0;
}
parameters {
  vector[p] beta;
  real<lower=0> sigma2_b;
  vector[J] b;
}
transformed parameters {
  real<lower=0> sigma_b;
  sigma_b = sqrt(sigma2_b);
}
model {
  beta    ~ normal(0, 10);
  sigma2_b ~ inv_gamma(a0, b0);
  b       ~ normal(0, sigma_b);
  y       ~ bernoulli_logit(X * beta + b[patient]);
}
"

stan_dat_toenail <- list(
  N       = n,
  p       = p,
  J       = J,
  patient = patient_index,
  X       = X,
  y       = y,
  a0      = a0,
  b0      = b0
)

sm_toenail <- stan_model(model_code = stan_code_toenail)

time_stan_toe <- system.time({
  fit_stan_toe <- sampling(
    sm_toenail,
    data   = stan_dat_toenail,
    iter   = 3000,
    warmup = 500,
    chains = 1,
    refresh = 0
  )
})

post_stan_toe <- as.matrix(
  fit_stan_toe,
  pars = c("beta", "sigma_b", paste0("b[", 1:J, "]"))
)
n_stan_toe <- nrow(post_stan_toe)

####################################################
## 3. PG-Gibbs (beta, b | sigma2_b via PG, 
##              sigma2_b | b via Inv-Gamma)
####################################################

X_full <- cbind(X, Z)       # N x (p+J)
K      <- ncol(X_full)

n_iter_pg <- 3000
burn_pg   <- 500

theta_pg_toe  <- matrix(NA, nrow = n_iter_pg, ncol = K)
sigma2_pg_toe <- numeric(n_iter_pg)

colnames(theta_pg_toe) <- c(colnames(X), colnames(Z))

theta_curr  <- rep(0, K)   # (beta, b) 초기값
sigma2_curr <- 1.0

time_pg_toe <- system.time({
  for (it in 1:n_iter_pg) {
    eta_pg <- as.numeric(X_full %*% theta_curr)
    omega  <- BayesLogit::rpg(n, h = 1, z = eta_pg)
    WX     <- X_full * omega
    XtOmegaX <- crossprod(X_full, WX)
    
    prior_var_vec <- c(rep(prior_var_beta, p), rep(sigma2_curr, J))
    V0_inv        <- diag(1 / prior_var_vec)
    
    Sigma_post <- solve(XtOmegaX + V0_inv)
    mu_post    <- Sigma_post %*% crossprod(X_full, y - 0.5)
    
    theta_curr <- as.numeric(MASS::mvrnorm(1, mu = mu_post, Sigma = Sigma_post))
    
    b_curr     <- theta_curr[(p + 1):(p + J)]
    shape_post <- a0 + J / 2
    rate_post  <- b0 + sum(b_curr^2) / 2
    sigma2_curr <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
    
    theta_pg_toe[it, ]  <- theta_curr
    sigma2_pg_toe[it]   <- sigma2_curr
  }
})

theta_pg_keep_toe  <- theta_pg_toe[(burn_pg + 1):n_iter_pg, , drop = FALSE]
sigma2_pg_keep_toe <- sigma2_pg_toe[(burn_pg + 1):n_iter_pg]

post_pg_toe <- cbind(
  beta    = theta_pg_keep_toe[, 1:p, drop = FALSE],
  sigma_b = sqrt(sigma2_pg_keep_toe),
  theta_pg_keep_toe[, (p + 1):(p + J), drop = FALSE]
)
colnames(post_pg_toe) <- c(
  paste0("beta[", 1:p, "]"),
  "sigma_b",
  paste0("b[", 1:J, "]")
)
n_pg_toe <- nrow(post_pg_toe)

####################################################
## 4. MH (beta, b random-walk + sigma2_b Gibbs)
####################################################

log_post_toe <- function(theta, sigma2_b, y, X, Z,
                         prior_sd_beta = 10) {
  p <- ncol(X)
  J <- ncol(Z)
  
  beta <- theta[1:p]
  b    <- theta[(p + 1):(p + J)]
  
  eta <- as.numeric(X %*% beta + Z %*% b)
  ll  <- sum(y * eta - log1p(exp(eta)))
  
  lp_beta <- sum(dnorm(beta, 0, prior_sd_beta, log = TRUE))
  lp_b    <- sum(dnorm(b, 0, sqrt(sigma2_b), log = TRUE))
  
  ll + lp_beta + lp_b
}

n_iter_mh <- 3000
burn_mh   <- 500
theta_dim <- p + J

theta_mh_toe  <- matrix(NA, nrow = n_iter_mh, ncol = theta_dim)
sigma2_mh_toe <- numeric(n_iter_mh)

colnames(theta_mh_toe) <- c(paste0("beta[",1:p,"]"), paste0("b[",1:J,"]"))

theta_curr  <- rep(0, theta_dim)
sigma2_curr <- 1.0
lp_curr     <- log_post_toe(theta_curr, sigma2_curr, y, X, Z)

prop_sd_beta <- 0.25
prop_sd_b    <- 0.40
prop_sd_vec  <- c(rep(prop_sd_beta, p), rep(prop_sd_b, J))

acc_cnt <- 0

time_mh_toe <- system.time({
  for (it in 1:n_iter_mh) {
    theta_prop <- theta_curr + rnorm(theta_dim, 0, prop_sd_vec)
    lp_prop    <- log_post_toe(theta_prop, sigma2_curr, y, X, Z)
    
    log_alpha <- lp_prop - lp_curr
    if (log(runif(1)) < log_alpha) {
      theta_curr <- theta_prop
      lp_curr    <- lp_prop
      acc_cnt    <- acc_cnt + 1
    }
    
    b_curr     <- theta_curr[(p + 1):(p + J)]
    shape_post <- a0 + J / 2
    rate_post  <- b0 + sum(b_curr^2) / 2
    sigma2_curr <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
    
    theta_mh_toe[it, ]  <- theta_curr
    sigma2_mh_toe[it]   <- sigma2_curr
  }
})

accept_rate_mh_toe <- acc_cnt / n_iter_mh
cat("MH (toenail) acceptance rate:", accept_rate_mh_toe, "\n")

theta_mh_keep_toe  <- theta_mh_toe[(burn_mh + 1):n_iter_mh, , drop = FALSE]
sigma2_mh_keep_toe <- sigma2_mh_toe[(burn_mh + 1):n_iter_mh]

post_mh_toe <- cbind(
  beta    = theta_mh_keep_toe[, 1:p, drop = FALSE],
  sigma_b = sqrt(sigma2_mh_keep_toe),
  theta_mh_keep_toe[, (p + 1):(p + J), drop = FALSE]
)
colnames(post_mh_toe) <- c(
  paste0("beta[", 1:p, "]"),
  "sigma_b",
  paste0("b[", 1:J, "]")
)
n_mh_toe <- nrow(post_mh_toe)


