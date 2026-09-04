library(MASS)     # mvrnorm
library(BayesLogit)  # PG sampler
library(microbenchmark)  # for benchmarking
library(mcmcse)   # for effectiveSize (ESS) 계산
# install.packages(c("BayesLogit","mcmcse","microbenchmark"))
library(HSAUR2)
library(dplyr)
#install.packages("coda")
library(coda)

set.seed(1015)

data("toenail")

toenail_bin <- toenail %>%
  mutate(
    y = ifelse(outcome == "moderate or severe", 1, 0),
    patientID = factor(patientID),
    visit = as.numeric(visit),
    treatment = factor(treatment)
  )

y  <- toenail_bin$y
X  <- model.matrix(~ treatment * visit, toenail_bin)
Z  <- model.matrix(~ 0 + patientID, toenail_bin)

n <- length(y)
p <- ncol(X)
J <- ncol(Z)

####################################################
log_post_glmm <- function(theta, X, Z, y, p, J, sigma2_beta, a_tau, b_tau) {
  
  beta <- theta[1:p]
  u    <- theta[(p+1):(p+J)]
  tau2 <- theta[p+J+1]
  
  # logit likelihood
  eta <- X %*% beta + Z %*% u
  ll  <- sum(y * eta - log1p(exp(eta)))
  
  # priors
  lp_beta <- - sum(beta^2)/(2*sigma2_beta)
  lp_u    <- - sum(u^2)/(2*tau2)
  lp_tau  <- (a_tau-1)*log(tau2) - b_tau*tau2
  
  return(ll + lp_beta + lp_u + lp_tau)
}

mh_glmm_acceptN <- function(target_accept = 5000,
                            X, Z, y, p, J,
                            sigma2_beta = 100,
                            a_tau = 1, b_tau = 1,
                            proposal_SD = 0.02) {
  
  # 초기값
  beta <- rep(0, p)
  u    <- rep(0, J)
  tau2 <- 1
  theta <- c(beta, u, tau2)
  
  accepted <- list()
  n_acc <- 0
  iter <- 0
  
  while (n_acc < target_accept) {
    
    iter <- iter + 1
    
    theta_prop <- as.numeric(theta + proposal_SD * rnorm(p+J+1))
    
    log_r <- log_post_glmm(theta_prop, X, Z, y, p, J,
                           sigma2_beta, a_tau, b_tau) -
      log_post_glmm(theta,      X, Z, y, p, J,
                    sigma2_beta, a_tau, b_tau)
    
    if (log(runif(1)) < log_r) {
      theta <- theta_prop
      n_acc <- n_acc + 1
      accepted[[n_acc]] <- theta
    }
  }
  
  samples <- do.call(rbind, accepted)
  list(samples = samples,
       total_iterations = iter,
       acceptance_rate = target_accept/iter)
}
######################################################
n_iter <- 5000
burn   <- 1000

### PG-Gibbs ###
time_pg <- system.time({
  pg_out <- pg_glmm(n_iter, X, Z, y)
})
post_pg <- pg_out$beta[(burn+1):n_iter, ]  # β만 비교 예시

### MH (while) ###
time_mh <- system.time({
  mh_out <- mh_glmm_acceptN(
    target_accept = n_iter,
    X=X, Z=Z, y=y, p=p, J=J
  )
})

samp_mh <- mh_out$samples
mh_total_iter <- mh_out$total_iterations
mh_accept <- mh_out$acceptance_rate

post_mh <- samp_mh[(burn+1):n_iter, 1:p]  # β 부분만 ESS 비교

### ESS
ess_pg <- apply(post_pg,  2, coda::effectiveSize)
ess_mh <- apply(post_mh,  2, coda::effectiveSize)

### ESR (two types)
ESR_pg_sec  <- ess_pg / time_pg[3]
ESR_pg_iter <- ess_pg / n_iter

ESR_mh_sec  <- ess_mh / time_mh[3]
ESR_mh_iter <- ess_mh / mh_total_iter

print(ess_pg)
print(ess_mh)
print(ESR_pg_sec)
print(ESR_mh_sec)
print(ESR_pg_iter)
print(ESR_mh_iter)



#################################################################
compare_GLMM_PG_MH <- function(y, X, Z, 
                               n_iter = 5000, burn = 1000,
                               sigma2_beta = 100,
                               a_tau = 1, b_tau = 1,
                               proposal_SD = 0.02) {
  
  n <- length(y)
  p <- ncol(X)
  J <- ncol(Z)
  
  ## ---------------------------------------------------------
  ## 1. PG-Gibbs GLMM
  ## ---------------------------------------------------------
  pg_glmm <- function(n_iter, X, Z, y, sigma2_beta, a_tau, b_tau) {
    
    n <- length(y); p <- ncol(X); J <- ncol(Z)
    
    beta <- rep(0, p)
    u    <- rep(0, J)
    tau2 <- 1
    
    out_beta <- matrix(NA, n_iter, p)
    out_u    <- matrix(NA, n_iter, J)
    out_tau2 <- numeric(n_iter)
    
    for (iter in 1:n_iter) {
      eta   <- X %*% beta + Z %*% u
      omega <- rpg(n, 1, eta)
      
      W        <- diag(omega)
      y_tilde  <- y - 0.5
      Xjoint   <- cbind(X, Z)
      Prec_prior <- diag(c(rep(1/sigma2_beta, p), rep(1/tau2, J)))
      
      Prec_post <- t(Xjoint) %*% W %*% Xjoint + Prec_prior
      V_post    <- solve(Prec_post)
      m_post    <- V_post %*% (t(Xjoint) %*% y_tilde)
      
      theta <- as.numeric(m_post + chol(V_post) %*% rnorm(p + J))
      
      beta <- theta[1:p]
      u    <- theta[(p+1):(p+J)]
      
      ## update tau2
      shape <- a_tau + J/2
      rate  <- b_tau + sum(u^2)/2
      tau2  <- 1 / rgamma(1, shape = shape, rate = rate)
      
      out_beta[iter, ] <- beta
      out_u[iter, ]    <- u
      out_tau2[iter]   <- tau2
      
      ## ★ 진행상황 출력 (예: 500번마다)
      if (iter %% 500 == 0) {
        cat("[PG] iter:", iter, "/", n_iter, "\n")
        flush.console()
      }
    }
    
    ## ★ 루프 끝난 다음에 리스트 반환
    list(beta = out_beta,
         u    = out_u,
         tau2 = out_tau2)
  }
  
  ## PG run
  time_pg <- system.time({
    pg_out <- pg_glmm(n_iter, X, Z, y, sigma2_beta, a_tau, b_tau)
  })
  post_pg <- pg_out$beta[(burn + 1):n_iter, ]
  
  
  ## ---------------------------------------------------------
  ## 2. MH (while) GLMM
  ## ---------------------------------------------------------
  log_post_glmm <- function(theta) {
    beta <- theta[1:p]
    u    <- theta[(p+1):(p+J)]
    tau2 <- theta[p+J+1]
    
    eta <- X %*% beta + Z %*% u
    ll  <- sum(y * eta - log1p(exp(eta)))
    
    lp_b <- - sum(beta^2) / (2 * sigma2_beta)
    lp_u <- - sum(u^2)    / (2 * tau2)
    lp_t <- (a_tau - 1) * log(tau2) - b_tau * tau2
    
    ll + lp_b + lp_u + lp_t
  }
  
  mh_sampler_acceptN <- function(target_accept) {
    theta <- c(rep(0, p), rep(0, J), 1)  # beta, u, tau2
    n_acc <- 0
    iter  <- 0
    accepted <- list()
    
    while (n_acc < target_accept) {
      iter <- iter + 1
      theta_prop <- theta + proposal_SD * rnorm(p + J + 1)
      
      log_r <- log_post_glmm(theta_prop) - log_post_glmm(theta)
      
      if (log(runif(1)) < log_r) {
        theta <- theta_prop
        n_acc <- n_acc + 1
        accepted[[n_acc]] <- theta
      }
      
      ## ★ 예: 10,000번 proposal마다 진행상황 출력
      if (iter %% 10000 == 0) {
        cat("[MH] iter:", iter, "accepted:", n_acc, "/", target_accept, "\n")
        flush.console()
      }
    }
    
    list(samples          = do.call(rbind, accepted),
         total_iterations = iter,
         acc_rate         = target_accept / iter)
  }
  
  ## MH run
  time_mh <- system.time({
    mh_out <- mh_sampler_acceptN(n_iter)
  })
  post_mh <- mh_out$samples[(burn + 1):n_iter, 1:p]  # beta만 ESS/ESR 비교
  
  
  ## ---------------------------------------------------------
  ## 3. Diagnostics: ESS / ESR
  ## ---------------------------------------------------------
  ess_pg <- coda::effectiveSize(as.mcmc(post_pg))
  ess_mh <- coda::effectiveSize(as.mcmc(post_mh))
  
  ESR_pg_sec <- ess_pg / time_pg[3]
  ESR_mh_sec <- ess_mh / time_mh[3]
  
  ESR_pg_iter <- ess_pg / n_iter
  ESR_mh_iter <- ess_mh / mh_out$total_iterations
  
  list(
    ESS      = list(PG = ess_pg, MH = ess_mh),
    ESR_sec  = list(PG = ESR_pg_sec, MH = ESR_mh_sec),
    ESR_iter = list(PG = ESR_pg_iter, MH = ESR_mh_iter),
    PG_time  = time_pg[3],
    MH_time  = time_mh[3],
    MH_total_iter = mh_out$total_iterations,
    PG_beta_draws = post_pg,
    MH_beta_draws = post_mh
  )
}

res <- compare_GLMM_PG_MH(y, X, Z)

res$ESS
res$ESR_sec
res$ESR_iter
res$PG_time
res$MH_time
res$MH_total_iter

plot(as.mcmc(res$PG_beta_draws))
plot(as.mcmc(res$MH_beta_draws))

acf(as.mcmc(res$PG_beta_draws))
acf(as.mcmc(res$MH_beta_draws))

##########################################################################

compare_GLMM_3methods_rstan <- function(y, X, Z, 
                                        n_iter = 5000, burn = 1000,
                                        iter_stan = 2000,
                                        sigma2_beta = 100,
                                        a_tau = 1, b_tau = 1,
                                        proposal_SD = 0.02) {
  
  library(BayesLogit)
  library(MASS)
  library(coda)
  library(brms)
  library(rstan)
  
  n <- length(y)
  p <- ncol(X)
  J <- ncol(Z)
  
  ###########################################################
  ## 1) PG-Gibbs GLMM
  ###########################################################
  pg_glmm <- function(n_iter, X, Z, y, sigma2_beta, a_tau, b_tau) {
    
    n <- length(y); p <- ncol(X); J <- ncol(Z)
    
    beta <- rep(0, p)
    u    <- rep(0, J)
    tau2 <- 1
    
    out_beta <- matrix(NA, n_iter, p)
    out_u    <- matrix(NA, n_iter, J)
    
    for (iter in 1:n_iter) {
      eta <- X %*% beta + Z %*% u
      omega <- rpg(n, 1, eta)
      
      W <- diag(omega)
      y_tilde <- y - 0.5
      Xjoint <- cbind(X, Z)
      
      Prec_prior <- diag(c(rep(1/sigma2_beta, p), rep(1/tau2, J)))
      Prec_post  <- t(Xjoint) %*% W %*% Xjoint + Prec_prior
      V_post     <- solve(Prec_post)
      m_post     <- V_post %*% (t(Xjoint) %*% y_tilde)
      
      theta <- as.numeric(m_post + chol(V_post) %*% rnorm(p+J))
      
      beta <- theta[1:p]
      u    <- theta[(p+1):(p+J)]
      
      shape <- a_tau + J/2
      rate  <- b_tau + sum(u^2)/2
      tau2  <- 1 / rgamma(1, shape=shape, rate=rate)
      
      out_beta[iter,] <- beta
      out_u[iter,]    <- u
      
      if (iter %% 500 == 0) {
        cat("[PG] iter:", iter, "/", n_iter, "\n")
        flush.console()
      }
    }
    
    list(beta = out_beta, u = out_u)
  }
  
  time_pg <- system.time({
    pg_out <- pg_glmm(n_iter, X, Z, y, sigma2_beta, a_tau, b_tau)
  })
  post_pg <- pg_out$beta[(burn+1):n_iter, ]
  
  
  ###########################################################
  ## 2) MH-while GLMM
  ###########################################################
  log_post_glmm <- function(theta) {
    beta <- theta[1:p]
    u    <- theta[(p+1):(p+J)]
    tau2 <- theta[p+J+1]
    
    eta <- X %*% beta + Z %*% u
    ll  <- sum(y * eta - log1p(exp(eta)))
    
    lp_b <- - sum(beta^2)/(2*sigma2_beta)
    lp_u <- - sum(u^2)/(2*tau2)
    lp_t <- (a_tau-1)*log(tau2) - b_tau*tau2
    
    ll + lp_b + lp_u + lp_t
  }
  
  mh_sampler_acceptN <- function(target_accept) {
    
    theta <- c(rep(0,p), rep(0,J), 1)
    n_acc <- 0
    iter  <- 0
    accepted <- list()
    
    while (n_acc < target_accept) {
      iter <- iter + 1
      
      theta_prop <- theta + proposal_SD * rnorm(p+J+1)
      
      log_r <- log_post_glmm(theta_prop) - log_post_glmm(theta)
      
      if (log(runif(1)) < log_r) {
        theta <- theta_prop
        n_acc <- n_acc + 1
        accepted[[n_acc]] <- theta
      }
      
      if (iter %% 10000 == 0) {
        cat("[MH] iter:", iter, " accepted:", n_acc, "/", target_accept, "\n")
        flush.console()
      }
    }
    
    list(
      samples = do.call(rbind, accepted),
      total_iterations = iter,
      acc_rate = target_accept / iter
    )
  }
  
  time_mh <- system.time({
    mh_out <- mh_sampler_acceptN(n_iter)
  })
  post_mh <- mh_out$samples[(burn+1):n_iter, 1:p]
  
  
  ###########################################################
  ## 3) Stan(brms) with backend="rstan"
  ###########################################################
  df_stan <- data.frame(
    y = y,
    X1 = X[,2],
    X2 = X[,3],
    X3 = X[,4],
    patientID = factor(apply(Z, 1, function(row) which(row==1)))
  )
  
  form <- bf(y ~ X1 + X2 + X3 + (1|patientID), family=bernoulli())
  
  time_stan <- system.time({
    fit_stan <- brm(
      form, data=df_stan,
      chains = 2,
      iter   = iter_stan,
      warmup = burn,
      cores  = 2,
      backend = "rstan",   # ★ cmdstanr → rstan
      refresh = 50
    )
  })
  
  post_stan <- as.matrix(fit_stan)[ , grep("^b_", colnames(as.matrix(fit_stan))) ]
  
  
  ###########################################################
  ## 4) Diagnostics
  ###########################################################
  ess_pg <- effectiveSize(as.mcmc(post_pg))
  ess_mh <- effectiveSize(as.mcmc(post_mh))
  ess_stan <- effectiveSize(as.mcmc(post_stan))
  
  ESR_pg_sec  <- ess_pg / time_pg[3]
  ESR_mh_sec  <- ess_mh / time_mh[3]
  ESR_stan_sec <- ess_stan / time_stan[3]
  
  ESR_pg_iter <- ess_pg / n_iter
  ESR_mh_iter <- ess_mh / mh_out$total_iterations
  ESR_stan_iter <- ess_stan / iter_stan
  
  
  list(
    ESS = list(PG = ess_pg, MH = ess_mh, Stan = ess_stan),
    ESR_sec = list(PG = ESR_pg_sec, MH = ESR_mh_sec, Stan = ESR_stan_sec),
    ESR_iter = list(PG = ESR_pg_iter, MH = ESR_mh_iter, Stan = ESR_stan_iter),
    Rhat_stan = rhat(fit_stan),
    PG_time = time_pg[3],
    MH_time = time_mh[3],
    Stan_time = time_stan[3],
    PG_beta_draws = post_pg,
    MH_beta_draws = post_mh,
    Stan_beta_draws = post_stan
  )
}

result <- compare_GLMM_3methods_rstan(y, X, Z)

result$ESS
result$ESR_sec
result$ESR_iter
result$Rhat_stan
result$Stan_time

par(mfrow=c(2,2))

plot(as.mcmc(result$PG_beta_draws[,1]), main="PG Trace: beta1")
plot(as.mcmc(result$PG_beta_draws[,2]), main="PG Trace: beta2")
plot(as.mcmc(result$PG_beta_draws[,3]), main="PG Trace: beta3")
plot(as.mcmc(result$PG_beta_draws[,4]), main="PG Trace: beta4")

plot(as.mcmc(result$MH_beta_draws[,1]), main="MH Trace: beta1")
plot(as.mcmc(result$MH_beta_draws[,2]), main="MH Trace: beta2")
plot(as.mcmc(result$MH_beta_draws[,3]), main="MH Trace: beta3")
plot(as.mcmc(result$MH_beta_draws[,4]), main="MH Trace: beta4")

plot(as.mcmc(result$Stan_beta_draws[,1]), main="Stan Trace: beta1")
plot(as.mcmc(result$Stan_beta_draws[,2]), main="Stan Trace: beta2")
plot(as.mcmc(result$Stan_beta_draws[,3]), main="Stan Trace: beta3")
plot(as.mcmc(result$Stan_beta_draws[,4]), main="Stan Trace: beta4")

###

acf(as.mcmc(result$PG_beta_draws[,1]), main="PG ACF: beta1")
acf(as.mcmc(result$PG_beta_draws[,2]), main="PG ACF: beta2")
acf(as.mcmc(result$PG_beta_draws[,3]), main="PG ACF: beta3")
acf(as.mcmc(result$PG_beta_draws[,4]), main="PG ACF: beta4")

acf(as.mcmc(result$MH_beta_draws[,1]), main="MH ACF: beta1")
acf(as.mcmc(result$MH_beta_draws[,2]), main="MH ACF: beta2")
acf(as.mcmc(result$MH_beta_draws[,3]), main="MH ACF: beta3")
acf(as.mcmc(result$MH_beta_draws[,4]), main="MH ACF: beta4")

acf(as.matrix(result$Stan_beta_draws[,1]), main="Stan ACF: beta1")
acf(as.matrix(result$Stan_beta_draws[,2]), main="Stan ACF: beta2")
acf(as.matrix(result$Stan_beta_draws[,3]), main="Stan ACF: beta3")
acf(as.matrix(result$Stan_beta_draws[,4]), main="Stan ACF: beta4")

### 

plot(density(result$PG_beta_draws[,1]), col="blue", lwd=2, main="beta1 Density")
lines(density(result$MH_beta_draws[,1]), col="red", lwd=2)
lines(density(result$Stan_beta_draws[,1]), col="green", lwd=2)
legend("topright", legend=c("PG","MH","Stan"), col=c("blue","red","green"), lwd=2)

plot(density(result$PG_beta_draws[,2]), col="blue", lwd=2, main="beta2 Density")
lines(density(result$MH_beta_draws[,2]), col="red", lwd=2)
lines(density(result$Stan_beta_draws[,2]), col="green", lwd=2)
legend("topright", legend=c("PG","MH","Stan"), col=c("blue","red","green"), lwd=2)

plot(density(result$PG_beta_draws[,3]), col="blue", lwd=2, main="beta3 Density")
lines(density(result$MH_beta_draws[,3]), col="red", lwd=2)
lines(density(result$Stan_beta_draws[,3]), col="green", lwd=2)
legend("topright", legend=c("PG","MH","Stan"), col=c("blue","red","green"), lwd=2)

plot(density(result$PG_beta_draws[,4]), col="blue", lwd=2, main="beta4 Density")
lines(density(result$MH_beta_draws[,4]), col="red", lwd=2)
lines(density(result$Stan_beta_draws[,4]), col="green", lwd=2)
legend("topright", legend=c("PG","MH","Stan"), col=c("blue","red","green"), lwd=2)

####

comp_table <- data.frame(
  Method = c("PG", "MH", "Stan"),
  Time_sec = c(result$PG_time, result$MH_time, result$Stan_time),
  ESS_min = c(min(result$ESS$PG), min(result$ESS$MH), min(result$ESS$Stan)),
  ESR_sec_min = c(min(result$ESR_sec$PG),
                  min(result$ESR_sec$MH),
                  min(result$ESR_sec$Stan)),
  ESR_iter_min = c(min(result$ESR_iter$PG),
                   min(result$ESR_iter$MH),
                   min(result$ESR_iter$Stan))
)

print(comp_table)




