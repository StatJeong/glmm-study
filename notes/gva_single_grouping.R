# 감마분포의 링크 함수가 뭔지도 확인해보고 
# 아마 shape, rate 두 모수를 잡는게 어려울거야 


suppressMessages(library(statmod))

## ============================================================
## 1. 데이터 생성
## ============================================================
simulate_data <- function(m, n, beta0, beta1, sigma2, family = c("poisson","gamma"), nu = 2) {
  family <- match.arg(family)
  X <- matrix(rnorm(m * n), nrow = m, ncol = n)
  U <- rnorm(m, 0, sqrt(sigma2))
  eta <- beta0 + beta1 * X + U   # U는 행(row)마다 재활용됨 -> 그룹 i에 U_i 더해짐
  mu  <- exp(eta)
  if (family == "poisson") {
    Y <- matrix(rpois(m * n, lambda = mu), nrow = m, ncol = n)
  } else {
    Y <- matrix(rgamma(m * n, shape = nu, rate = nu / mu), nrow = m, ncol = n)
  }
  list(X = X, Y = Y)
}

## ============================================================
## 2. GVA 추정량 (variational EM)
## ============================================================
.obj_grad_hess <- function(mu, lam, family, Bi_or_Di, c1, sigma2, nu = NULL) {
  if (family == "poisson") {
    E    <- Bi_or_Di * exp(mu + lam/2)
    gmu  <- c1 - E - mu/sigma2
    glam <- -0.5*E - 1/(2*sigma2) + 1/(2*lam)
    H11  <- -E - 1/sigma2
    H12  <- -0.5*E
    H22  <- -0.25*E - 1/(2*lam^2)
  } else {
    E    <- nu * Bi_or_Di * exp(-mu + lam/2)
    gmu  <- c1 + E - mu/sigma2
    glam <- -0.5*E - 1/(2*sigma2) + 1/(2*lam)
    H11  <- -E - 1/sigma2
    H12  <- 0.5*E
    H22  <- -0.25*E - 1/(2*lam^2)
  }
  list(gmu=gmu, glam=glam, H11=H11, H12=H12, H22=H22, E=E)
}

.newton_backtrack <- function(mu, lam, family, stat_i, c1, sigma2, nu=NULL, n_iter = 50) {
  obj_fun <- function(mu, lam) {
    if (family == "poisson") {
      E <- stat_i * exp(mu + lam/2)
      c1*mu - E - (mu^2+lam)/(2*sigma2) + 0.5*log(lam)
    } else {
      E <- nu * stat_i * exp(-mu + lam/2)
      c1*mu - E - (mu^2+lam)/(2*sigma2) + 0.5*log(lam)
    }
  }
  cur_obj <- obj_fun(mu, lam)
  for (it in 1:n_iter) {
    gh <- .obj_grad_hess(mu, lam, family, stat_i, c1, sigma2, nu)
    gmu <- gh$gmu
    glam <- gh$glam
    H11<-gh$H11; H12<-gh$H12; H22<-gh$H22
    det <- H11*H22 - H12^2
    use_newton <- is.finite(det) & (det > 1e-12) & (H11 < 0)
    dmu  <- numeric(length(mu)); dlam <- numeric(length(mu))
    dmu[use_newton]  <- -( H22[use_newton]*gmu[use_newton] - H12[use_newton]*glam[use_newton]) / det[use_newton]
    dlam[use_newton] <- -(-H12[use_newton]*gmu[use_newton] + H11[use_newton]*glam[use_newton]) / det[use_newton]
    step0 <- 0.05
    dmu[!use_newton]  <- step0*gmu[!use_newton]
    dlam[!use_newton] <- step0*glam[!use_newton]
    
    step <- rep(1, length(mu))
    for (bt in 1:30) {
      mu_try  <- mu  + step*dmu
      lam_try <- lam + step*dlam
      lam_try[lam_try <= 1e-8] <- NA
      new_obj <- obj_fun(mu_try, lam_try)
      improved <- is.finite(new_obj) & (new_obj >= cur_obj - 1e-10)
      if (all(improved)) break
      step[!improved] <- step[!improved] / 2
      if (max(step) < 1e-10) break
    }
    mu_try[is.na(mu_try)] <- mu[is.na(mu_try)]
    lam_try[is.na(lam_try) | lam_try<=1e-8] <- lam[is.na(lam_try) | lam_try<=1e-8]
    new_obj <- obj_fun(mu_try, lam_try)
    keep <- is.finite(new_obj) & (new_obj >= cur_obj - 1e-10)
    mu[keep]  <- mu_try[keep]
    lam[keep] <- lam_try[keep]
    cur_obj[keep] <- new_obj[keep]
    if (max(abs(step*dmu)) < 1e-9 && max(abs(step*dlam)) < 1e-9) break
  }
  list(mu = mu, lam = lam)
}

e_step_poisson <- function(mu, lam, Bi, Ybar, sigma2, n_iter = 50) {
  .newton_backtrack(mu, lam, "poisson", stat_i = Bi, c1 = Ybar, sigma2 = sigma2, n_iter = n_iter)
}

e_step_gamma <- function(mu, lam, Di, sigma2, nu, n, n_iter = 50) {
  .newton_backtrack(mu, lam, "gamma", stat_i = Di, c1 = -nu*n, sigma2 = sigma2, nu = nu, n_iter = n_iter)
}

fit_gva <- function(X, Y, family = c("poisson","gamma"), nu = 2, max_iter = 150, tol = 1e-7) {
  family <- match.arg(family)
  m <- nrow(X); n <- ncol(X)
  Ybar <- rowSums(Y)
  Xvec <- as.vector(X)
  Yvec <- as.vector(Y)
  
  fam_glm <- if (family == "poisson") poisson(link="log") else Gamma(link="log")
  init_fit <- suppressWarnings(glm(Yvec ~ Xvec, family = fam_glm))
  beta0 <- unname(coef(init_fit)[1]); beta1 <- unname(coef(init_fit)[2])
  sigma2 <- 0.3
  mu  <- rep(0, m)
  lam <- rep(0.1, m)
  
  for (it in 1:max_iter) {
    eta_fix <- beta0 + beta1 * X
    if (family == "poisson") {
      Bi <- rowSums(exp(eta_fix))
      es <- e_step_poisson(mu, lam, Bi, Ybar, sigma2)
    } else {
      Di <- rowSums(Y * exp(-eta_fix))
      es <- e_step_gamma(mu, lam, Di, sigma2, nu, n)
    }
    mu <- es$mu; lam <- es$lam
    
    offset_i <- if (family == "poisson") mu + lam/2 else -mu + lam/2
    offset_full <- rep(offset_i, times = n)
    
    bnew <- robust_beta_update(Xvec, Yvec, offset_full, family, nu, beta0, beta1)
    beta0_new <- unname(bnew["beta0"]); beta1_new <- unname(bnew["beta1"])
    sigma2_target <- mean(mu^2 + lam)
    damp <- 0.5
    sigma2_new <- sigma2 + damp*(sigma2_target - sigma2)
    
    delta <- abs(beta0_new-beta0) + abs(beta1_new-beta1) + abs(sigma2_new-sigma2)
    beta0 <- beta0_new; beta1 <- beta1_new; sigma2 <- sigma2_new
    if (delta < tol) break
  }
  c(beta0 = unname(beta0), beta1 = unname(beta1), sigma2 = sigma2, iters = it)
}

## ============================================================
## 3. Exact MLE (Adaptive Gauss-Hermite quadrature)
## ============================================================
exact_negloglik <- function(par, X, Y, family, nu = 2, n_quad = 11) {
  beta0 <- par[1]; beta1 <- par[2]; sigma2 <- exp(par[3])
  gh <- gauss.quad(n_quad, kind = "hermite")
  xk <- gh$nodes; wk <- gh$weights
  m <- nrow(X); n <- ncol(X)
  eta_fix <- beta0 + beta1 * X
  ETA_CLIP <- 30
  
  grad_hess_vec <- function(u_hat) {
    eta <- pmin(pmax(eta_fix + u_hat, -ETA_CLIP), ETA_CLIP)
    if (family == "poisson") {
      muv <- exp(eta)
      g <- rowSums(Y - muv) - u_hat/sigma2
      h <- -rowSums(muv) - 1/sigma2
    } else {
      ev <- exp(-eta)
      g <- rowSums(-nu + nu*Y*ev) - u_hat/sigma2
      h <- -rowSums(nu*Y*ev) - 1/sigma2
    }
    list(g=g, h=h)
  }
  
  u_hat <- rep(0, m)
  for (nt in 1:25) {
    gh_ <- grad_hess_vec(u_hat)
    g <- gh_$g; h <- gh_$h
    bad <- !is.finite(g) | !is.finite(h) | h >= -1e-10
    h[bad] <- -1; g[bad] <- 0
    step <- -g/h
    step[!is.finite(step)] <- 0
    step <- pmax(pmin(step, 5), -5)
    u_hat <- u_hat + step
    if (max(abs(step)) < 1e-10) break
  }
  h_at_mode <- grad_hess_vec(u_hat)$h
  h_at_mode[!is.finite(h_at_mode) | h_at_mode >= 0] <- -1
  s <- sqrt(-1/h_at_mode)
  
  ll_mat <- matrix(0, nrow = m, ncol = n_quad)
  for (k in 1:n_quad) {
    u_k <- u_hat + sqrt(2)*s*xk[k]
    eta_k <- pmin(pmax(eta_fix + u_k, -ETA_CLIP), ETA_CLIP)
    loglik_k <- if (family == "poisson") rowSums(dpois(Y, lambda = exp(eta_k), log = TRUE))
    else rowSums(dgamma(Y, shape = nu, rate = nu/exp(eta_k), log = TRUE))
    ll_mat[, k] <- loglik_k - u_k^2/(2*sigma2) + xk[k]^2
  }
  mx <- apply(ll_mat, 1, max)
  wsum <- as.vector(exp(ll_mat - mx) %*% wk)
  ll_total <- log(sqrt(2)*s) + mx + log(wsum) - 0.5*log(2*pi*sigma2)
  out <- -sum(ll_total)
  if (!is.finite(out)) out <- 1e10
  out
}

fit_exact <- function(X, Y, family, nu = 2, n_quad = 11, start = c(0,0,0)) {
  opt <- optim(start, exact_negloglik, X=X, Y=Y, family=family, nu=nu, n_quad=n_quad, method = "BFGS")
  c(beta0 = opt$par[1], beta1 = opt$par[2], sigma2 = exp(opt$par[3]))
}

## ============================================================
## 4. 시뮬레이션 드라이버
## ============================================================
run_pilot_sim <- function(m_grid, n_grid, reps, beta0=1, beta1=0.5, sigma2=1,
                          family = c("poisson","gamma"), nu = 2, n_quad = 11, verbose = TRUE) {
  family <- match.arg(family)
  results <- list()
  idx <- 1
  for (n in n_grid) for (m in m_grid) {
    for (r in 1:reps) {
      d <- simulate_data(m, n, beta0, beta1, sigma2, family, nu)
      g <- tryCatch(fit_gva(d$X, d$Y, family, nu), error = function(e) rep(NA,4))
      e <- tryCatch(fit_exact(d$X, d$Y, family, nu, n_quad), error = function(e) rep(NA,3))
      results[[idx]] <- data.frame(m=m, n=n, rep=r,
                                   gva_beta0=g["beta0"], gva_beta1=g["beta1"], gva_sigma2=g["sigma2"],
                                   exact_beta0=e["beta0"], exact_beta1=e["beta1"], exact_sigma2=e["sigma2"])
      idx <- idx + 1
    }
    if (verbose) cat(sprintf("[%s] m=%d, n=%d 완료\n", family, m, n))
  }
  out <- do.call(rbind, results)
  attr(out, "true_par") <- c(beta0=beta0, beta1=beta1, sigma2=sigma2)
  out
}

summarize_bias <- function(sim_out) {
  tp <- attr(sim_out, "true_par")
  agg <- aggregate(cbind(gva_beta0,gva_beta1,gva_sigma2,exact_beta0,exact_beta1,exact_sigma2) ~ m+n,
                   data = sim_out, FUN = function(x) mean(x, na.rm=TRUE))
  agg$gva_bias_beta0   <- agg$gva_beta0   - tp["beta0"]
  agg$gva_bias_sigma2  <- agg$gva_sigma2  - tp["sigma2"]
  agg$exact_bias_beta0 <- agg$exact_beta0 - tp["beta0"]
  agg$exact_bias_sigma2<- agg$exact_sigma2- tp["sigma2"]
  agg[order(agg$n, agg$m), ]
}

diagnose_convergence <- function(sim_out) {
  n_total <- nrow(sim_out)
  n_bad <- sum(!is.finite(sim_out$gva_beta0) | !is.finite(sim_out$exact_beta0))
  cat(sprintf("전체 %d개 중 fit 실패(NA/Inf) %d개 (%.1f%%)\n", n_total, n_bad, 100*n_bad/n_total))
  agg <- aggregate(cbind(gva_ok=is.finite(gva_beta0), exact_ok=is.finite(exact_beta0)) ~ m+n,
                   data=sim_out, FUN=function(x) sum(!x))
  names(agg)[3:4] <- c("gva_fail","exact_fail")
  agg
}

run_pilot_sim_parallel <- function(m_grid, n_grid, reps, beta0=1, beta1=0.5, sigma2=1,
                                   family = c("poisson","gamma"), nu = 2, n_quad = 11,
                                   mc.cores = parallel::detectCores(), base_seed = 1) {
  family <- match.arg(family)
  is_windows <- .Platform$OS.type == "windows"
  if (is_windows) {
    cl <- parallel::makeCluster(mc.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, varlist = c("simulate_data","fit_gva","fit_exact",
                                            "e_step_poisson","e_step_gamma",".newton_backtrack",".obj_grad_hess","exact_negloglik",
                                            "robust_beta_update"), envir = globalenv())
  }
  
  one_rep <- function(rr, mm, nn) {
    set.seed(base_seed*100000 + mm*1000 + nn*100 + rr)
    d <- simulate_data(mm, nn, beta0, beta1, sigma2, family, nu)
    g <- tryCatch(fit_gva(d$X, d$Y, family, nu), error = function(e) rep(NA,4))
    e <- tryCatch(fit_exact(d$X, d$Y, family, nu, n_quad), error = function(e) rep(NA,3))
    data.frame(m=mm, n=nn, rep=rr,
               gva_beta0=g["beta0"], gva_beta1=g["beta1"], gva_sigma2=g["sigma2"],
               exact_beta0=e["beta0"], exact_beta1=e["beta1"], exact_sigma2=e["sigma2"])
  }
  
  cat(sprintf("코어 %d개 사용. (m,n) 조합 %d개 x reps %d개\n", mc.cores, length(m_grid)*length(n_grid), reps))
  all_results <- list(); grp_idx <- 1
  t_start_all <- Sys.time()
  for (n in n_grid) for (m in m_grid) {
    t0 <- Sys.time()
    if (is_windows) {
      res_list <- parallel::parLapply(cl, 1:reps, one_rep, mm=m, nn=n)
    } else {
      res_list <- parallel::mclapply(1:reps, one_rep, mm=m, nn=n, mc.cores = mc.cores,
                                     mc.preschedule = FALSE)
    }
    all_results[[grp_idx]] <- do.call(rbind, res_list)
    grp_idx <- grp_idx + 1
    elapsed <- as.numeric(Sys.time()-t0, units="secs")
    total_elapsed <- as.numeric(Sys.time()-t_start_all, units="secs")
    cat(sprintf("[%s] m=%-7d n=%-3d 완료  (이 조합 %.1fs, 누적 %.1fs)\n",
                family, m, n, elapsed, total_elapsed))
  }
  out <- do.call(rbind, all_results)
  attr(out, "true_par") <- c(beta0=beta0, beta1=beta1, sigma2=sigma2)
  out
}

## glm()의 IRLS가 sigma2가 클 때 가끔 죽는 문제를 우회하는, 직접 짠 안전한 beta 업데이트
robust_beta_update <- function(Xvec, Yvec, offset_full, family, nu, beta0, beta1, n_iter=50) {
  beta0 <- unname(beta0); beta1 <- unname(beta1)
  obj_and_gh <- function(b0, b1) {
    eta <- b0 + b1*Xvec + offset_full
    if (family == "poisson") {
      mu <- exp(eta)
      obj <- sum(Yvec*eta - mu)
      g0 <- sum(Yvec - mu);           g1 <- sum((Yvec - mu)*Xvec)
      h00<- -sum(mu);                 h01<- -sum(mu*Xvec);  h11<- -sum(mu*Xvec^2)
    } else {
      ev <- exp(-eta)
      obj <- sum(-eta - Yvec*ev)
      g0 <- sum(-1 + Yvec*ev);        g1 <- sum((-1 + Yvec*ev)*Xvec)
      h00<- -sum(Yvec*ev);            h01<- -sum(Yvec*ev*Xvec); h11<- -sum(Yvec*ev*Xvec^2)
    }
    list(obj=obj, g0=g0, g1=g1, h00=h00, h01=h01, h11=h11)
  }
  cur <- obj_and_gh(beta0, beta1); cur_obj <- cur$obj
  for (it in 1:n_iter) {
    det <- cur$h00*cur$h11 - cur$h01^2
    if (!is.finite(det) || det <= 1e-12 || cur$h00 >= 0) break
    d0 <- -( cur$h11*cur$g0 - cur$h01*cur$g1) / det
    d1 <- -(-cur$h01*cur$g0 + cur$h00*cur$g1) / det
    step <- 1
    for (bt in 1:40) {
      b0_try <- beta0+step*d0; b1_try <- beta1+step*d1
      cand <- obj_and_gh(b0_try, b1_try)
      if (is.finite(cand$obj) && cand$obj >= cur_obj - 1e-10) break
      step <- step/2
    }
    beta0 <- beta0+step*d0; beta1 <- beta1+step*d1
    cur <- obj_and_gh(beta0, beta1)
    if (abs(step*d0) < 1e-10 && abs(step*d1) < 1e-10) break
    cur_obj <- cur$obj
  }
  c(beta0=beta0, beta1=beta1)
}

## 완전한 population ELBO -- 여러 시작점 중 진짜 승자를 가리는 데 씀
full_elbo <- function(beta0, beta1, sigma2, X, Y, family, nu) {
  m <- nrow(X); n <- ncol(X)
  eta_fix <- beta0 + beta1*X
  if (family == "poisson") {
    stat_i <- rowSums(exp(eta_fix)); c1 <- rowSums(Y)
  } else {
    stat_i <- rowSums(Y*exp(-eta_fix)); c1 <- -nu*n
  }
  mu <- rep(0,m); lam <- rep(0.1,m)
  es <- .newton_backtrack(mu, lam, family, stat_i, c1, sigma2, nu=nu, n_iter=100)
  mu<-es$mu; lam<-es$lam
  E <- if (family=="poisson") stat_i*exp(mu+lam/2) else nu*stat_i*exp(-mu+lam/2)
  per_group <- c1*mu - E - (mu^2+lam)/(2*sigma2) - 0.5*log(sigma2) + 0.5*log(lam)
  extra <- if (family=="poisson") sum(Y*eta_fix) else -nu*sum(eta_fix)
  sum(per_group) + extra
}

## fit_gva와 같은 EM 구조인데, 여러 시작점에서 돌려서 ELBO가 진짜 제일 높은 지점을 골라줌
fit_gva_multistart <- function(X, Y, family = c("poisson","gamma"), nu = 2,
                               n_starts = 6, seed_true = NULL, max_iter = 200, tol = 1e-8) {
  family <- match.arg(family)
  m <- nrow(X); n <- ncol(X)
  Ybar <- rowSums(Y)
  Xvec <- as.vector(X); Yvec <- as.vector(Y)
  fam_glm <- if (family == "poisson") poisson(link="log") else Gamma(link="log")
  init_fit <- suppressWarnings(glm(Yvec ~ Xvec, family = fam_glm))
  b0n <- unname(coef(init_fit)[1]); b1n <- unname(coef(init_fit)[2])
  
  starts <- list(
    c(beta0=b0n,      beta1=b1n,      sigma2=0.3),
    c(beta0=b0n*0.5,  beta1=b1n,      sigma2=0.1),
    c(beta0=b0n*1.5,  beta1=b1n,      sigma2=0.5),
    c(beta0=b0n,      beta1=b1n*0.5,  sigma2=1.0),
    c(beta0=b0n,      beta1=b1n,      sigma2=2.0),
    c(beta0=0,        beta1=0,        sigma2=1.0)
  )
  if (!is.null(seed_true)) starts[[length(starts)+1]] <- c(beta0=seed_true[1], beta1=seed_true[2], sigma2=seed_true[3])
  
  run_one <- function(st) {
    beta0<-unname(st["beta0"]); beta1<-unname(st["beta1"]); sigma2<-unname(st["sigma2"])
    mu <- rep(0, m); lam <- rep(0.1, m)
    for (it in 1:max_iter) {
      eta_fix <- beta0 + beta1 * X
      if (family == "poisson") { Bi<-rowSums(exp(eta_fix)); es<-e_step_poisson(mu,lam,Bi,Ybar,sigma2) }
      else { Di<-rowSums(Y*exp(-eta_fix)); es<-e_step_gamma(mu,lam,Di,sigma2,nu,n) }
      mu<-es$mu; lam<-es$lam
      offset_i <- if (family=="poisson") mu+lam/2 else -mu+lam/2
      offset_full <- rep(offset_i, times=n)
      bnew <- robust_beta_update(Xvec, Yvec, offset_full, family, nu, beta0, beta1)
      beta0_new<-unname(bnew["beta0"]); beta1_new<-unname(bnew["beta1"])
      sigma2_new <- sigma2 + 0.5*(mean(mu^2+lam)-sigma2)
      delta <- abs(beta0_new-beta0)+abs(beta1_new-beta1)+abs(sigma2_new-sigma2)
      beta0<-beta0_new; beta1<-beta1_new; sigma2<-sigma2_new
      if (delta<tol) break
    }
    elbo_val <- full_elbo(beta0, beta1, sigma2, X, Y, family, nu)
    c(beta0=beta0, beta1=beta1, sigma2=sigma2, iters=it, elbo=elbo_val)
  }
  
  results <- t(sapply(starts, run_one))
  out <- as.data.frame(results)
  out$start_id <- seq_len(nrow(out))
  out[order(-out$elbo), c("start_id","beta0","beta1","sigma2","iters","elbo")]
}

## 기존 단일경로 grid 결과를 multi-start로 검증할 때 씀
verify_with_multistart <- function(m_grid, n_grid, reps, beta0=1, beta1=0.5, sigma2=1,
                                   family = c("poisson","gamma"), nu = 2, n_starts_seed_true = TRUE) {
  family <- match.arg(family)
  results <- list(); idx <- 1
  for (n in n_grid) for (m in m_grid) {
    t0 <- Sys.time()
    for (r in 1:reps) {
      set.seed(900000 + m*1000 + n*100 + r)
      d <- simulate_data(m, n, beta0, beta1, sigma2, family, nu)
      ms <- tryCatch(
        fit_gva_multistart(d$X, d$Y, family, nu,
                           seed_true = if (n_starts_seed_true) c(beta0,beta1,sigma2) else NULL),
        error = function(e) NULL)
      if (is.null(ms)) {
        best <- list(beta0=NA, beta1=NA, sigma2=NA)
      } else {
        best <- list(beta0=ms$beta0[1], beta1=ms$beta1[1], sigma2=ms$sigma2[1])
      }
      results[[idx]] <- data.frame(m=m, n=n, rep=r,
                                   ms_beta0=best$beta0, ms_beta1=best$beta1, ms_sigma2=best$sigma2)
      idx <- idx + 1
    }
    cat(sprintf("[%s multistart] m=%d n=%d 완료 (%.1fs)\n", family, m, n,
                as.numeric(Sys.time()-t0, units="secs")))
  }
  out <- do.call(rbind, results)
  attr(out, "true_par") <- c(beta0=beta0, beta1=beta1, sigma2=sigma2)
  out
}

summarize_multistart_bias <- function(sim_out) {
  tp <- attr(sim_out, "true_par")
  agg <- aggregate(cbind(ms_beta0,ms_beta1,ms_sigma2) ~ m+n, data=sim_out,
                   FUN=function(x) mean(x, na.rm=TRUE))
  agg$ms_bias_beta0  <- agg$ms_beta0  - tp["beta0"]
  agg$ms_bias_sigma2 <- agg$ms_sigma2 - tp["sigma2"]
  agg[order(agg$n, agg$m), ]
}

res <- run_pilot_sim_parallel(m_grid=c(1000,10000), n_grid=c(5), reps=50, family="gamma", nu=2)
print(summarize_bias(res))




res_s2 <- run_pilot_sim_parallel(m_grid=c(500,1000,10000), n_grid=c(5), reps=50,
                                 family="gamma", nu=2, beta0=1, beta1=0.5, sigma2=2)
print(summarize_bias(res_s2))
diagnose_convergence(res_s2)


## m을 더 키워서 진짜 plateau인지 확인
res_s2_deep <- run_pilot_sim_parallel(m_grid=c(10000,30000,100000), n_grid=c(5), reps=30,
                                      family="gamma", nu=2, beta0=1, beta1=0.5, sigma2=2)
print(summarize_bias(res_s2_deep))

