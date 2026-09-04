###############################################################################
##  glmm_multi.R
##  다차원 랜덤효과 Gamma/Poisson GLMM: GVA vs AGHQ
##
##  모형:  eta_ij = x_ij' beta + z_ij' U_i,   U_i ~ N(0, Sigma),  dim(U_i) = d
##         Y_ij | U_i ~ Gamma(shape = nu, rate = nu/exp(eta_ij))   (log link)
##
##  d = 1 이면 랜덤절편만인 기존 모형과 동일하다.
##  d = 2 이면 랜덤절편 + 랜덤기울기.
##
##  ---------------------------------------------------------------------------
##  GVA 유도 (d차원)
##
##  q(U_i) = N(m_i, Lam_i) 로 두면, 다변량 정규 MGF
##      E_q[exp(-z'U_i)] = exp(-z'm_i + z'Lam_i z / 2)
##  덕분에 ELBO가 닫힌 형태로 남는다(수치적분 불필요):
##
##    ELBO_i = sum_j [ -nu*(eta_fix_ij + z_ij'm_i) - nu*Y_ij*E_ij ]
##             - 0.5*( tr(Sinv %*% Lam_i) + m_i'Sinv m_i - d
##                     - log|Lam_i| + log|Sigma| )
##    E_ij   = exp( -eta_fix_ij - z_ij'm_i + z_ij'Lam_i z_ij / 2 )
##
##  1차조건:
##    dELBO/dLam_i = 0  ->  Lam_i^{-1} = Sinv + sum_j nu*Y_ij*E_ij* z_ij z_ij'
##    dELBO/dm_i   = 0  ->  sum_j z_ij*(-nu + nu*Y_ij*E_ij) - Sinv m_i = 0
##    (m_i 의 헤시안은 정확히 -Lam_i^{-1} 이라 뉴턴 스텝이 자연스럽다)
##
##  M-step:
##    Sigma = (1/m) sum_i ( m_i m_i' + Lam_i )                     (닫힌 형태)
##    beta  = offset-GLM,  offset_ij = z_ij'm_i - z_ij'Lam_i z_ij/2
##
##  d=1 이면 offset은 mu_i - lam_i/2 로 환원된다.
###############################################################################

if (!requireNamespace("statmod", quietly = TRUE))
  stop("statmod 패키지가 필요합니다:  install.packages('statmod')")

ETA_CLIP <- 30

## ===========================================================================
## 1. 데이터 생성
## ===========================================================================
## d = 1: z_ij = 1                       (랜덤절편)
## d = 2: z_ij = (1, X1_ij)              (절편 + 기울기)
## d = 3: z_ij = (1, X1_ij, X2_ij)       ...
## 고정효과도 같은 공변량을 쓴다: eta = beta0 + sum_k beta_k X_k + z'U
simulate_data_multi <- function(m, n, beta, Sigma, nu = 2,
                                family = c("gamma", "poisson")) {
  family <- match.arg(family)
  d <- nrow(Sigma)
  p <- length(beta)                 # beta = (beta0, beta1, ..., beta_{p-1})
  ncov <- max(p - 1, d - 1)         # 필요한 공변량 개수

  Xl <- if (ncov > 0)
    lapply(seq_len(ncov), function(k) matrix(rnorm(m * n), m, n)) else list()

  ## 고정부: beta0 + beta1*X1 + ...
  eta <- matrix(beta[1], m, n)
  if (p > 1) for (k in seq_len(p - 1)) eta <- eta + beta[k + 1] * Xl[[k]]

  ## 랜덤부: z = (1, X1, ..., X_{d-1})
  U <- matrix(rnorm(m * d), m, d) %*% chol(Sigma)      # m x d,  행별 N(0,Sigma)
  eta <- eta + U[, 1]
  if (d > 1) for (k in seq_len(d - 1)) eta <- eta + U[, k + 1] * Xl[[k]]

  eta <- pmin(pmax(eta, -ETA_CLIP), ETA_CLIP)
  mu <- exp(eta)
  Y <- if (family == "poisson") matrix(rpois(m * n, mu), m, n)
       else matrix(rgamma(m * n, shape = nu, rate = nu / mu), m, n)

  list(Y = Y, Xl = Xl, U = U, d = d, family = family, nu = nu)
}

## 랜덤효과 설계행렬 목록 Zl[[k]] (각각 m x n) 을 만든다: (1, X1, ..., X_{d-1})
make_Z <- function(Xl, d, m, n) {
  Zl <- vector("list", d)
  Zl[[1]] <- matrix(1, m, n)
  if (d > 1) for (k in seq_len(d - 1)) Zl[[k + 1]] <- Xl[[k]]
  Zl
}
make_Xdesign <- function(Xl, p, m, n) {
  Dl <- vector("list", p)
  Dl[[1]] <- matrix(1, m, n)
  if (p > 1) for (k in seq_len(p - 1)) Dl[[k + 1]] <- Xl[[k]]
  Dl
}

## ===========================================================================
## 2. 배치 선형대수 헬퍼
##    그룹별 d x d 대칭행렬 m개를 (m x d x d) 배열로 다룬다.
## ===========================================================================
## 대칭 양정부호 배열의 역행렬 (d=1,2는 닫힌 형태, d>=3은 루프)
batch_solve_sym <- function(A) {
  m <- dim(A)[1]; d <- dim(A)[2]
  if (d == 1) {
    out <- A; out[, 1, 1] <- 1 / A[, 1, 1]; return(out)
  }
  if (d == 2) {
    a <- A[, 1, 1]; b <- A[, 1, 2]; c <- A[, 2, 2]
    det <- a * c - b * b
    det[!is.finite(det) | det <= 0] <- NA_real_
    out <- array(0, dim(A))
    out[, 1, 1] <-  c / det; out[, 2, 2] <-  a / det
    out[, 1, 2] <- -b / det; out[, 2, 1] <- -b / det
    return(out)
  }
  out <- array(0, dim(A))
  for (i in seq_len(m))
    out[i, , ] <- tryCatch(chol2inv(chol(A[i, , ])),
                           error = function(e) matrix(NA_real_, d, d))
  out
}

batch_logdet_sym <- function(A) {
  m <- dim(A)[1]; d <- dim(A)[2]
  if (d == 1) return(log(A[, 1, 1]))
  if (d == 2) return(log(A[, 1, 1] * A[, 2, 2] - A[, 1, 2]^2))
  vapply(seq_len(m), function(i)
    tryCatch(2 * sum(log(diag(chol(A[i, , ])))), error = function(e) NA_real_),
    numeric(1))
}

## 배열 A (m x d x d) 와 벡터모음 v (m x d) 의 그룹별 곱
batch_mv <- function(A, v) {
  d <- dim(A)[2]
  out <- matrix(0, nrow(v), d)
  for (k in seq_len(d)) {
    s <- 0
    for (l in seq_len(d)) s <- s + A[, k, l] * v[, l]
    out[, k] <- s
  }
  out
}

## ===========================================================================
## 3. GVA
## ===========================================================================
## E-step: 그룹별 (m_i, Lam_i) 를 1차조건으로 푼다.
##   Lam_i <- (Sinv + sum_j c_ij z z')^{-1}   (자동으로 양정부호)
##   m_i   <- 뉴턴 스텝 (헤시안 = -Lam_i^{-1})
gva_estep <- function(Mi, Lam, Zl, Y, eta_fix, Sinv, nu, family,
                      n_iter = 50, tol = 1e-10) {
  m <- nrow(Y); d <- length(Zl)
  for (it in seq_len(n_iter)) {
    ## zm_ij = z_ij' m_i,   zLz_ij = z_ij' Lam_i z_ij
    zm <- 0; for (k in seq_len(d)) zm <- zm + Zl[[k]] * Mi[, k]
    zLz <- 0
    for (k in seq_len(d)) for (l in seq_len(d))
      zLz <- zLz + Zl[[k]] * Zl[[l]] * Lam[, k, l]

    if (family == "poisson") {
      ## Poisson: E[exp(+z'U)],  data term  Y*eta - exp(eta)
      A <- exp(pmin(eta_fix + zm + zLz / 2, ETA_CLIP))     # = mu_ij 의 기댓값
      cij <- A                                             # z z' 계수
      gcoef <- Y - A                                       # dE/d(z'm) 계수
    } else {
      E <- exp(pmin(-eta_fix - zm + zLz / 2, ETA_CLIP))
      cij <- nu * Y * E
      gcoef <- -nu + nu * Y * E
    }

    ## Lam_i^{-1} = Sinv + sum_j cij z z'
    Ainv <- array(0, c(m, d, d))
    for (k in seq_len(d)) for (l in k:d) {
      s <- rowSums(cij * Zl[[k]] * Zl[[l]])
      Ainv[, k, l] <- Sinv[k, l] + s
      if (l != k) Ainv[, l, k] <- Ainv[, k, l]
    }
    Lam_new <- batch_solve_sym(Ainv)
    ok <- !is.na(Lam_new[, 1, 1])
    if (any(ok)) Lam[ok, , ] <- Lam_new[ok, , ]

    ## 기울기: sum_j z_ij * gcoef_ij  -  Sinv m_i
    G <- matrix(0, m, d)
    for (k in seq_len(d)) G[, k] <- rowSums(gcoef * Zl[[k]])
    G <- G - batch_mv(array(rep(Sinv, each = m), c(m, d, d)), Mi)

    ## 뉴턴: step = Lam_i %*% G   (헤시안 = -Lam_i^{-1})
    step <- batch_mv(Lam, G)
    step[!is.finite(step)] <- 0
    step <- pmax(pmin(step, 5), -5)
    Mi <- Mi + step
    if (max(abs(step)) < tol) break
  }
  list(Mi = Mi, Lam = Lam)
}

## beta M-step: offset-GLM 을 2p-변수 뉴턴으로 푼다 (백트래킹 포함)
beta_update_multi <- function(Dl, Yv, offset_v, beta, family, nu, n_iter = 50) {
  p <- length(beta)
  Dv <- lapply(Dl, as.vector)
  obj_gh <- function(b) {
    eta <- offset_v
    for (k in seq_len(p)) eta <- eta + b[k] * Dv[[k]]
    eta <- pmin(pmax(eta, -ETA_CLIP), ETA_CLIP)
    if (family == "poisson") {
      mu <- exp(eta); r <- Yv - mu; w <- mu
      obj <- sum(Yv * eta - mu)
    } else {
      ev <- exp(-eta); r <- -1 + Yv * ev; w <- Yv * ev
      obj <- sum(-eta - Yv * ev)
    }
    g <- vapply(seq_len(p), function(k) sum(r * Dv[[k]]), numeric(1))
    H <- matrix(0, p, p)
    for (k in seq_len(p)) for (l in k:p) {
      H[k, l] <- -sum(w * Dv[[k]] * Dv[[l]]); H[l, k] <- H[k, l]
    }
    list(obj = obj, g = g, H = H)
  }
  cur <- obj_gh(beta); cur_obj <- cur$obj
  for (it in seq_len(n_iter)) {
    step <- tryCatch(solve(cur$H, -cur$g), error = function(e) rep(0, p))
    if (!all(is.finite(step)) || all(step == 0)) break
    s <- 1; moved <- FALSE
    for (bt in 1:40) {
      cand <- obj_gh(beta + s * step)
      if (is.finite(cand$obj) && cand$obj >= cur_obj - 1e-10) { moved <- TRUE; break }
      s <- s / 2
    }
    if (!moved) break
    beta <- beta + s * step
    cur <- obj_gh(beta)
    if (max(abs(s * step)) < 1e-10) break
    cur_obj <- cur$obj
  }
  beta
}

## GVA 전체 적합
fit_gva_multi <- function(dat, p = NULL, max_iter = 3000, tol = 1e-9,
                          verbose = FALSE) {
  Y <- dat$Y; family <- dat$family; nu <- dat$nu; d <- dat$d
  m <- nrow(Y); n <- ncol(Y)
  if (is.null(p)) p <- length(dat$Xl) + 1
  Zl <- make_Z(dat$Xl, d, m, n)
  Dl <- make_Xdesign(dat$Xl, p, m, n)
  Yv <- as.vector(Y); Dv <- lapply(Dl, as.vector)

  ## 초기값: 랜덤효과를 무시한 주변 적합.
  ## Gamma+log 의 glm 은 시작값 문제로 "cannot correct step size" 로 자주 죽는다.
  ## log(Y) 에 대한 OLS 는 log link 하에서 항상 안정적이고 초기값으로 충분하다.
  Dmat <- do.call(cbind, Dv)
  beta <- tryCatch({
    z <- if (family == "poisson") log(pmax(Yv, 0.5)) else log(pmax(Yv, 1e-8))
    unname(qr.solve(Dmat, z))
  }, error = function(e) rep(0, p))
  beta[!is.finite(beta)] <- 0

  Sigma <- diag(0.3, d)
  Mi <- matrix(0, m, d)
  Lam <- array(0, c(m, d, d)); for (k in seq_len(d)) Lam[, k, k] <- 0.1

  for (it in seq_len(max_iter)) {
    Sinv <- tryCatch(chol2inv(chol(Sigma)), error = function(e) diag(1e6, d))
    eta_fix <- 0; for (k in seq_len(p)) eta_fix <- eta_fix + beta[k] * Dl[[k]]

    es <- gva_estep(Mi, Lam, Zl, Y, eta_fix, Sinv, nu, family)
    Mi <- es$Mi; Lam <- es$Lam

    ## offset_ij = z'm_i - z'Lam_i z/2   (poisson 은 +z'm + z'Lam z/2)
    zm <- 0; for (k in seq_len(d)) zm <- zm + Zl[[k]] * Mi[, k]
    zLz <- 0
    for (k in seq_len(d)) for (l in seq_len(d))
      zLz <- zLz + Zl[[k]] * Zl[[l]] * Lam[, k, l]
    off <- if (family == "poisson") zm + zLz / 2 else zm - zLz / 2

    beta_new <- beta_update_multi(Dl, Yv, as.vector(off), beta, family, nu)

    ## Sigma = mean_i ( m_i m_i' + Lam_i )
    S_new <- matrix(0, d, d)
    for (k in seq_len(d)) for (l in seq_len(d))
      S_new[k, l] <- mean(Mi[, k] * Mi[, l] + Lam[, k, l])
    S_new <- (S_new + t(S_new)) / 2

    delta <- max(abs(beta_new - beta), abs(S_new - Sigma))
    beta <- beta_new; Sigma <- S_new
    if (verbose && it %% 100 == 0)
      cat(sprintf("  it %4d  delta %.3e  beta0 %.4f\n", it, delta, beta[1]))
    if (delta < tol) break
  }
  list(beta = beta, Sigma = Sigma, iters = it, converged = it < max_iter,
       Mi = Mi, Lam = Lam)
}

## GVA의 ELBO (수렴 품질 확인용) -- 주어진 (beta, Sigma)에서 E-step을 다시 풀어 평가
elbo_multi <- function(beta, Sigma, dat, p = NULL, n_iter = 200) {
  Y <- dat$Y; family <- dat$family; nu <- dat$nu; d <- dat$d
  m <- nrow(Y); n <- ncol(Y)
  if (is.null(p)) p <- length(dat$Xl) + 1
  Zl <- make_Z(dat$Xl, d, m, n); Dl <- make_Xdesign(dat$Xl, p, m, n)
  Sinv <- tryCatch(chol2inv(chol(Sigma)), error = function(e) return(NA_real_))
  if (length(Sinv) == 1 && is.na(Sinv)) return(NA_real_)
  eta_fix <- 0; for (k in seq_len(p)) eta_fix <- eta_fix + beta[k] * Dl[[k]]
  Mi <- matrix(0, m, d); Lam <- array(0, c(m, d, d))
  for (k in seq_len(d)) Lam[, k, k] <- 0.1
  es <- gva_estep(Mi, Lam, Zl, Y, eta_fix, Sinv, nu, family, n_iter = n_iter)
  Mi <- es$Mi; Lam <- es$Lam

  zm <- 0; for (k in seq_len(d)) zm <- zm + Zl[[k]] * Mi[, k]
  zLz <- 0
  for (k in seq_len(d)) for (l in seq_len(d))
    zLz <- zLz + Zl[[k]] * Zl[[l]] * Lam[, k, l]

  data_term <- if (family == "poisson") {
    A <- exp(pmin(eta_fix + zm + zLz / 2, ETA_CLIP))
    sum(Y * (eta_fix + zm) - A)
  } else {
    E <- exp(pmin(-eta_fix - zm + zLz / 2, ETA_CLIP))
    sum(-nu * (eta_fix + zm) - nu * Y * E)
  }
  ## KL 항
  trSL <- 0
  for (k in seq_len(d)) for (l in seq_len(d)) trSL <- trSL + Sinv[k, l] * Lam[, l, k]
  mSm <- rowSums(batch_mv(array(rep(Sinv, each = m), c(m, d, d)), Mi) * Mi)
  logdetLam <- batch_logdet_sym(Lam)
  kl <- 0.5 * sum(trSL + mSm - d - logdetLam + determinant(Sigma, logarithm = TRUE)$modulus[1])
  data_term - kl
}

## ===========================================================================
## 4. AGHQ  (adaptive Gauss-Hermite quadrature, d차원)
##    그룹별로 최빈값 u_hat 과 헤시안을 찾고, 그 주위에 q^d 격자를 깐다.
##    비용이 q^d 로 지수적으로 늘어나는 것이 GVA와의 핵심 차이다.
## ===========================================================================
aghq_negloglik <- function(par, dat, p, d, n_quad = 11, gh = NULL) {
  Y <- dat$Y; family <- dat$family; nu <- dat$nu
  m <- nrow(Y); n <- ncol(Y)
  beta <- par[seq_len(p)]
  ## Sigma 는 log-Cholesky 로 모수화 (양정부호 보장)
  L <- matrix(0, d, d); idx <- p
  for (k in seq_len(d)) for (l in seq_len(k)) {
    idx <- idx + 1
    L[k, l] <- if (k == l) exp(par[idx]) else par[idx]
  }
  Sigma <- L %*% t(L)
  Sinv <- tryCatch(chol2inv(chol(Sigma)), error = function(e) NULL)
  if (is.null(Sinv)) return(1e10)
  logdetS <- 2 * sum(log(diag(L)))

  Zl <- make_Z(dat$Xl, d, m, n); Dl <- make_Xdesign(dat$Xl, p, m, n)
  eta_fix <- 0; for (k in seq_len(p)) eta_fix <- eta_fix + beta[k] * Dl[[k]]

  ## --- 그룹별 최빈값 찾기 (뉴턴) ---
  Uh <- matrix(0, m, d)
  gh_of <- function(Uh) {
    zu <- 0; for (k in seq_len(d)) zu <- zu + Zl[[k]] * Uh[, k]
    eta <- pmin(pmax(eta_fix + zu, -ETA_CLIP), ETA_CLIP)
    if (family == "poisson") { A <- exp(eta); gc_ <- Y - A; wc <- A }
    else { ev <- exp(-eta); gc_ <- -nu + nu * Y * ev; wc <- nu * Y * ev }
    G <- matrix(0, m, d)
    for (k in seq_len(d)) G[, k] <- rowSums(gc_ * Zl[[k]])
    G <- G - batch_mv(array(rep(Sinv, each = m), c(m, d, d)), Uh)
    H <- array(0, c(m, d, d))
    for (k in seq_len(d)) for (l in k:d) {
      s <- rowSums(wc * Zl[[k]] * Zl[[l]])
      H[, k, l] <- Sinv[k, l] + s
      if (l != k) H[, l, k] <- H[, k, l]
    }
    list(G = G, Hinv_neg = H)   # -Hessian (양정부호)
  }
  for (nt in 1:30) {
    r <- gh_of(Uh)
    Cov <- batch_solve_sym(r$Hinv_neg)
    step <- batch_mv(Cov, r$G)
    step[!is.finite(step)] <- 0
    step <- pmax(pmin(step, 5), -5)
    Uh <- Uh + step
    if (max(abs(step)) < 1e-10) break
  }
  r <- gh_of(Uh)
  Cov <- batch_solve_sym(r$Hinv_neg)
  if (any(is.na(Cov))) return(1e10)

  ## Cov_i = L_i L_i'  (그룹별 촐레스키)
  Lc <- array(0, c(m, d, d))
  for (i in seq_len(m)) {
    ci <- tryCatch(t(chol(Cov[i, , ])), error = function(e) NULL)
    if (is.null(ci)) return(1e10)
    Lc[i, , ] <- ci
  }
  logdetL <- rowSums(log(vapply(seq_len(d), function(k) Lc[, k, k], numeric(m))))

  ## --- q^d 격자 ---
  if (is.null(gh)) gh <- statmod::gauss.quad(n_quad, kind = "hermite")
  grid <- as.matrix(expand.grid(rep(list(seq_len(n_quad)), d)))
  ll_mat <- matrix(-Inf, m, nrow(grid))
  for (g in seq_len(nrow(grid))) {
    tk <- gh$nodes[grid[g, ]]                  # 길이 d
    lw <- sum(log(gh$weights[grid[g, ]]))
    ## u = u_hat + sqrt(2) * L_i %*% t
    uk <- Uh
    for (k in seq_len(d)) {
      s <- 0
      for (l in seq_len(d)) s <- s + Lc[, k, l] * tk[l]
      uk[, k] <- Uh[, k] + sqrt(2) * s
    }
    zu <- 0; for (k in seq_len(d)) zu <- zu + Zl[[k]] * uk[, k]
    eta <- pmin(pmax(eta_fix + zu, -ETA_CLIP), ETA_CLIP)
    ll <- if (family == "poisson") rowSums(dpois(Y, exp(eta), log = TRUE))
          else rowSums(dgamma(Y, shape = nu, rate = nu / exp(eta), log = TRUE))
    quad <- rowSums(batch_mv(array(rep(Sinv, each = m), c(m, d, d)), uk) * uk)
    ll_mat[, g] <- ll - 0.5 * quad + sum(tk^2) + lw
  }
  mx <- apply(ll_mat, 1, max)
  wsum <- rowSums(exp(ll_mat - mx))
  ## log ∫ = log( 2^{d/2} |L_i| ) + mx + log(wsum) - 0.5*log|Sigma| - (d/2)log(2pi)
  ll_total <- (d / 2) * log(2) + logdetL + mx + log(wsum) -
              0.5 * logdetS - (d / 2) * log(2 * pi)
  out <- -sum(ll_total)
  if (!is.finite(out)) out <- 1e10
  out
}

## AGHQ 적합. 시작값을 GVA 결과로 주면 발산을 크게 줄일 수 있다.
fit_aghq_multi <- function(dat, p = NULL, n_quad = 11, start = NULL,
                           gva_start = TRUE, maxit = 300) {
  d <- dat$d
  if (is.null(p)) p <- length(dat$Xl) + 1
  gh <- statmod::gauss.quad(n_quad, kind = "hermite")

  if (is.null(start)) {
    if (gva_start) {
      g <- tryCatch(fit_gva_multi(dat, p, max_iter = 300, tol = 1e-6),
                    error = function(e) NULL)
      if (!is.null(g)) {
        S <- g$Sigma
        Lc <- tryCatch(t(chol(S)), error = function(e) diag(sqrt(0.5), d))
        par <- g$beta
        for (k in seq_len(d)) for (l in seq_len(k))
          par <- c(par, if (k == l) log(max(Lc[k, k], 1e-3)) else Lc[k, l])
        start <- par
      }
    }
    if (is.null(start)) {
      par <- rep(0, p)
      for (k in seq_len(d)) for (l in seq_len(k))
        par <- c(par, if (k == l) log(0.7) else 0)
      start <- par
    }
  }

  ## 주의: optim(par, fn, gr, ...) 에 p = p 를 넘기면 R의 부분매칭 때문에
  ## p 가 par 로 매칭돼 버린다("p"는 "par"의 접두사). 클로저로 감싸서
  ## optim 의 ... 로 아무것도 넘기지 않는다.
  nll <- function(par) aghq_negloglik(par, dat = dat, p = p, d = d,
                                      n_quad = n_quad, gh = gh)
  opt <- optim(start, nll, method = "BFGS", control = list(maxit = maxit))
  beta <- opt$par[seq_len(p)]
  L <- matrix(0, d, d); idx <- p
  for (k in seq_len(d)) for (l in seq_len(k)) {
    idx <- idx + 1
    L[k, l] <- if (k == l) exp(opt$par[idx]) else opt$par[idx]
  }
  list(beta = beta, Sigma = L %*% t(L), negll = opt$value,
       converged = opt$convergence == 0, n_eval = opt$counts[1])
}
