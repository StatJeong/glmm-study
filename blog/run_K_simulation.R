###############################################################################
##  run_K_simulation.R
##  랜덤효과 차원 K = 1, 2, 3 에서 GVA 와 AGHQ 의 편향·시간을 비교한다.
##  (Poisson, log link -- epilepsy 와 같은 반응 유형)
##
##  ---------------------------------------------------------------------------
##  두 가지 Z 설계를 나란히 본다. 난이도가 다르다.
##
##   design = "poly"  (종단형, 주 관심)
##       z_ij = (1, t_j, t_j^2)      t 는 중심화된 시간
##       1, t, t^2 가 강하게 상관되어 Z 가 악조건이다.
##       epilepsy 처럼 방문 시점이 고정된 종단자료에서 실제로 마주치는 구조.
##
##   design = "indep" (대조군)
##       z_ij = (1, X1_ij, X2_ij)    X 는 독립 N(0,1)
##       Z 열이 무상관이라 조건이 좋다. 앞선 d-스터디가 쓰던 설계.
##
##  참값을 알고 있으므로 여기서는 **편향**을 잴 수 있다.
##  (실데이터에서는 참값을 모르므로 근사오차 |GVA - MLE| 만 측정 가능했다.)
###############################################################################

source("glmm_multi.R")

CFG <- list(
  K_grid  = c(1, 2, 3),
  m_grid  = c(200, 500),
  n       = 8,                 # 그룹당 시점 수 (K=3 이면 최소 3 이상 필요)
  reps    = 30,

  beta    = c(0.5, 0.3, -0.2), # 고정효과: (절편, t, t^2) 또는 (절편, X1, X2)
  sig     = c(0.50, 0.25, 0.15),  # 랜덤효과 표준편차 (K 개까지 사용)
  rho     = 0.2,               # 랜덤효과 간 상관

  n_quad  = 9,                 # AGHQ 노드. K=3 이면 9^3 = 729
  gva_tol = 1e-9, gva_max_iter = 1000,

  n_cores = max(1, parallel::detectCores() - 2),
  base_seed = 500,
  checkpoint = "Ksim.rds"
)

## 참 Sigma: 대각 sig^2, 상관 rho
make_Sigma <- function(K, cfg) {
  s <- cfg$sig[seq_len(K)]
  S <- outer(s, s) * cfg$rho
  diag(S) <- s^2
  S
}

## ---------------------------------------------------------------------------
## 데이터 생성. Xl 의 순서가 곧 Z 와 X 의 열 순서가 된다:
##   Z = (1, Xl[[1]], ..., Xl[[K-1]]),  X = (1, Xl[[1]], ..., Xl[[p-1]])
## poly  : Xl = (t, t^2)
## indep : Xl = (X1, X2)
## ---------------------------------------------------------------------------
sim_K <- function(m, n, K, cfg, design = c("poly", "indep")) {
  design <- match.arg(design)
  p <- K                                  # 고정효과 개수 = 랜덤효과 개수로 맞춘다
  beta <- cfg$beta[seq_len(p)]
  Sigma <- make_Sigma(K, cfg)

  if (design == "poly") {
    tt <- seq(-1, 1, length.out = n)      # 중심화·척도화된 시간
    Xl <- list(matrix(rep(tt, each = m), m, n),
               matrix(rep(tt^2, each = m), m, n))
  } else {
    Xl <- list(matrix(rnorm(m * n), m, n), matrix(rnorm(m * n), m, n))
  }
  Xl <- Xl[seq_len(max(K - 1, 1))]
  if (K == 1) Xl <- list(matrix(0, m, n))  # 자리만 채움 (p=1 이라 안 쓰임)

  eta <- matrix(beta[1], m, n)
  if (p > 1) for (k in seq_len(p - 1)) eta <- eta + beta[k + 1] * Xl[[k]]

  U <- matrix(rnorm(m * K), m, K) %*% chol(Sigma)
  eta <- eta + U[, 1]
  if (K > 1) for (k in seq_len(K - 1)) eta <- eta + U[, k + 1] * Xl[[k]]

  eta <- pmin(pmax(eta, -ETA_CLIP), ETA_CLIP)
  list(Y = matrix(rpois(m * n, exp(eta)), m, n), Xl = Xl, d = K,
       family = "poisson", nu = NULL, trials = 1,
       true_beta = beta, true_Sigma = Sigma, p = p)
}

## ---------------------------------------------------------------------------
run_one <- function(tk, cfg) {
  K <- tk$K
  set.seed(cfg$base_seed + tk$di * 1e5 + tk$mi * 1e4 + tk$rep +
             (tk$design == "indep") * 1e6)
  dat <- sim_K(tk$m, cfg$n, K, cfg, design = tk$design)
  p <- dat$p

  t0 <- Sys.time()
  g <- tryCatch(fit_gva_multi(dat, p = p, accel = "squarem",
                              max_iter = cfg$gva_max_iter, tol = cfg$gva_tol),
                error = function(e) NULL)
  tg <- as.numeric(Sys.time() - t0, units = "secs")

  t0 <- Sys.time()
  a <- tryCatch(fit_aghq_multi(dat, p = p, n_quad = cfg$n_quad, gva_start = TRUE),
                error = function(e) NULL)
  ta <- as.numeric(Sys.time() - t0, units = "secs")

  ok <- !is.null(g) && !is.null(a)
  data.frame(
    design = tk$design, K = K, m = tk$m, rep = tk$rep,
    true_b0 = dat$true_beta[1], true_s1 = sqrt(dat$true_Sigma[1, 1]),
    ## beta0 와 Sigma 대각 평균을 대표 지표로
    gva_b0  = if (ok) g$beta[1] else NA_real_,
    aghq_b0 = if (ok) a$beta[1] else NA_real_,
    gva_sv  = if (ok) mean(diag(g$Sigma)) else NA_real_,
    aghq_sv = if (ok) mean(diag(a$Sigma)) else NA_real_,
    true_sv = mean(diag(dat$true_Sigma)),
    gva_time = tg, aghq_time = ta,
    gva_nG = if (!is.null(g)) g$n_G else NA_real_,
    gva_conv = !is.null(g) && g$converged,
    aghq_conv = !is.null(a) && a$converged,
    aghq_nodes = cfg$n_quad^K,
    key = tk$key, stringsAsFactors = FALSE)
}

run_sim <- function(cfg = CFG, batch = 30) {
  tasks <- expand.grid(rep = seq_len(cfg$reps), mi = seq_along(cfg$m_grid),
                       di = seq_along(cfg$K_grid),
                       design = c("poly", "indep"),
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  tasks$K <- cfg$K_grid[tasks$di]; tasks$m <- cfg$m_grid[tasks$mi]
  tasks$cost <- tasks$m * cfg$n_quad^tasks$K
  tasks$key <- sprintf("%s_K%d_m%d_r%d", tasks$design, tasks$K, tasks$m, tasks$rep)

  res <- list(); done <- character(0)
  if (file.exists(cfg$checkpoint)) {
    prev <- readRDS(cfg$checkpoint); res <- prev$res; done <- prev$done
    message(sprintf("체크포인트: %d개 완료", length(done)))
  }
  tasks <- tasks[!(tasks$key %in% done), ]
  if (!nrow(tasks)) return(do.call(rbind, res))
  tasks <- tasks[order(tasks$cost), ]
  t0all <- Sys.time(); nok <- 0L
  message(sprintf("남은 %d개 | 코어 %d | reps=%d", nrow(tasks), cfg$n_cores, cfg$reps))

  for (ck in unique(tasks$cost)) {
    idx <- which(tasks$cost == ck)
    cores <- if (ck > 5e5) max(1L, min(cfg$n_cores, 4L)) else cfg$n_cores
    message(sprintf("--- %s K=%d m=%d (노드 %d) : %d개, 코어 %d ---",
                    tasks$design[idx[1]], tasks$K[idx[1]], tasks$m[idx[1]],
                    cfg$n_quad^tasks$K[idx[1]], length(idx), cores))
    cl <- NULL
    if (cores > 1) {
      cl <- parallel::makeCluster(cores)
      parallel::clusterEvalQ(cl, { source("glmm_multi.R"); NULL })
      parallel::clusterExport(cl, c("run_one","sim_K","make_Sigma"), envir = environment())
    }
    for (st in seq(1, length(idx), by = batch)) {
      b <- idx[st:min(st + batch - 1, length(idx))]
      tl <- split(tasks[b, ], seq_along(b))
      out <- if (!is.null(cl)) parallel::parLapplyLB(cl, tl, run_one, cfg = cfg)
             else lapply(tl, run_one, cfg = cfg)
      out <- out[vapply(out, is.data.frame, logical(1))]
      if (!length(out)) next
      ch <- do.call(rbind, out)
      res[[length(res)+1]] <- ch; done <- c(done, ch$key); nok <- nok + nrow(ch)
      saveRDS(list(res = res, done = done, cfg = cfg), cfg$checkpoint)
      message(sprintf("   [%d] 누적 %.1f분", nok, as.numeric(Sys.time()-t0all, units="mins")))
    }
    if (!is.null(cl)) parallel::stopCluster(cl)
    gc()
  }
  do.call(rbind, res)
}

## ==== 요약 =================================================================
summarize_K <- function(df) {
  ok <- df[df$gva_conv & df$aghq_conv, ]
  do.call(rbind, lapply(split(ok, list(ok$design, ok$K, ok$m), drop = TRUE), function(z) {
    db <- z$gva_b0 - z$aghq_b0; ds <- z$gva_sv - z$aghq_sv
    data.frame(
      design = z$design[1], K = z$K[1], m = z$m[1], reps = nrow(z),
      nodes = z$aghq_nodes[1],
      ## 참값 대비 편향 (참값을 알기에 측정 가능)
      gva_b0_bias  = mean(z$gva_b0  - z$true_b0),
      aghq_b0_bias = mean(z$aghq_b0 - z$true_b0),
      gva_sv_bias  = mean(z$gva_sv  - z$true_sv),
      aghq_sv_bias = mean(z$aghq_sv - z$true_sv),
      ## 두 방법의 대응차 (같은 데이터라 정밀)
      gap_b0 = mean(db), gap_b0_t = mean(db)/(sd(db)/sqrt(length(db))),
      gap_sv = mean(ds), gap_sv_t = mean(ds)/(sd(ds)/sqrt(length(ds))),
      ## 시간
      gva_time = mean(z$gva_time), aghq_time = mean(z$aghq_time),
      speedup = mean(z$aghq_time)/mean(z$gva_time),
      gva_nG = mean(z$gva_nG), row.names = NULL)
  }))
}

## ==== 실행 =================================================================
raw <- run_sim()
s <- summarize_K(raw)
s <- s[order(s$design, s$K, s$m), ]
write.csv(s, "Ksim_summary.csv", row.names = FALSE)

cat("\n=========== 참값 대비 편향 ===========\n")
print(s[, c("design","K","m","reps","gva_b0_bias","aghq_b0_bias",
            "gva_sv_bias","aghq_sv_bias")], row.names = FALSE, digits = 3)

cat("\n=========== GVA - AGHQ 대응차 (t 는 유의성) ===========\n")
print(s[, c("design","K","m","gap_b0","gap_b0_t","gap_sv","gap_sv_t")],
      row.names = FALSE, digits = 3)

cat("\n=========== 시간 ===========\n")
print(s[, c("design","K","m","nodes","gva_time","aghq_time","speedup","gva_nG")],
      row.names = FALSE, digits = 3)

cat("\n=========== 수렴 실패 ===========\n")
cat(sprintf("GVA %d/%d, AGHQ %d/%d\n", sum(!raw$gva_conv), nrow(raw),
            sum(!raw$aghq_conv), nrow(raw)))
