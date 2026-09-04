###############################################################################
##  run_dimension_study.R
##  랜덤효과 차원 d 를 늘려가며 GVA vs AGHQ 의 편향과 비용을 비교한다.
##
##  실행: RStudio에서 Source (Cmd/Ctrl + Shift + S)
##        working directory 에 glmm_multi.R 이 있어야 한다.
##
##  ---------------------------------------------------------------------------
##  왜 차원인가
##    AGHQ 비용은 그룹당 q^d (q = 노드 수, d = 랜덤효과 차원) 으로 지수 증가한다.
##    GVA 는 변분모수가 d + d(d+1)/2 개로 다항 증가한다.
##    d = 1 (랜덤절편만) 은 하필 AGHQ 가 가장 유리한 조건이라 두 방법이 갈리지
##    않는다. 이 스크립트는 d 를 늘렸을 때 어디서 갈리는지를 측정한다.
##
##  모형
##    eta_ij = beta0 + sum_k beta_k X_kij + U_i0 + sum_k U_ik X_kij
##    U_i ~ N(0, Sigma),  dim = d
##    Y_ij | U_i ~ Gamma(shape = nu, rate = nu/exp(eta_ij))
###############################################################################

## ==== 설정 ==================================================================
CFG <- list(
  d_grid   = c(1, 2, 3),      # 랜덤효과 차원. 4 이상은 AGHQ가 급격히 느려진다.
  m_grid   = c(200, 500),     # 그룹 수
  n_grid   = c(20, 40),       # 그룹당 관측치
  reps     = 10,

  ## 아래 beta0/beta_slope 기본값은 로지스틱 기준이다(0 근처가 정보량이 많다).
  ## family 를 gamma/poisson 으로 바꿀 때는 set_family() 를 쓰면 함께 조정된다.
  beta0    = 0.5,             # 고정절편
  beta_slope = 0.8,           # 고정기울기 (모든 공변량 공통)
  sig_var  = 1.0,             # Sigma 대각 (분산)
  sig_cor  = 0.3,             # Sigma 비대각 상관

  ## family: "binomial"(로지스틱) / "gamma" / "poisson"
  ## 로지스틱은 ELBO가 닫힌 형태가 아니어서 GVA 안에서 관측치별 1차원 GH 구적을
  ## 쓴다(비용 O(q), 차원 무관). AGHQ 의 q^d 와 대비되는 지점.
  family   = "gamma",
  trials   = 1,               # 1이면 베르누이. 크게 하면 정보량이 늘어난다.
  nu       = 2,               # gamma 일 때만 사용

  n_quad   = 7,               # AGHQ 노드 수. d가 커지면 q^d 라 작게 잡아야 한다.
  n_quad_gva = 11,            # GVA 내부 1차원 구적 노드 수 (binomial 전용)

  ## GVA 최적화 방식
  ##  "squarem": SQUAREM 가속 (기본). beta0 와 m_i 의 aliasing 때문에 생기는
  ##             선형 수렴을 외삽으로 건너뛴다. 감마에서 G호출이 700회 -> 20~35회.
  ##  "none"   : 가속 없는 단순 좌표상승. 비교용.
  accel        = "squarem",
  gva_max_iter = 3000,        # accel="squarem" 이면 외곽 반복 수. 500이면 충분.
  gva_tol      = 1e-8,

  n_cores    = max(1, parallel::detectCores() - 2),
  base_seed  = 1,

  ## 저장 파일명은 family 에 따라 자동으로 붙는다(아래 참조).
  ## 직접 지정하고 싶으면 여기에 문자열을 넣으면 그게 우선한다.
  checkpoint = NULL,
  csv_out    = NULL,
  tag        = "squarem"             # 같은 family 로 여러 설정을 돌릴 때 구분용 꼬리표
)

## --- 저장 파일명 자동 생성 -------------------------------------------------
## family 를 바꿔 돌릴 때 결과가 섞이지 않도록 파일명에 family 를 박는다.
##   binomial -> dim_study_binomial.rds / dim_study_binomial_summary.csv
auto_names <- function(cfg) {
  sfx <- paste0(cfg$family, if (nzchar(cfg$tag)) paste0("_", cfg$tag) else "")
  if (is.null(cfg$checkpoint)) cfg$checkpoint <- sprintf("dim_study_%s.rds", sfx)
  if (is.null(cfg$csv_out))    cfg$csv_out    <- sprintf("dim_study_%s_summary.csv", sfx)
  cfg
}
CFG <- auto_names(CFG)

## --- family 전환 헬퍼 -------------------------------------------------------
## family 마다 정보량이 크게 달라서 같은 (beta, n) 이 적절하지 않다.
##   로지스틱: 베르누이라 관측치당 정보가 적고, eta 가 0에서 멀면 포화된다
##   gamma/poisson: eta 가 로그평균이라 절편을 1 정도로 두는 게 자연스럽다
## 파일명도 함께 다시 계산한다.
##
## 사용:  CFG <- set_family(CFG, "gamma")   그리고 Source
set_family <- function(cfg, family = c("binomial", "gamma", "poisson"),
                       n_grid = NULL) {
  family <- match.arg(family)
  cfg$family <- family
  if (family == "binomial") {
    cfg$beta0 <- 0.5; cfg$beta_slope <- 0.8
    if (is.null(n_grid)) n_grid <- c(20, 40)      # 이항은 n 이 커야 한다
  } else {
    cfg$beta0 <- 1.0; cfg$beta_slope <- 0.5
    if (is.null(n_grid)) n_grid <- c(10, 30)
  }
  cfg$n_grid <- n_grid
  cfg$checkpoint <- NULL; cfg$csv_out <- NULL     # 파일명 다시 생성
  auto_names(cfg)
}

## --- 설정 지문 -------------------------------------------------------------
## 파일명만으로는 beta0 나 n_quad 같은 걸 바꿨을 때를 못 잡는다.
## 결과에 영향을 주는 설정을 지문으로 남겨두고, 이어받기 전에 대조한다.
config_fingerprint <- function(cfg) {
  keys <- c("family", "trials", "nu", "beta0", "beta_slope", "sig_var", "sig_cor",
            "n_quad", "n_quad_gva", "base_seed", "gva_max_iter", "gva_tol", "accel")
  paste(vapply(keys, function(k) paste0(k, "=", paste(cfg[[k]], collapse = ",")),
               character(1)), collapse = "; ")
}

## ==== 준비 ==================================================================
if (!file.exists("glmm_multi.R"))
  stop("glmm_multi.R 를 찾을 수 없습니다. working directory: ", getwd())
source("glmm_multi.R")

## 차원 d 의 참 Sigma: 대각 sig_var, 비대각 sig_cor*sig_var
make_Sigma <- function(d, cfg) {
  S <- matrix(cfg$sig_cor * cfg$sig_var, d, d)
  diag(S) <- cfg$sig_var
  S
}
make_beta <- function(d, cfg) c(cfg$beta0, rep(cfg$beta_slope, max(d - 1, 1)))

make_seed <- function(base, di, mi, ni, rr)
  base * 1e7 + di * 1e5 + mi * 1e4 + ni * 1e3 + rr

## ==== 한 건 ================================================================
run_one <- function(tk, cfg) {
  d <- tk$d
  Sigma <- make_Sigma(d, cfg); beta <- make_beta(d, cfg)
  p <- length(beta)

  set.seed(make_seed(cfg$base_seed, tk$di, tk$mi, tk$ni, tk$rep))
  dat <- simulate_data_multi(tk$m, tk$n, beta, Sigma, nu = cfg$nu,
                             family = cfg$family,
                             trials = if (is.null(cfg$trials)) 1 else cfg$trials)

  ## --- GVA ---
  t0 <- Sys.time()
  g <- tryCatch(fit_gva_multi(dat, p = p, max_iter = cfg$gva_max_iter,
                              tol = cfg$gva_tol,
                              n_quad_gva = if (is.null(cfg$n_quad_gva)) 11
                                           else cfg$n_quad_gva,
                              accel = if (is.null(cfg$accel)) "squarem"
                                      else cfg$accel),
                error = function(e) NULL)
  gva_time <- as.numeric(Sys.time() - t0, units = "secs")

  ## --- AGHQ (시작값은 GVA 결과 -- c(0,...,0) 시작은 기울기가 크면 발산한다) ---
  t0 <- Sys.time()
  a <- tryCatch(fit_aghq_multi(dat, p = p, n_quad = cfg$n_quad,
                               gva_start = TRUE),
                error = function(e) NULL)
  aghq_time <- as.numeric(Sys.time() - t0, units = "secs")

  ## 요약 지표: beta0, 첫 기울기, Sigma 대각 평균
  gv <- if (is.null(g)) rep(NA_real_, 3) else
    c(g$beta[1], if (p > 1) g$beta[2] else NA_real_, mean(diag(g$Sigma)))
  av <- if (is.null(a)) rep(NA_real_, 3) else
    c(a$beta[1], if (p > 1) a$beta[2] else NA_real_, mean(diag(a$Sigma)))

  data.frame(
    d = d, m = tk$m, n = tk$n, rep = tk$rep,
    true_beta0 = beta[1],
    true_slope = if (p > 1) beta[2] else NA_real_,
    true_sigvar = cfg$sig_var,
    gva_beta0 = gv[1], gva_slope = gv[2], gva_sigvar = gv[3],
    aghq_beta0 = av[1], aghq_slope = av[2], aghq_sigvar = av[3],
    gva_time = gva_time, aghq_time = aghq_time,
    gva_iters = if (is.null(g)) NA_real_ else g$iters,
    gva_nG    = if (is.null(g) || is.null(g$n_G)) NA_real_ else g$n_G,
    gva_converged = if (is.null(g)) FALSE else g$converged,
    aghq_converged = if (is.null(a)) FALSE else a$converged,
    aghq_nodes = cfg$n_quad^d,
    key = tk$key, stringsAsFactors = FALSE)
}

## ==== 메인 ==================================================================
run_study <- function(cfg = CFG, batch_size = 20) {
  tasks <- expand.grid(rep = seq_len(cfg$reps),
                       ni  = seq_along(cfg$n_grid),
                       mi  = seq_along(cfg$m_grid),
                       di  = seq_along(cfg$d_grid),
                       KEEP.OUT.ATTRS = FALSE)
  tasks$d <- cfg$d_grid[tasks$di]
  tasks$m <- cfg$m_grid[tasks$mi]
  tasks$n <- cfg$n_grid[tasks$ni]
  tasks$cost <- tasks$m * tasks$n * cfg$n_quad^tasks$d   # 대략적 비용 지표
  tasks$key <- sprintf("d%d_m%d_n%d_r%d", tasks$d, tasks$m, tasks$n, tasks$rep)

  fp <- config_fingerprint(cfg)
  results <- list(); done <- character(0)
  if (file.exists(cfg$checkpoint)) {
    prev <- readRDS(cfg$checkpoint)
    ## 설정이 바뀌었는데 이어받으면 다른 조건의 결과가 한 파일에 섞인다.
    ## 예전 파일(지문 없음)도 여기서 걸러진다.
    if (is.null(prev$fingerprint) || !identical(prev$fingerprint, fp)) {
      stop("체크포인트 '", cfg$checkpoint, "' 는 다른 설정으로 만들어졌습니다.\n",
           "  파일: ", if (is.null(prev$fingerprint)) "(지문 없음 - 구버전)"
                        else prev$fingerprint, "\n",
           "  현재: ", fp, "\n",
           "파일을 지우거나 CFG$tag 를 다르게 주고 다시 실행하세요.")
    }
    results <- prev$results; done <- prev$done
    message(sprintf("체크포인트 '%s': %d개 완료됨, 이어서 진행합니다.",
                    cfg$checkpoint, length(done)))
  }
  tasks <- tasks[!(tasks$key %in% done), ]
  if (!nrow(tasks)) { message("모든 작업 완료됨."); return(do.call(rbind, results)) }

  tasks <- tasks[order(tasks$cost), ]     # 싼 것부터
  n_all <- nrow(tasks); n_ok <- 0L; t_all <- Sys.time()
  message(sprintf("family=%s | 남은 작업 %d개 | 코어 %d | reps=%d | d=%s",
                  cfg$family, n_all, cfg$n_cores, cfg$reps,
                  paste(cfg$d_grid, collapse = ",")))
  message(sprintf("저장: %s , %s | GVA 최적화: %s",
                  cfg$checkpoint, cfg$csv_out,
                  if (is.null(cfg$accel)) "squarem" else cfg$accel))

  for (ck in unique(tasks$cost)) {
    idx <- which(tasks$cost == ck)
    tk1 <- tasks[idx[1], ]
    cores <- if (ck > 5e7) max(1L, min(cfg$n_cores, 2L))
             else if (ck > 5e6) max(1L, min(cfg$n_cores, 4L))
             else cfg$n_cores
    message(sprintf("\n--- d=%d m=%d n=%d (AGHQ 노드 %d) : %d개, 코어 %d ---",
                    tk1$d, tk1$m, tk1$n, cfg$n_quad^tk1$d, length(idx), cores))

    cl <- NULL
    if (cores > 1) {
      cl <- parallel::makeCluster(cores)
      parallel::clusterEvalQ(cl, { source("glmm_multi.R"); NULL })
      parallel::clusterExport(cl, c("run_one", "make_seed", "make_Sigma", "make_beta"),
                              envir = environment())
    }

    for (st in seq(1, length(idx), by = batch_size)) {
      bidx <- idx[st:min(st + batch_size - 1, length(idx))]
      tk_list <- split(tasks[bidx, ], seq_along(bidx))
      t0 <- Sys.time()
      res <- if (!is.null(cl)) parallel::parLapplyLB(cl, tk_list, run_one, cfg = cfg)
             else lapply(tk_list, run_one, cfg = cfg)

      bad <- !vapply(res, is.data.frame, logical(1))
      if (any(bad)) { for (bi in which(bad))
        warning(sprintf("실패 %s", tasks$key[bidx[bi]]), call. = FALSE)
        res <- res[!bad] }
      if (!length(res)) next

      chunk <- do.call(rbind, res)
      results[[length(results) + 1L]] <- chunk
      done <- c(done, chunk$key); n_ok <- n_ok + nrow(chunk)
      saveRDS(list(results = results, done = done, cfg = cfg, fingerprint = fp),
              cfg$checkpoint)
      message(sprintf("  [%3d/%3d] %2d건 | %5.1f초 | GVA %6.2f s/fit, AGHQ %6.2f s/fit | 누적 %.1f분",
                      n_ok, n_all, nrow(chunk),
                      as.numeric(Sys.time() - t0, units = "secs"),
                      mean(chunk$gva_time), mean(chunk$aghq_time),
                      as.numeric(Sys.time() - t_all, units = "mins")))
    }
    if (!is.null(cl)) parallel::stopCluster(cl)
    gc()
  }

  out <- do.call(rbind, results)
  message(sprintf("\n완료: %d행, 총 %.1f분", nrow(out),
                  as.numeric(Sys.time() - t_all, units = "mins")))
  out
}

## ==== 요약 ==================================================================
summarize_study <- function(df) {
  ok <- df[df$gva_converged & df$aghq_converged, ]
  do.call(rbind, lapply(split(ok, list(ok$d, ok$m, ok$n), drop = TRUE), function(z)
    data.frame(
      d = z$d[1], m = z$m[1], n = z$n[1], reps = nrow(z),
      aghq_nodes = z$aghq_nodes[1],
      gva_b0_bias   = mean(z$gva_beta0  - z$true_beta0),
      aghq_b0_bias  = mean(z$aghq_beta0 - z$true_beta0),
      gva_sl_bias   = mean(z$gva_slope  - z$true_slope),
      aghq_sl_bias  = mean(z$aghq_slope - z$true_slope),
      gva_sv_bias   = mean(z$gva_sigvar - z$true_sigvar),
      aghq_sv_bias  = mean(z$aghq_sigvar- z$true_sigvar),
      gva_time  = mean(z$gva_time),
      aghq_time = mean(z$aghq_time),
      speedup   = mean(z$aghq_time) / mean(z$gva_time),   # >1 이면 GVA가 빠름
      gva_iters = mean(z$gva_iters),
      gva_nG    = mean(z$gva_nG),
      row.names = NULL)))
}

## 차원별 비용 곡선 (m, n 에 걸쳐 평균)
cost_curve <- function(df) {
  ok <- df[df$gva_converged & df$aghq_converged, ]
  do.call(rbind, lapply(split(ok, ok$d), function(z) data.frame(
    d = z$d[1], aghq_nodes = z$aghq_nodes[1], reps = nrow(z),
    gva_time = mean(z$gva_time), aghq_time = mean(z$aghq_time),
    speedup = mean(z$aghq_time) / mean(z$gva_time),
    gva_nG = mean(z$gva_nG),
    max_abs_bias_diff = max(abs(c(
      mean(z$gva_beta0 - z$true_beta0)  - mean(z$aghq_beta0 - z$true_beta0),
      mean(z$gva_sigvar- z$true_sigvar) - mean(z$aghq_sigvar- z$true_sigvar)))),
    row.names = NULL)))
}

load_study <- function(f = CFG$checkpoint) do.call(rbind, readRDS(f)$results)

## 여러 family 결과를 한 번에 읽어 비교하고 싶을 때
load_all_families <- function(pattern = "^dim_study_.*\\.rds$") {
  fs <- list.files(".", pattern = pattern)
  do.call(rbind, lapply(fs, function(f) {
    ck <- readRDS(f); z <- do.call(rbind, ck$results)
    z$family <- ck$cfg$family; z
  }))
}

## ==== 실행 ==================================================================
raw <- run_study()
smry <- summarize_study(raw)
write.csv(smry, CFG$csv_out, row.names = FALSE)
message(sprintf("요약 저장: %s", CFG$csv_out))

cat("\n=========== d, m, n 별 편향과 시간 ===========\n")
print(smry, row.names = FALSE, digits = 3)

cat("\n=========== 차원별 비용 곡선 (핵심) ===========\n")
cat("speedup > 1 이면 GVA가 빠름. d가 커질수록 커지는지가 관건.\n\n")
print(cost_curve(raw), row.names = FALSE, digits = 3)

cat("\n=========== 수렴 실패 ===========\n")
f1 <- sum(!raw$gva_converged); f2 <- sum(!raw$aghq_converged)
cat(sprintf("GVA 미수렴 %d/%d, AGHQ 미수렴 %d/%d\n", f1, nrow(raw), f2, nrow(raw)))
if (f1 > 0) print(table(d = raw$d[!raw$gva_converged], n = raw$n[!raw$gva_converged]))
