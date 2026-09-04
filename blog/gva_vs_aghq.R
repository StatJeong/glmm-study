###############################################################################
##  GVA vs AGHQ  --  같은 조건에서 시간과 편향 비교
##
##  실행: RStudio에서 이 파일 열고 Source (Cmd/Ctrl + Shift + S)
##        (working directory에 pilot_sim.R 이 있어야 함. mywork.Rproj 로 열면 됨)
##
##  ---------------------------------------------------------------------------
##  주의 1. pilot_sim.R 의 offset 부호 버그가 수정된 상태를 전제한다.
##          fit_gva 안에 다음 줄이 있어야 한다:  else mu - lam/2
##          (수정 전은 -mu + lam/2 였고, 이 경우 GVA 결과가 전부 틀린다.)
##          아래 sanity check 가 자동으로 확인한다.
##
##  주의 2. 이전 체크포인트(param_sweep.rds)는 버그가 있던 코드로 만든 것이라
##          재사용하면 안 된다. 이 스크립트는 새 파일명을 쓴다.
##
##  주의 3. GVA는 "ELBO를 최대화하는 추정량"이므로 수렴할 때까지 돌려야 한다.
##          max_iter 가 작으면 수렴 전에 멈춰 가짜 편향이 잡힌다(150회로는
##          한참 모자란다). 기본값 5000 으로 두고 수렴 여부를 기록한다.
###############################################################################

## ==== 설정 ==================================================================
CFG <- list(
  ## --- 모수 그리드 ---
  beta0_grid  = c(1, 2, 3),        # 고정절편
  beta1_grid  = c(1, 2, 3),        # 고정기울기
  sigma2_grid = c(1, 2, 4),        # 랜덤절편 분산

  ## --- 표본 크기 ---
  ## m = 10000 은 매우 오래 걸린다(아래 "예상 비용" 참고).
  ## 먼저 c(100, 1000) 으로 돌려보고, 나중에 10000 을 추가해 다시 Source 하면
  ## 체크포인트 덕분에 m=10000 만 이어서 계산한다.
  m_grid      = c(100, 1000),
  n_grid      = c(5, 25, 100),

  reps        = 5,
  family      = "gamma",
  nu          = 2,                 # Gamma shape (고정 상수)
  n_quad      = 11,                # AGHQ 노드 수

  ## --- GVA 수렴 설정 ---
  gva_max_iter = 5000,
  gva_tol      = 1e-9,

  n_cores     = max(1, parallel::detectCores() - 2),
  base_seed   = 1,
  checkpoint  = "gva_aghq.rds",
  csv_out     = "gva_aghq_summary.csv"
)

## ==== 준비 ==================================================================
if (!file.exists("pilot_sim.R"))
  stop("pilot_sim.R 를 찾을 수 없습니다. working directory 를 확인하세요: ", getwd())
source("pilot_sim.R")
if (!requireNamespace("statmod", quietly = TRUE))
  stop("statmod 패키지가 필요합니다:  install.packages('statmod')")

## --- sanity check: offset 부호 버그가 고쳐져 있는지 ---
if (grepl("-mu + lam/2", paste(deparse(fit_gva), collapse = " "), fixed = TRUE))
  stop("pilot_sim.R 의 fit_gva 에 offset 부호 버그가 남아 있습니다.\n",
       "  'else -mu + lam/2'  ->  'else mu - lam/2'  로 고쳐야 합니다.")

PARAM_GRID <- expand.grid(beta0  = CFG$beta0_grid,
                          beta1  = CFG$beta1_grid,
                          sigma2 = CFG$sigma2_grid,
                          KEEP.OUT.ATTRS = FALSE)

## 인덱스 기반 시드 -- (조합, m, n, rep)이 다르면 시드도 반드시 다르다
make_seed <- function(base_seed, pi_, mi, ni, rr)
  base_seed * 1e7 + pi_ * 1e5 + mi * 1e4 + ni * 1e3 + rr

## ==== 한 건: 같은 데이터에 GVA와 AGHQ를 둘 다 적합 ==========================
run_one_task <- function(tk, param_grid, cfg, gh) {
  b0 <- param_grid$beta0[tk$pi_]
  b1 <- param_grid$beta1[tk$pi_]
  s2 <- param_grid$sigma2[tk$pi_]

  set.seed(make_seed(cfg$base_seed, tk$pi_, tk$mi, tk$ni, tk$rep))
  d <- simulate_data(tk$m, tk$n, b0, b1, s2, cfg$family, cfg$nu)
  ctx <- sprintf("b0=%g b1=%g s2=%g m=%d n=%d rep=%d", b0, b1, s2, tk$m, tk$n, tk$rep)

  ## --- GVA (ELBO 최대화) ---
  t0 <- Sys.time()
  g <- .fit_or_na(fit_gva(d$X, d$Y, cfg$family, cfg$nu,
                          max_iter = cfg$gva_max_iter, tol = cfg$gva_tol),
                  4, "GVA", ctx)
  gva_time <- as.numeric(Sys.time() - t0, units = "secs")

  ## --- AGHQ (adaptive Gauss-Hermite quadrature) ---
  t0 <- Sys.time()
  a <- .fit_or_na(fit_exact(d$X, d$Y, cfg$family, cfg$nu, cfg$n_quad, gh = gh),
                  3, "AGHQ", ctx)
  aghq_time <- as.numeric(Sys.time() - t0, units = "secs")

  ## GVA가 실제로 수렴했는지 (반복 상한에 걸렸으면 FALSE)
  gva_conv <- !is.na(g["iters"]) && g["iters"] < cfg$gva_max_iter

  ## 도달한 ELBO (수렴 품질 확인용)
  elbo <- if (any(is.na(g[1:3]))) NA_real_ else
    tryCatch(full_elbo(g["beta0"], g["beta1"], g["sigma2"],
                       d$X, d$Y, cfg$family, cfg$nu), error = function(e) NA_real_)

  data.frame(
    true_beta0 = b0, true_beta1 = b1, true_sigma2 = s2,
    m = tk$m, n = tk$n, rep = tk$rep,
    gva_beta0  = unname(g["beta0"]), gva_beta1  = unname(g["beta1"]),
    gva_sigma2 = unname(g["sigma2"]),
    aghq_beta0 = unname(a["beta0"]), aghq_beta1 = unname(a["beta1"]),
    aghq_sigma2 = unname(a["sigma2"]),
    gva_time = gva_time, aghq_time = aghq_time,
    gva_iters = unname(g["iters"]), gva_converged = unname(gva_conv),
    gva_elbo = unname(elbo),
    gva_err = .errmsg_of(g), aghq_err = .errmsg_of(a),
    key = tk$key, stringsAsFactors = FALSE)
}

## ==== 메인 ==================================================================
run_comparison <- function(cfg = CFG, param_grid = PARAM_GRID, batch_size = 40) {

  gh <- statmod::gauss.quad(cfg$n_quad, kind = "hermite")

  ## (모수조합 x m x n x rep) 을 전부 평탄화해서 코어에 흩뿌린다.
  ## reps 가 코어 수보다 작아도 코어가 놀지 않게 하려면 이렇게 해야 한다.
  tasks <- expand.grid(rep = seq_len(cfg$reps),
                       ni  = seq_along(cfg$n_grid),
                       mi  = seq_along(cfg$m_grid),
                       pi_ = seq_len(nrow(param_grid)),
                       KEEP.OUT.ATTRS = FALSE)
  tasks$m   <- cfg$m_grid[tasks$mi]
  tasks$n   <- cfg$n_grid[tasks$ni]
  tasks$obs <- tasks$m * tasks$n
  tasks$key <- sprintf("p%d_m%d_n%d_r%d", tasks$pi_, tasks$m, tasks$n, tasks$rep)

  results <- list(); done <- character(0)
  if (file.exists(cfg$checkpoint)) {
    prev <- readRDS(cfg$checkpoint)
    results <- prev$results; done <- prev$done
    message(sprintf("체크포인트: %d개 완료됨, 이어서 진행합니다.", length(done)))
  }
  tasks <- tasks[!(tasks$key %in% done), ]
  if (!nrow(tasks)) { message("모든 작업 완료됨."); return(do.call(rbind, results)) }

  tasks <- tasks[order(tasks$obs), ]     # 작은 칸부터 -> 중간 결과를 일찍 봄
  n_all <- nrow(tasks); n_ok <- 0L; t_all <- Sys.time()
  message(sprintf("남은 작업 %d개 | 코어 %d개 | reps=%d | GVA max_iter=%d",
                  n_all, cfg$n_cores, cfg$reps, cfg$gva_max_iter))

  for (ob in unique(tasks$obs)) {
    idx <- which(tasks$obs == ob)
    ## 워커당 m*n 크기 행렬을 여러 개 잡으므로 큰 칸은 코어를 줄인다
    cores <- if (ob > 5e6) max(1L, min(cfg$n_cores, 3L))
             else if (ob > 5e5) max(1L, min(cfg$n_cores, 6L))
             else cfg$n_cores
    message(sprintf("\n--- m=%d n=%d (%.0e obs) : %d개, 코어 %d ---",
                    tasks$m[idx[1]], tasks$n[idx[1]], ob, length(idx), cores))

    ## RStudio(macOS)에서 fork(mclapply)는 불안정하므로 PSOCK 을 쓴다
    cl <- NULL
    if (cores > 1) {
      cl <- parallel::makeCluster(cores)
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
      parallel::clusterEvalQ(cl, { suppressMessages(source("pilot_sim.R")); NULL })
      parallel::clusterExport(cl, c("run_one_task", "make_seed"), envir = environment())
    }

    for (st in seq(1, length(idx), by = batch_size)) {
      bidx <- idx[st:min(st + batch_size - 1, length(idx))]
      tk_list <- split(tasks[bidx, ], seq_along(bidx))
      t0 <- Sys.time()

      res <- if (!is.null(cl))
        parallel::parLapplyLB(cl, tk_list, run_one_task,
                              param_grid = param_grid, cfg = cfg, gh = gh)
      else lapply(tk_list, run_one_task, param_grid = param_grid, cfg = cfg, gh = gh)

      bad <- !vapply(res, is.data.frame, logical(1))
      if (any(bad)) {
        for (bi in which(bad))
          warning(sprintf("작업 실패 %s: %s", tasks$key[bidx[bi]],
                          paste(as.character(res[[bi]]), collapse = " ")), call. = FALSE)
        res <- res[!bad]
      }
      if (!length(res)) next

      chunk <- do.call(rbind, res)
      results[[length(results) + 1L]] <- chunk
      done <- c(done, chunk$key); n_ok <- n_ok + nrow(chunk)
      saveRDS(list(results = results, done = done,
                   param_grid = param_grid, cfg = cfg), cfg$checkpoint)

      message(sprintf("  [%4d/%4d] %2d건 | %5.1f초 | GVA %6.2f s/fit (수렴 %d/%d), AGHQ %5.2f s/fit | 누적 %.1f분",
                      n_ok, n_all, nrow(chunk),
                      as.numeric(Sys.time() - t0, units = "secs"),
                      mean(chunk$gva_time, na.rm = TRUE),
                      sum(chunk$gva_converged, na.rm = TRUE), nrow(chunk),
                      mean(chunk$aghq_time, na.rm = TRUE),
                      as.numeric(Sys.time() - t_all, units = "mins")))
    }
    if (!is.null(cl)) { try(parallel::stopCluster(cl), silent = TRUE); cl <- NULL }
    gc()
  }

  out <- do.call(rbind, results)
  message(sprintf("\n완료: %d행, 총 %.1f분", nrow(out),
                  as.numeric(Sys.time() - t_all, units = "mins")))
  out
}

## ==== 요약표 ================================================================
## 조합 x (m,n) 별로 두 방법의 편향과 시간을 나란히 놓는다.
summarize_comparison <- function(df) {
  key <- list(df$true_beta0, df$true_beta1, df$true_sigma2, df$m, df$n)
  out <- do.call(rbind, lapply(split(df, key, drop = TRUE), function(z) {
    mb <- function(v, tv) mean(z[[v]] - z[[tv]], na.rm = TRUE)   # 편향
    data.frame(
      beta0 = z$true_beta0[1], beta1 = z$true_beta1[1], sigma2 = z$true_sigma2[1],
      m = z$m[1], n = z$n[1], reps = nrow(z),
      ## --- 편향 (추정 평균 - 참값) ---
      gva_b0_bias  = mb("gva_beta0","true_beta0"),
      aghq_b0_bias = mb("aghq_beta0","true_beta0"),
      gva_b1_bias  = mb("gva_beta1","true_beta1"),
      aghq_b1_bias = mb("aghq_beta1","true_beta1"),
      gva_s2_bias  = mb("gva_sigma2","true_sigma2"),
      aghq_s2_bias = mb("aghq_sigma2","true_sigma2"),
      ## --- 시간 ---
      gva_time  = mean(z$gva_time,  na.rm = TRUE),
      aghq_time = mean(z$aghq_time, na.rm = TRUE),
      time_ratio = mean(z$gva_time, na.rm = TRUE) / mean(z$aghq_time, na.rm = TRUE),
      ## --- 품질 ---
      gva_iters = mean(z$gva_iters, na.rm = TRUE),
      gva_conv_rate = mean(z$gva_converged, na.rm = TRUE),
      n_fail = sum(is.na(z$gva_beta0)) + sum(is.na(z$aghq_beta0)),
      row.names = NULL)
  }))
  out[order(out$sigma2, out$beta0, out$beta1, out$n, out$m), ]
}

## m,n 별로만 뭉뚱그린 요약 (모수조합에 걸쳐 평균) -- 발표용 핵심 표
summarize_by_size <- function(df) {
  do.call(rbind, lapply(split(df, list(df$m, df$n), drop = TRUE), function(z) data.frame(
    m = z$m[1], n = z$n[1], reps = nrow(z),
    gva_b0_bias  = mean(z$gva_beta0  - z$true_beta0,  na.rm = TRUE),
    aghq_b0_bias = mean(z$aghq_beta0 - z$true_beta0,  na.rm = TRUE),
    gva_s2_bias  = mean(z$gva_sigma2 - z$true_sigma2, na.rm = TRUE),
    aghq_s2_bias = mean(z$aghq_sigma2- z$true_sigma2, na.rm = TRUE),
    gva_time  = mean(z$gva_time,  na.rm = TRUE),
    aghq_time = mean(z$aghq_time, na.rm = TRUE),
    time_ratio = mean(z$gva_time, na.rm = TRUE)/mean(z$aghq_time, na.rm = TRUE),
    gva_conv_rate = mean(z$gva_converged, na.rm = TRUE),
    row.names = NULL)))
}

load_results <- function(file = CFG$checkpoint) do.call(rbind, readRDS(file)$results)

## ==== 실행 ==================================================================
raw <- run_comparison()
sum_full <- summarize_comparison(raw)
sum_size <- summarize_by_size(raw)

write.csv(sum_full, CFG$csv_out, row.names = FALSE)
message(sprintf("요약 저장: %s (%d행)", CFG$csv_out, nrow(sum_full)))

cat("\n=========== m, n 별 요약 (모수조합 평균) ===========\n")
print(sum_size, row.names = FALSE, digits = 3)

cat("\n=========== GVA 수렴 실패한 칸 (있으면 max_iter 를 늘려야 함) ===========\n")
bad <- sum_full[sum_full$gva_conv_rate < 1, ]
if (nrow(bad)) {
  print(head(bad[, c("beta0","beta1","sigma2","m","n","gva_iters","gva_conv_rate")], 20),
        row.names = FALSE)
} else {
  cat("없음 -- 모든 적합이 수렴했습니다.\n")
}
