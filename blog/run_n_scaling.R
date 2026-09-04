###############################################################################
##  run_n_scaling.R
##  GVA 와 AGHQ 의 격차가 n(그룹당 관측치)에 따라 어떻게 줄어드는지 측정한다.
##
##  실행: RStudio 에서 Source. working directory 에 glmm_multi.R 이 있어야 한다.
##
##  ---------------------------------------------------------------------------
##  배경
##    앞선 실험(n = 20, 40 두 점)에서 로지스틱의 GVA-AGHQ 격차가
##    n 을 2배 하면 약 절반이 되는 패턴이 d=1,2,3 모두에서 나왔다(비율 0.42~0.49).
##    m 은 2.5배로 늘려도 격차가 거의 안 변했다(비율 1.01~1.29).
##
##    다만 n 이 두 점뿐이라 "격차 ~ n^(-alpha)" 의 alpha 를 추정할 수 없다.
##    이 스크립트는 n 을 5점으로 넓혀 지수를 추정한다.
##
##  핵심 지표
##    두 방법이 **같은 데이터**를 쓰므로 대응(paired) 차이를 본다.
##    대응 차이는 두 추정치를 따로 평균내어 빼는 것보다 분산이 훨씬 작아,
##    적은 reps 로도 미세한 체계적 차이를 잡아낸다.
###############################################################################

## ==== 설정 ==================================================================
CFG <- list(
  d_grid   = c(1, 2, 3),
  m_grid   = c(200, 500),              # m 무관함을 재확인하는 용도
  n_grid   = c(10, 20, 40, 80, 160),   # 지수 추정용 5점
  reps     = 10,

  beta0 = 0.5, beta_slope = 0.8,
  sig_var = 1.0, sig_cor = 0.3,

  family = "binomial",                 # 격차가 뚜렷한 쪽. gamma 는 값이 너무 작다
  trials = 1, nu = 2,

  n_quad = 7, n_quad_gva = 11,
  accel = "squarem", gva_max_iter = 500, gva_tol = 1e-8,

  n_cores   = max(1, parallel::detectCores() - 2),
  base_seed = 20,                      # 앞선 실험과 겹치지 않게
  tag = "nscale",
  checkpoint = NULL, csv_out = NULL
)

auto_names <- function(cfg) {
  sfx <- paste0(cfg$family, if (nzchar(cfg$tag)) paste0("_", cfg$tag) else "")
  if (is.null(cfg$checkpoint)) cfg$checkpoint <- sprintf("nscale_%s.rds", sfx)
  if (is.null(cfg$csv_out))    cfg$csv_out    <- sprintf("nscale_%s_summary.csv", sfx)
  cfg
}
CFG <- auto_names(CFG)

config_fingerprint <- function(cfg) {
  keys <- c("family","trials","nu","beta0","beta_slope","sig_var","sig_cor",
            "n_quad","n_quad_gva","base_seed","gva_max_iter","gva_tol","accel")
  paste(vapply(keys, function(k) paste0(k,"=",paste(cfg[[k]],collapse=",")),
               character(1)), collapse="; ")
}

## ==== 준비 ==================================================================
if (!file.exists("glmm_multi.R"))
  stop("glmm_multi.R 를 찾을 수 없습니다. working directory: ", getwd())
source("glmm_multi.R")

make_Sigma <- function(d, cfg) { S <- matrix(cfg$sig_cor*cfg$sig_var, d, d)
                                 diag(S) <- cfg$sig_var; S }
make_beta  <- function(d, cfg) c(cfg$beta0, rep(cfg$beta_slope, max(d-1, 1)))
make_seed  <- function(base, di, mi, ni, rr)
  base*1e7 + di*1e5 + mi*1e4 + ni*1e3 + rr

## ==== 한 건 ================================================================
run_one <- function(tk, cfg) {
  d <- tk$d; Sigma <- make_Sigma(d, cfg); beta <- make_beta(d, cfg); p <- length(beta)
  set.seed(make_seed(cfg$base_seed, tk$di, tk$mi, tk$ni, tk$rep))
  dat <- simulate_data_multi(tk$m, tk$n, beta, Sigma, nu = cfg$nu,
                             family = cfg$family, trials = cfg$trials)

  t0 <- Sys.time()
  g <- tryCatch(fit_gva_multi(dat, p = p, max_iter = cfg$gva_max_iter,
                              tol = cfg$gva_tol, n_quad_gva = cfg$n_quad_gva,
                              accel = cfg$accel), error = function(e) NULL)
  gva_time <- as.numeric(Sys.time()-t0, units="secs")

  t0 <- Sys.time()
  a <- tryCatch(fit_aghq_multi(dat, p = p, n_quad = cfg$n_quad, gva_start = TRUE),
                error = function(e) NULL)
  aghq_time <- as.numeric(Sys.time()-t0, units="secs")

  gv <- if (is.null(g)) rep(NA_real_,3) else
    c(g$beta[1], if(p>1) g$beta[2] else NA_real_, mean(diag(g$Sigma)))
  av <- if (is.null(a)) rep(NA_real_,3) else
    c(a$beta[1], if(p>1) a$beta[2] else NA_real_, mean(diag(a$Sigma)))

  data.frame(d=d, m=tk$m, n=tk$n, rep=tk$rep,
             true_beta0=beta[1], true_slope=if(p>1) beta[2] else NA_real_,
             true_sigvar=cfg$sig_var,
             gva_beta0=gv[1], gva_slope=gv[2], gva_sigvar=gv[3],
             aghq_beta0=av[1], aghq_slope=av[2], aghq_sigvar=av[3],
             gva_time=gva_time, aghq_time=aghq_time,
             gva_nG = if(is.null(g)) NA_real_ else g$n_G,
             gva_converged = if(is.null(g)) FALSE else g$converged,
             aghq_converged = if(is.null(a)) FALSE else a$converged,
             key=tk$key, stringsAsFactors=FALSE)
}

## ==== 메인 ==================================================================
run_scaling <- function(cfg = CFG, batch_size = 20) {
  tasks <- expand.grid(rep=seq_len(cfg$reps), ni=seq_along(cfg$n_grid),
                       mi=seq_along(cfg$m_grid), di=seq_along(cfg$d_grid),
                       KEEP.OUT.ATTRS=FALSE)
  tasks$d <- cfg$d_grid[tasks$di]; tasks$m <- cfg$m_grid[tasks$mi]
  tasks$n <- cfg$n_grid[tasks$ni]
  tasks$cost <- tasks$m * tasks$n * cfg$n_quad^tasks$d
  tasks$key <- sprintf("d%d_m%d_n%d_r%d", tasks$d, tasks$m, tasks$n, tasks$rep)

  fp <- config_fingerprint(cfg)
  results <- list(); done <- character(0)
  if (file.exists(cfg$checkpoint)) {
    prev <- readRDS(cfg$checkpoint)
    if (is.null(prev$fingerprint) || !identical(prev$fingerprint, fp))
      stop("체크포인트 '", cfg$checkpoint, "' 는 다른 설정으로 만들어졌습니다.\n",
           "  파일: ", prev$fingerprint, "\n  현재: ", fp,
           "\n파일을 지우거나 CFG$tag 를 다르게 주세요.")
    results <- prev$results; done <- prev$done
    message(sprintf("체크포인트: %d개 완료됨, 이어서 진행합니다.", length(done)))
  }
  tasks <- tasks[!(tasks$key %in% done), ]
  if (!nrow(tasks)) { message("모든 작업 완료됨."); return(do.call(rbind, results)) }

  tasks <- tasks[order(tasks$cost), ]
  n_all <- nrow(tasks); n_ok <- 0L; t_all <- Sys.time()
  message(sprintf("family=%s | 남은 %d개 | 코어 %d | n=%s",
                  cfg$family, n_all, cfg$n_cores, paste(cfg$n_grid, collapse=",")))
  message(sprintf("저장: %s", cfg$checkpoint))

  for (ck in unique(tasks$cost)) {
    idx <- which(tasks$cost == ck); tk1 <- tasks[idx[1], ]
    cores <- if (ck > 5e7) max(1L, min(cfg$n_cores, 2L))
             else if (ck > 5e6) max(1L, min(cfg$n_cores, 4L)) else cfg$n_cores
    message(sprintf("\n--- d=%d m=%d n=%d : %d개, 코어 %d ---",
                    tk1$d, tk1$m, tk1$n, length(idx), cores))
    cl <- NULL
    if (cores > 1) {
      cl <- parallel::makeCluster(cores)
      parallel::clusterEvalQ(cl, { source("glmm_multi.R"); NULL })
      parallel::clusterExport(cl, c("run_one","make_seed","make_Sigma","make_beta"),
                              envir = environment())
    }
    for (st in seq(1, length(idx), by = batch_size)) {
      bidx <- idx[st:min(st+batch_size-1, length(idx))]
      tk_list <- split(tasks[bidx, ], seq_along(bidx)); t0 <- Sys.time()
      res <- if (!is.null(cl)) parallel::parLapplyLB(cl, tk_list, run_one, cfg=cfg)
             else lapply(tk_list, run_one, cfg=cfg)
      res <- res[vapply(res, is.data.frame, logical(1))]
      if (!length(res)) next
      chunk <- do.call(rbind, res)
      results[[length(results)+1L]] <- chunk
      done <- c(done, chunk$key); n_ok <- n_ok + nrow(chunk)
      saveRDS(list(results=results, done=done, cfg=cfg, fingerprint=fp), cfg$checkpoint)
      message(sprintf("  [%3d/%3d] %2d건 | %5.1f초 | GVA %5.2f s/fit, AGHQ %5.2f s/fit | 누적 %.1f분",
                      n_ok, n_all, nrow(chunk), as.numeric(Sys.time()-t0,units="secs"),
                      mean(chunk$gva_time), mean(chunk$aghq_time),
                      as.numeric(Sys.time()-t_all, units="mins")))
    }
    if (!is.null(cl)) parallel::stopCluster(cl)
    gc()
  }
  out <- do.call(rbind, results)
  message(sprintf("\n완료: %d행, 총 %.1f분", nrow(out),
                  as.numeric(Sys.time()-t_all, units="mins")))
  out
}

## ==== 분석 ==================================================================
## (d, m, n) 별 대응차이와 그 표준오차
gap_table <- function(df) {
  ok <- df[df$gva_converged & df$aghq_converged, ]
  do.call(rbind, lapply(split(ok, list(ok$d, ok$m, ok$n), drop=TRUE), function(z) {
    g <- z$gva_sigvar - z$aghq_sigvar
    se <- sd(g)/sqrt(length(g))
    data.frame(d=z$d[1], m=z$m[1], n=z$n[1], reps=nrow(z),
               gap = mean(g), se = se, t = mean(g)/se,
               gva_time=mean(z$gva_time), aghq_time=mean(z$aghq_time),
               speedup=mean(z$aghq_time)/mean(z$gva_time),
               row.names=NULL)
  }))
}

## log|gap| ~ log n 회귀로 지수 alpha 추정 (gap ~ n^(-alpha))
scaling_exponent <- function(df) {
  ok <- df[df$gva_converged & df$aghq_converged, ]
  do.call(rbind, lapply(split(ok, list(ok$d, ok$m), drop=TRUE), function(z) {
    tab <- do.call(rbind, lapply(split(z, z$n), function(w)
      data.frame(n=w$n[1], gap=mean(w$gva_sigvar - w$aghq_sigvar))))
    tab <- tab[abs(tab$gap) > 1e-10, ]
    if (nrow(tab) < 3) return(NULL)
    fit <- lm(log(abs(gap)) ~ log(n), data=tab)
    ci <- tryCatch(confint(fit)[2, ], error=function(e) c(NA,NA))
    data.frame(d=z$d[1], m=z$m[1], n_points=nrow(tab),
               alpha = -unname(coef(fit)[2]),
               alpha_lo = -unname(ci[2]), alpha_hi = -unname(ci[1]),
               r2 = summary(fit)$r.squared, row.names=NULL)
  }))
}

load_scaling <- function(f = CFG$checkpoint) do.call(rbind, readRDS(f)$results)

## ==== 실행 ==================================================================
raw <- run_scaling()
gt  <- gap_table(raw)
write.csv(gt, CFG$csv_out, row.names = FALSE)
message(sprintf("요약 저장: %s", CFG$csv_out))

cat("\n=========== (d, m, n) 별 GVA-AGHQ 격차 ===========\n")
print(gt[order(gt$d, gt$m, gt$n), ], row.names = FALSE, digits = 3)

cat("\n=========== n 스케일링 지수  (격차 ~ n^-alpha) ===========\n")
cat("alpha = 1 이면 격차가 1/n, alpha = 0.5 면 1/sqrt(n)\n\n")
print(scaling_exponent(raw), row.names = FALSE, digits = 3)

cat("\n=========== m 효과 (n 에 걸쳐 평균) ===========\n")
ok <- raw[raw$gva_converged & raw$aghq_converged, ]
mm <- do.call(rbind, lapply(split(ok, list(ok$d, ok$m), drop=TRUE), function(z)
  data.frame(d=z$d[1], m=z$m[1], gap=mean(z$gva_sigvar - z$aghq_sigvar))))
print(reshape(mm, idvar="d", timevar="m", direction="wide"), row.names=FALSE, digits=3)

cat("\n=========== 수렴 실패 ===========\n")
cat(sprintf("GVA %d/%d, AGHQ %d/%d\n", sum(!raw$gva_converged), nrow(raw),
            sum(!raw$aghq_converged), nrow(raw)))
