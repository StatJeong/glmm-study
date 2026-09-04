###############################################################################
##  epilepsy_benchmark.R   --  Phase 3
##  epilepsy 데이터에서 K=1 (랜덤절편) vs K=2 (절편+기울기) 로
##  GVA / AGHQ / Laplace 의 시간과 추정치를 비교한다.
##
##  ---------------------------------------------------------------------------
##  측정상의 함정과 대응
##
##  구현 언어가 섞여 있다:
##      GVA (ours)          순수 R
##      AGHQ (ours)         순수 R
##      GLMMadaptive        C++
##      glmmTMB (Laplace)   C++ / TMB 자동미분
##  wall-clock 만 비교하면 알고리즘 차이와 언어 차이가 섞인다. 그래서
##
##   (1) 동일 언어 대조군을 둔다: GVA(R) vs 우리 AGHQ(R) 가 공정한 비교.
##       GLMMadaptive / glmmTMB 는 "실무 성능" 참조로 따로 표기.
##   (2) 주 지표를 비율로 둔다: 각 구현 안에서 K=1 -> K=2 시간이 몇 배 되는지.
##       언어 상수항이 상쇄되므로 K 확장의 효과만 남는다.
##
##  고정효과는 K=1, K=2 에서 동일하게 6개로 두어 K 효과만 분리한다.
###############################################################################

source("epilepsy_fit.R")
suppressMessages({ library(GLMMadaptive); library(glmmTMB) })

REPS   <- 5     # 반복 측정 후 중앙값 (단발 측정은 노이즈가 크다)
## 수렴 감사(별도 실행) 결과에 따른 설정:
##   우리 AGHQ  : n_quad=11 에서 이미 수렴 (K=1 변화 2.5e-7, K=2 변화 1.2e-5)
##   GLMMadaptive: nAGQ 를 31까지 올려도 beta0 가 1e-2 수준으로 계속 배회한다.
##                 수렴 지점이 없으므로 21 로 두고 불안정성을 결과에 표기한다.
N_QUAD    <- 11
N_AGQ_GLM <- 21

E <- build_epilepsy()
dd <- E$raw
dd$lbase <- log(dd$base / 4)
dd$trt   <- as.numeric(dd$treatment == "Progabide")
dd$lage  <- log(dd$age)
dd$bXt   <- dd$lbase * dd$trt

## K=1 / K=2 데이터 (고정효과는 동일, 랜덤효과만 다름)
dat_K1 <- E$dat; dat_K1$d <- 1L
dat_K2 <- E$dat                      # d = 2

## 주의: R 은 지연평가라 promise 를 force() 하면 첫 회만 실제로 계산되고
## 이후에는 캐시된 값이 즉시 반환된다(-> 중앙값이 0 이 된다).
## substitute + eval 로 매 반복 다시 평가해야 한다.
timeit <- function(expr, reps = REPS) {
  e <- substitute(expr); pf <- parent.frame()
  ts <- numeric(reps); val <- NULL
  for (i in seq_len(reps)) {
    t0 <- Sys.time()
    val <- eval(e, pf)
    ts[i] <- as.numeric(Sys.time() - t0, units = "secs")
  }
  list(median = median(ts), min = min(ts), value = val)
}

fx <- seizure.rate ~ visit + lbase + trt + bXt + lage
rows <- list()
est  <- list()

for (K in c(1, 2)) {
  dat <- if (K == 1) dat_K1 else dat_K2

  ## --- GVA (순수 R) ---
  r <- timeit(fit_gva_multi(dat, p = E$p, accel = "squarem",
                            max_iter = 1000, tol = 1e-10))
  g <- r$value
  rows[[length(rows)+1]] <- data.frame(K=K, method="GVA (R)", lang="R",
                                       time=r$median, tmin=r$min)
  est[[paste0("GVA_K",K)]] <- list(beta=g$beta, S=g$Sigma)

  ## --- AGHQ, 우리 구현 (순수 R) ---
  r <- timeit(fit_aghq_multi(dat, p = E$p, n_quad = N_QUAD, gva_start = TRUE))
  a <- r$value
  rows[[length(rows)+1]] <- data.frame(K=K, method="AGHQ (R, ours)", lang="R",
                                       time=r$median, tmin=r$min)
  est[[paste0("AGHQr_K",K)]] <- list(beta=a$beta, S=a$Sigma)

  ## --- GLMMadaptive (C++) ---
  rnd <- if (K == 1) ~ 1 | subject else ~ visit | subject
  r <- timeit(mixed_model(fx, random = rnd, data = dd,
                          family = poisson(), nAGQ = N_AGQ_GLM))
  m <- r$value
  rows[[length(rows)+1]] <- data.frame(K=K, method="GLMMadaptive (AGHQ)",
                                       lang="C++", time=r$median, tmin=r$min)
  est[[paste0("GLMMad_K",K)]] <- list(beta=fixef(m), S=m$D, ll=as.numeric(logLik(m)))

  ## --- glmmTMB (Laplace, C++) ---
  ff <- if (K == 1) update(fx, . ~ . + (1 | subject))
        else        update(fx, . ~ . + (visit | subject))
  r <- timeit(glmmTMB(ff, data = dd, family = poisson()))
  t <- r$value
  rows[[length(rows)+1]] <- data.frame(K=K, method="glmmTMB (Laplace)",
                                       lang="C++", time=r$median, tmin=r$min)
  est[[paste0("TMB_K",K)]] <- list(beta=fixef(t)$cond,
                                   S=VarCorr(t)$cond$subject,
                                   ll=as.numeric(logLik(t)))
}

bench <- do.call(rbind, rows)

## K=1 -> K=2 배율 (언어 상수항이 상쇄되는 주 지표)
w <- reshape(bench[, c("K","method","lang","time")], idvar=c("method","lang"),
             timevar="K", direction="wide")
names(w)[3:4] <- c("t_K1","t_K2")
w$ratio_K2_K1 <- w$t_K2 / w$t_K1

cat("\n=========== 절대 시간 (중앙값 ", REPS, "회) ===========\n", sep="")
print(bench[order(bench$K, bench$time), ], row.names = FALSE, digits = 3)

cat("\n=========== 주 지표: K=1 -> K=2 시간 배율 ===========\n")
cat("언어 상수항이 상쇄되므로 K 확장 비용만 남는다.\n\n")
print(w[order(w$ratio_K2_K1), ], row.names = FALSE, digits = 3)

cat("\n=========== 같은 언어(R) 안에서의 GVA vs AGHQ ===========\n")
for (K in c(1,2)) {
  g <- bench$time[bench$K==K & bench$method=="GVA (R)"]
  a <- bench$time[bench$K==K & bench$method=="AGHQ (R, ours)"]
  cat(sprintf("  K=%d : GVA %.4fs, AGHQ %.4fs  ->  AGHQ/GVA = %.2f배\n", K, g, a, a/g))
}

cat("\n=========== 추정치 일치 확인 (고정효과) ===========\n")
for (K in c(1,2)) {
  b <- rbind(GVA = est[[paste0("GVA_K",K)]]$beta,
             `AGHQ(R)` = est[[paste0("AGHQr_K",K)]]$beta,
             GLMMadaptive = est[[paste0("GLMMad_K",K)]]$beta,
             glmmTMB = est[[paste0("TMB_K",K)]]$beta)
  colnames(b) <- E$names
  cat(sprintf("\n--- K=%d ---\n", K)); print(round(b, 4))
  cat(sprintf("  최대 |GVA - GLMMadaptive| = %.5f\n",
              max(abs(b["GVA",] - b["GLMMadaptive",]))))
}

cat("\n=========== 분산성분 ===========\n")
for (K in c(1,2)) {
  cat(sprintf("\n--- K=%d ---\n", K))
  cat("GVA          :", round(as.numeric(est[[paste0("GVA_K",K)]]$S), 5), "\n")
  cat("AGHQ (R)     :", round(as.numeric(est[[paste0("AGHQr_K",K)]]$S), 5), "\n")
  cat("GLMMadaptive :", round(as.numeric(est[[paste0("GLMMad_K",K)]]$S), 5), "\n")
  cat("glmmTMB      :", round(as.numeric(est[[paste0("TMB_K",K)]]$S), 5), "\n")
}

saveRDS(list(bench = bench, wide = w, est = est), "epilepsy_benchmark2.rds")
cat("\n저장: epilepsy_benchmark.rds\n")
