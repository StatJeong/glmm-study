###############################################################################
##  epilepsy_fit.R
##  Ormerod & Wand (2012) Section 7.2 의 Poisson random intercept + slope 모형
##  (= Breslow & Clayton 1993 Model IV) 을 GVA 로 적합하고 AGHQ 와 비교한다.
##
##  논문 원문 모형 (preprint p.14):
##    y_ij | u0i, u1i ~ Poisson[ exp{ (b0 + u0i) + (b_visit + u1i) visit_j
##                                    + b_base  log(base_i/4)
##                                    + b_trt   trt_i
##                                    + b_bXt   log(base_i/4) x trt_i
##                                    + b_age   log(age_i) } ]
##    (u0i, u1i)' ~ N(0, Sigma),  Sigma = [s0^2, rho s0 s1; rho s0 s1, s1^2]
##
##    visit 코딩은 논문대로  visit1=-3, visit2=-1, visit3=1, visit4=3
###############################################################################

source("glmm_multi.R")

## ---------------------------------------------------------------------------
## 데이터 구성: 59명 x 4방문 완전균형 -> m x n 행렬
## ---------------------------------------------------------------------------
build_epilepsy <- function() {
  data(epilepsy, package = "HSAUR")
  d <- epilepsy
  d$visit <- 2 * as.integer(as.character(d$period)) - 5   # 1,2,3,4 -> -3,-1,1,3
  d$subject <- factor(d$subject)
  d <- d[order(d$subject, d$visit), ]

  m <- nlevels(d$subject); n <- 4L
  stopifnot(nrow(d) == m * n, all(table(d$subject) == n))

  mat <- function(v) matrix(v, nrow = m, ncol = n, byrow = TRUE)

  Y     <- mat(d$seizure.rate)
  visit <- mat(d$visit)
  lbase <- mat(log(d$base / 4))
  trt   <- mat(as.numeric(d$treatment == "Progabide"))
  lage  <- mat(log(d$age))
  bXt   <- lbase * trt

  ## Xl 순서가 곧 설계행렬 순서다.
  ##   Z = (1, Xl[[1]])          -> (1, visit)             K = 2
  ##   D = (1, Xl[[1]], ..., Xl[[5]]) -> 고정효과 6개
  list(dat = list(Y = Y, Xl = list(visit, lbase, trt, bXt, lage),
                  d = 2L, family = "poisson", nu = NULL, trials = 1),
       p = 6L,
       names = c("(Intercept)", "visit", "log(base/4)", "trt",
                 "log(base/4):trt", "log(age)"),
       raw = d)
}

## Sigma -> (sigma0, sigma1, rho) 로 보기 좋게
sigma_report <- function(S)
  c(sigma0 = sqrt(S[1,1]), sigma1 = sqrt(S[2,2]),
    rho = S[1,2] / sqrt(S[1,1] * S[2,2]))

## ---------------------------------------------------------------------------
## multi-start: 시작값을 흔들어 전역해인지 확인
## fit_gva_multi 는 내부에서 초기값을 고정 규칙으로 잡으므로,
## 여기서는 Sigma 초기값과 beta 교란으로 시작점을 바꾼다.
## ---------------------------------------------------------------------------
fit_gva_multistart <- function(dat, p, n_starts = 8, seed = 1, ...) {
  set.seed(seed)
  base <- fit_gva_multi(dat, p = p, accel = "squarem", ...)
  best <- base; best_elbo <- elbo_multi(base$beta, base$Sigma, dat, p = p)
  tab <- data.frame(start = 0, elbo = best_elbo, beta0 = base$beta[1],
                    s0 = sqrt(base$Sigma[1,1]), s1 = sqrt(base$Sigma[2,2]),
                    nG = base$n_G, conv = base$converged)

  for (k in seq_len(n_starts)) {
    ## 시작점 교란: beta 를 흔들고 Sigma 초기 스케일을 바꾼다
    pert <- rnorm(p, 0, 0.5)
    scl  <- exp(runif(1, log(0.05), log(3)))
    f <- tryCatch(
      fit_gva_multi(dat, p = p, accel = "squarem",
                    beta_start = base$beta + pert,
                    Sigma_start = diag(scl, 2), ...),
      error = function(e) NULL)
    if (is.null(f)) next
    e <- elbo_multi(f$beta, f$Sigma, dat, p = p)
    tab <- rbind(tab, data.frame(start = k, elbo = e, beta0 = f$beta[1],
                                 s0 = sqrt(f$Sigma[1,1]), s1 = sqrt(f$Sigma[2,2]),
                                 nG = f$n_G, conv = f$converged))
    if (is.finite(e) && e > best_elbo) { best <- f; best_elbo <- e }
  }
  list(fit = best, elbo = best_elbo, table = tab)
}
