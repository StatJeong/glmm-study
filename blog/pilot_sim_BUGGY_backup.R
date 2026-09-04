## pilot_sim.R -- Gamma/Poisson random-intercept GLMM: GVA vs exact MLE
## (functions recovered from .RData workspace)

library(statmod)   # gauss.quad

## ---------------------------------------------------------------------------
## Error handling for fits.
##
## The original code used  tryCatch(fit(...), error = function(e) rep(NA, k))
## which discarded the error text entirely, so a fatal, trivially-fixable error
## ("could not find function gauss.quad") was indistinguishable from a genuine
## numerical non-convergence. Both just showed up as NA.
##
## .fit_or_na() keeps the same NA-on-failure contract but (a) prints the real
## message to stderr -- which survives mclapply forking -- and (b) attaches it
## to the returned vector as attr(, "errmsg") so the caller can record it.
## Set options(pilot_sim.verbose_errors = FALSE) to silence the printing.
## ---------------------------------------------------------------------------
.fit_or_na <- function(expr, k, what = "fit", ctx = "") {
  tryCatch(expr, error = function(e) {
    msg <- conditionMessage(e)
    if (!identical(getOption("pilot_sim.verbose_errors", TRUE), FALSE))
      message(sprintf("[%s FAILED%s] %s", what,
                      if (nzchar(ctx)) paste0(" ", ctx) else "", msg))
    structure(rep(NA_real_, k), errmsg = msg)
  })
}

.errmsg_of <- function(x) {
  m <- attr(x, "errmsg")
  if (is.null(m)) NA_character_ else m
}

diagnose_convergence <- 
function (sim_out) 
{
    n_total <- nrow(sim_out)
    n_bad <- sum(!is.finite(sim_out$gva_beta0) | !is.finite(sim_out$exact_beta0))
    cat(sprintf("전체 %d개 중 fit 실패(NA/Inf) %d개 (%.1f%%)\n", 
        n_total, n_bad, 100 * n_bad/n_total))
    agg <- aggregate(cbind(gva_ok = is.finite(gva_beta0), exact_ok = is.finite(exact_beta0)) ~ 
        m + n, data = sim_out, FUN = function(x) sum(!x))
    names(agg)[3:4] <- c("gva_fail", "exact_fail")
    agg
}

e_step_gamma <- 
function (mu, lam, Di, sigma2, nu, n, n_iter = 50) 
{
    .newton_backtrack(mu, lam, "gamma", stat_i = Di, c1 = -nu * 
        n, sigma2 = sigma2, nu = nu, n_iter = n_iter)
}

e_step_poisson <- 
function (mu, lam, Bi, Ybar, sigma2, n_iter = 50) 
{
    .newton_backtrack(mu, lam, "poisson", stat_i = Bi, c1 = Ybar, 
        sigma2 = sigma2, n_iter = n_iter)
}

exact_negloglik <-
function (par, X, Y, family, nu = 2, n_quad = 11, gh = NULL)
{
    beta0 <- par[1]
    beta1 <- par[2]
    sigma2 <- exp(par[3])
    ## FIX: namespace-qualify gauss.quad. Unqualified, this silently required
    ## statmod to be attached; in any session without library(statmod) every
    ## call died instantly with "could not find function", and the tryCatch in
    ## one_rep turned that into NA. Also allow a precomputed `gh` so the nodes
    ## are not recomputed on every single optim() evaluation.
    if (is.null(gh)) gh <- statmod::gauss.quad(n_quad, kind = "hermite")
    xk <- gh$nodes
    wk <- gh$weights
    m <- nrow(X)
    n <- ncol(X)
    eta_fix <- beta0 + beta1 * X
    ETA_CLIP <- 30
    grad_hess_vec <- function(u_hat) {
        eta <- pmin(pmax(eta_fix + u_hat, -ETA_CLIP), ETA_CLIP)
        if (family == "poisson") {
            muv <- exp(eta)
            g <- rowSums(Y - muv) - u_hat/sigma2
            h <- -rowSums(muv) - 1/sigma2
        }
        else {
            ev <- exp(-eta)
            g <- rowSums(-nu + nu * Y * ev) - u_hat/sigma2
            h <- -rowSums(nu * Y * ev) - 1/sigma2
        }
        list(g = g, h = h)
    }
    u_hat <- rep(0, m)
    for (nt in 1:25) {
        gh_ <- grad_hess_vec(u_hat)
        g <- gh_$g
        h <- gh_$h
        bad <- !is.finite(g) | !is.finite(h) | h >= -1e-10
        h[bad] <- -1
        g[bad] <- 0
        step <- -g/h
        step[!is.finite(step)] <- 0
        step <- pmax(pmin(step, 5), -5)
        u_hat <- u_hat + step
        if (max(abs(step)) < 1e-10) 
            break
    }
    h_at_mode <- grad_hess_vec(u_hat)$h
    h_at_mode[!is.finite(h_at_mode) | h_at_mode >= 0] <- -1
    s <- sqrt(-1/h_at_mode)
    ll_mat <- matrix(0, nrow = m, ncol = n_quad)
    for (k in 1:n_quad) {
        u_k <- u_hat + sqrt(2) * s * xk[k]
        eta_k <- pmin(pmax(eta_fix + u_k, -ETA_CLIP), ETA_CLIP)
        loglik_k <- if (family == "poisson") 
            rowSums(dpois(Y, lambda = exp(eta_k), log = TRUE))
        else rowSums(dgamma(Y, shape = nu, rate = nu/exp(eta_k), 
            log = TRUE))
        ll_mat[, k] <- loglik_k - u_k^2/(2 * sigma2) + xk[k]^2
    }
    mx <- apply(ll_mat, 1, max)
    wsum <- as.vector(exp(ll_mat - mx) %*% wk)
    ll_total <- log(sqrt(2) * s) + mx + log(wsum) - 0.5 * log(2 * 
        pi * sigma2)
    out <- -sum(ll_total)
    if (!is.finite(out)) 
        out <- 1e+10
    out
}

fit_exact <-
function (X, Y, family, nu = 2, n_quad = 11, start = c(0, 0,
    0), gh = NULL)
{
    ## Compute the quadrature rule ONCE per fit instead of once per optim
    ## evaluation. Nodes depend only on n_quad, never on the parameters.
    if (is.null(gh)) gh <- statmod::gauss.quad(n_quad, kind = "hermite")
    opt <- optim(start, exact_negloglik, X = X, Y = Y, family = family,
        nu = nu, n_quad = n_quad, gh = gh, method = "BFGS")
    c(beta0 = opt$par[1], beta1 = opt$par[2], sigma2 = exp(opt$par[3]))
}

fit_gva <- 
function (X, Y, family = c("poisson", "gamma"), nu = 2, max_iter = 150, 
    tol = 1e-07) 
{
    family <- match.arg(family)
    m <- nrow(X)
    n <- ncol(X)
    Ybar <- rowSums(Y)
    Xvec <- as.vector(X)
    Yvec <- as.vector(Y)
    fam_glm <- if (family == "poisson") 
        poisson(link = "log")
    else Gamma(link = "log")
    init_fit <- suppressWarnings(glm(Yvec ~ Xvec, family = fam_glm))
    beta0 <- unname(coef(init_fit)[1])
    beta1 <- unname(coef(init_fit)[2])
    sigma2 <- 0.3
    mu <- rep(0, m)
    lam <- rep(0.1, m)
    for (it in 1:max_iter) {
        eta_fix <- beta0 + beta1 * X
        if (family == "poisson") {
            Bi <- rowSums(exp(eta_fix))
            es <- e_step_poisson(mu, lam, Bi, Ybar, sigma2)
        }
        else {
            Di <- rowSums(Y * exp(-eta_fix))
            es <- e_step_gamma(mu, lam, Di, sigma2, nu, n)
        }
        mu <- es$mu
        lam <- es$lam
        offset_i <- if (family == "poisson") 
            mu + lam/2
        else -mu + lam/2
        offset_full <- rep(offset_i, times = n)
        bnew <- robust_beta_update(Xvec, Yvec, offset_full, family, 
            nu, beta0, beta1)
        beta0_new <- unname(bnew["beta0"])
        beta1_new <- unname(bnew["beta1"])
        sigma2_target <- mean(mu^2 + lam)
        damp <- 0.5
        sigma2_new <- sigma2 + damp * (sigma2_target - sigma2)
        delta <- abs(beta0_new - beta0) + abs(beta1_new - beta1) + 
            abs(sigma2_new - sigma2)
        beta0 <- beta0_new
        beta1 <- beta1_new
        sigma2 <- sigma2_new
        if (delta < tol) 
            break
    }
    c(beta0 = unname(beta0), beta1 = unname(beta1), sigma2 = sigma2, 
        iters = it)
}

fit_gva_multistart <- 
function (X, Y, family = c("poisson", "gamma"), nu = 2, n_starts = 6, 
    seed_true = NULL, max_iter = 200, tol = 1e-08) 
{
    family <- match.arg(family)
    m <- nrow(X)
    n <- ncol(X)
    Ybar <- rowSums(Y)
    Xvec <- as.vector(X)
    Yvec <- as.vector(Y)
    fam_glm <- if (family == "poisson") 
        poisson(link = "log")
    else Gamma(link = "log")
    init_fit <- suppressWarnings(glm(Yvec ~ Xvec, family = fam_glm))
    b0n <- unname(coef(init_fit)[1])
    b1n <- unname(coef(init_fit)[2])
    starts <- list(c(beta0 = b0n, beta1 = b1n, sigma2 = 0.3), 
        c(beta0 = b0n * 0.5, beta1 = b1n, sigma2 = 0.1), c(beta0 = b0n * 
            1.5, beta1 = b1n, sigma2 = 0.5), c(beta0 = b0n, beta1 = b1n * 
            0.5, sigma2 = 1), c(beta0 = b0n, beta1 = b1n, sigma2 = 2), 
        c(beta0 = 0, beta1 = 0, sigma2 = 1))
    if (!is.null(seed_true)) 
        starts[[length(starts) + 1]] <- c(beta0 = seed_true[1], 
            beta1 = seed_true[2], sigma2 = seed_true[3])
    run_one <- function(st) {
        beta0 <- unname(st["beta0"])
        beta1 <- unname(st["beta1"])
        sigma2 <- unname(st["sigma2"])
        mu <- rep(0, m)
        lam <- rep(0.1, m)
        for (it in 1:max_iter) {
            eta_fix <- beta0 + beta1 * X
            if (family == "poisson") {
                Bi <- rowSums(exp(eta_fix))
                es <- e_step_poisson(mu, lam, Bi, Ybar, sigma2)
            }
            else {
                Di <- rowSums(Y * exp(-eta_fix))
                es <- e_step_gamma(mu, lam, Di, sigma2, nu, n)
            }
            mu <- es$mu
            lam <- es$lam
            offset_i <- if (family == "poisson") 
                mu + lam/2
            else -mu + lam/2
            offset_full <- rep(offset_i, times = n)
            bnew <- robust_beta_update(Xvec, Yvec, offset_full, 
                family, nu, beta0, beta1)
            beta0_new <- unname(bnew["beta0"])
            beta1_new <- unname(bnew["beta1"])
            sigma2_new <- sigma2 + 0.5 * (mean(mu^2 + lam) - 
                sigma2)
            delta <- abs(beta0_new - beta0) + abs(beta1_new - 
                beta1) + abs(sigma2_new - sigma2)
            beta0 <- beta0_new
            beta1 <- beta1_new
            sigma2 <- sigma2_new
            if (delta < tol) 
                break
        }
        elbo_val <- full_elbo(beta0, beta1, sigma2, X, Y, family, 
            nu)
        c(beta0 = beta0, beta1 = beta1, sigma2 = sigma2, iters = it, 
            elbo = elbo_val)
    }
    results <- t(sapply(starts, run_one))
    out <- as.data.frame(results)
    out$start_id <- seq_len(nrow(out))
    out[order(-out$elbo), c("start_id", "beta0", "beta1", "sigma2", 
        "iters", "elbo")]
}

full_elbo <- 
function (beta0, beta1, sigma2, X, Y, family, nu) 
{
    m <- nrow(X)
    n <- ncol(X)
    eta_fix <- beta0 + beta1 * X
    if (family == "poisson") {
        stat_i <- rowSums(exp(eta_fix))
        c1 <- rowSums(Y)
    }
    else {
        stat_i <- rowSums(Y * exp(-eta_fix))
        c1 <- -nu * n
    }
    mu <- rep(0, m)
    lam <- rep(0.1, m)
    es <- .newton_backtrack(mu, lam, family, stat_i, c1, sigma2, 
        nu = nu, n_iter = 100)
    mu <- es$mu
    lam <- es$lam
    E <- if (family == "poisson") 
        stat_i * exp(mu + lam/2)
    else nu * stat_i * exp(-mu + lam/2)
    per_group <- c1 * mu - E - (mu^2 + lam)/(2 * sigma2) - 0.5 * 
        log(sigma2) + 0.5 * log(lam)
    extra <- if (family == "poisson") 
        sum(Y * eta_fix)
    else -nu * sum(eta_fix)
    sum(per_group) + extra
}

robust_beta_update <- 
function (Xvec, Yvec, offset_full, family, nu, beta0, beta1, 
    n_iter = 50) 
{
    beta0 <- unname(beta0)
    beta1 <- unname(beta1)
    obj_and_gh <- function(b0, b1) {
        eta <- b0 + b1 * Xvec + offset_full
        if (family == "poisson") {
            mu <- exp(eta)
            obj <- sum(Yvec * eta - mu)
            g0 <- sum(Yvec - mu)
            g1 <- sum((Yvec - mu) * Xvec)
            h00 <- -sum(mu)
            h01 <- -sum(mu * Xvec)
            h11 <- -sum(mu * Xvec^2)
        }
        else {
            ev <- exp(-eta)
            obj <- sum(-eta - Yvec * ev)
            g0 <- sum(-1 + Yvec * ev)
            g1 <- sum((-1 + Yvec * ev) * Xvec)
            h00 <- -sum(Yvec * ev)
            h01 <- -sum(Yvec * ev * Xvec)
            h11 <- -sum(Yvec * ev * Xvec^2)
        }
        list(obj = obj, g0 = g0, g1 = g1, h00 = h00, h01 = h01, 
            h11 = h11)
    }
    cur <- obj_and_gh(beta0, beta1)
    cur_obj <- cur$obj
    for (it in 1:n_iter) {
        det <- cur$h00 * cur$h11 - cur$h01^2
        if (!is.finite(det) || det <= 1e-12 || cur$h00 >= 0) 
            break
        d0 <- -(cur$h11 * cur$g0 - cur$h01 * cur$g1)/det
        d1 <- -(-cur$h01 * cur$g0 + cur$h00 * cur$g1)/det
        step <- 1
        for (bt in 1:40) {
            b0_try <- beta0 + step * d0
            b1_try <- beta1 + step * d1
            cand <- obj_and_gh(b0_try, b1_try)
            if (is.finite(cand$obj) && cand$obj >= cur_obj - 
                1e-10) 
                break
            step <- step/2
        }
        beta0 <- beta0 + step * d0
        beta1 <- beta1 + step * d1
        cur <- obj_and_gh(beta0, beta1)
        if (abs(step * d0) < 1e-10 && abs(step * d1) < 1e-10) 
            break
        cur_obj <- cur$obj
    }
    c(beta0 = beta0, beta1 = beta1)
}

run_pilot_sim <- 
function (m_grid, n_grid, reps, beta0 = 1, beta1 = 0.5, sigma2 = 1, 
    family = c("poisson", "gamma"), nu = 2, n_quad = 11, verbose = TRUE) 
{
    family <- match.arg(family)
    results <- list()
    idx <- 1
    for (n in n_grid) for (m in m_grid) {
        for (r in 1:reps) {
            d <- simulate_data(m, n, beta0, beta1, sigma2, family, 
                nu)
            ctx <- sprintf("m=%d n=%d rep=%d", m, n, r)
            g <- .fit_or_na(fit_gva(d$X, d$Y, family, nu), 4, "GVA", ctx)
            e <- .fit_or_na(fit_exact(d$X, d$Y, family, nu, n_quad), 3, "exact", ctx)
            results[[idx]] <- data.frame(m = m, n = n, rep = r,
                gva_beta0 = g["beta0"], gva_beta1 = g["beta1"],
                gva_sigma2 = g["sigma2"], exact_beta0 = e["beta0"],
                exact_beta1 = e["beta1"], exact_sigma2 = e["sigma2"],
                gva_err = .errmsg_of(g), exact_err = .errmsg_of(e))
            idx <- idx + 1
        }
        if (verbose) 
            cat(sprintf("[%s] m=%d, n=%d 완료\n", family, m, 
                n))
    }
    out <- do.call(rbind, results)
    attr(out, "true_par") <- c(beta0 = beta0, beta1 = beta1, 
        sigma2 = sigma2)
    out
}

run_pilot_sim_parallel <- 
function (m_grid, n_grid, reps, beta0 = 1, beta1 = 0.5, sigma2 = 1, 
    family = c("poisson", "gamma"), nu = 2, n_quad = 11, mc.cores = parallel::detectCores(), 
    base_seed = 1, checkpoint_file = "checkpoint.rds") 
{
    family <- match.arg(family)
    is_windows <- .Platform$OS.type == "windows"
    if (is_windows) {
        cl <- parallel::makeCluster(mc.cores)
        on.exit(parallel::stopCluster(cl), add = TRUE)
        parallel::clusterExport(cl, varlist = c("simulate_data",
            "fit_gva", "fit_exact", "e_step_poisson", "e_step_gamma",
            ".newton_backtrack", ".obj_grad_hess", "exact_negloglik",
            "robust_beta_update", ".fit_or_na", ".errmsg_of"),
            envir = globalenv())
    }
    ## Quadrature nodes depend only on n_quad, so build them once here in the
    ## parent. Passing them down also means a worker never has to resolve
    ## gauss.quad itself.
    gh_shared <- statmod::gauss.quad(n_quad, kind = "hermite")
    one_rep <- function(rr, mm, nn) {
        set.seed(base_seed * 1e+05 + mm * 1000 + nn * 100 + rr)
        d <- simulate_data(mm, nn, beta0, beta1, sigma2, family,
            nu)
        ctx <- sprintf("m=%d n=%d rep=%d", mm, nn, rr)
        t0 <- Sys.time()
        g <- .fit_or_na(fit_gva(d$X, d$Y, family, nu), 4, "GVA", ctx)
        t_gva <- as.numeric(Sys.time() - t0, units = "secs")
        t0 <- Sys.time()
        e <- .fit_or_na(fit_exact(d$X, d$Y, family, nu, n_quad, gh = gh_shared),
            3, "exact", ctx)
        t_exact <- as.numeric(Sys.time() - t0, units = "secs")
        data.frame(m = mm, n = nn, rep = rr, gva_beta0 = g["beta0"],
            gva_beta1 = g["beta1"], gva_sigma2 = g["sigma2"],
            exact_beta0 = e["beta0"], exact_beta1 = e["beta1"],
            exact_sigma2 = e["sigma2"], gva_time = t_gva, exact_time = t_exact,
            gva_err = .errmsg_of(g), exact_err = .errmsg_of(e))
    }
    all_results <- list()
    done_cells <- character(0)
    if (file.exists(checkpoint_file)) {
        prev <- readRDS(checkpoint_file)
        all_results <- prev$all_results
        done_cells <- prev$done_cells
        cat(sprintf("체크포인트 발견: %d칸 이미 완료됨, 이어서 진행합니다.\n", 
            length(done_cells)))
    }
    n_cells <- length(m_grid) * length(n_grid)
    cat(sprintf("코어 %d개 사용(큰 셀은 자동으로 줄임). (m,n) 조합 %d개 x reps %d개 -- 전체 %d칸\n", 
        mc.cores, n_cells, reps, n_cells))
    grp_idx <- length(all_results) + 1
    t_start_all <- Sys.time()
    for (n in n_grid) for (m in m_grid) {
        cell_key <- paste0("m", m, "_n", n)
        if (cell_key %in% done_cells) 
            next
        cell_size <- m * n
        cell_cores <- if (cell_size > 2e+07) 
            max(1, min(mc.cores, 2))
        else if (cell_size > 5e+06) 
            max(1, min(mc.cores, 4))
        else mc.cores
        t0 <- Sys.time()
        if (is_windows) {
            res_list <- parallel::parLapply(cl, 1:reps, one_rep, 
                mm = m, nn = n)
        }
        else {
            res_list <- parallel::mclapply(1:reps, one_rep, mm = m, 
                nn = n, mc.cores = cell_cores, mc.preschedule = FALSE)
        }
        ## mclapply does NOT throw when a worker dies (segfault, OOM, killed):
        ## it returns a try-error in that slot. rbind-ing those silently
        ## corrupts the cell, so surface them and substitute an NA row.
        crashed <- vapply(res_list, function(z) inherits(z, "try-error") ||
                            !is.data.frame(z), logical(1))
        if (any(crashed)) {
            for (ci in which(crashed))
                message(sprintf("[WORKER CRASHED] m=%d n=%d rep=%d: %s",
                                m, n, ci, paste(as.character(res_list[[ci]]),
                                                collapse = " ")))
            template <- if (any(!crashed)) res_list[[which(!crashed)[1]]] else NULL
            for (ci in which(crashed)) {
                if (is.null(template)) { res_list[[ci]] <- NULL; next }
                row <- template[1, ]
                row[] <- NA
                row$m <- m; row$n <- n; row$rep <- ci
                row$exact_err <- "worker crashed"
                res_list[[ci]] <- row
            }
            res_list <- Filter(Negate(is.null), res_list)
        }
        cell_df <- do.call(rbind, res_list)
        all_results[[grp_idx]] <- cell_df
        done_cells <- c(done_cells, cell_key)
        elapsed <- as.numeric(Sys.time() - t0, units = "secs")
        total_elapsed <- as.numeric(Sys.time() - t_start_all, 
            units = "secs")
        avg_gva <- mean(cell_df$gva_time, na.rm = TRUE)
        avg_exact <- mean(cell_df$exact_time, na.rm = TRUE)
        cat(sprintf("[%d/%d 완료] %s m=%-7d n=%-5d (코어 %d개 사용) | 이 조합 %.1fs (GVA 평균 %.1fs/fit, exact 평균 %.1fs/fit) | 누적 %.1fs\n", 
            grp_idx, n_cells, family, m, n, cell_cores, elapsed, 
            avg_gva, avg_exact, total_elapsed))
        grp_idx <- grp_idx + 1
        saveRDS(list(all_results = all_results, done_cells = done_cells, 
            true_par = c(beta0 = beta0, beta1 = beta1, sigma2 = sigma2)), 
            checkpoint_file)
        cat(sprintf("  -> 체크포인트 저장됨 (%s)\n", 
            checkpoint_file))
        gc()
    }
    out <- do.call(rbind, all_results)
    attr(out, "true_par") <- c(beta0 = beta0, beta1 = beta1, 
        sigma2 = sigma2)
    out
}

simulate_data <- 
function (m, n, beta0, beta1, sigma2, family = c("poisson", "gamma"), 
    nu = 2) 
{
    family <- match.arg(family)
    X <- matrix(rnorm(m * n), nrow = m, ncol = n)
    U <- rnorm(m, 0, sqrt(sigma2))
    eta <- beta0 + beta1 * X + U
    mu <- exp(eta)
    if (family == "poisson") {
        Y <- matrix(rpois(m * n, lambda = mu), nrow = m, ncol = n)
    }
    else {
        Y <- matrix(rgamma(m * n, shape = nu, rate = nu/mu), 
            nrow = m, ncol = n)
    }
    list(X = X, Y = Y)
}

summarize_bias <- 
function (sim_out) 
{
    tp <- attr(sim_out, "true_par")
    agg <- aggregate(cbind(gva_beta0, gva_beta1, gva_sigma2, 
        exact_beta0, exact_beta1, exact_sigma2) ~ m + n, data = sim_out, 
        FUN = function(x) mean(x, na.rm = TRUE))
    agg$gva_bias_beta0 <- agg$gva_beta0 - tp["beta0"]
    agg$gva_bias_sigma2 <- agg$gva_sigma2 - tp["sigma2"]
    agg$exact_bias_beta0 <- agg$exact_beta0 - tp["beta0"]
    agg$exact_bias_sigma2 <- agg$exact_sigma2 - tp["sigma2"]
    agg[order(agg$n, agg$m), ]
}

summarize_full_estimates <- 
function (sim_out) 
{
    tp <- attr(sim_out, "true_par")
    agg <- aggregate(cbind(gva_beta0, gva_beta1, gva_sigma2, 
        exact_beta0, exact_beta1, exact_sigma2) ~ m + n, data = sim_out, 
        FUN = function(x) mean(x, na.rm = TRUE))
    agg <- agg[order(agg$n, agg$m), ]
    cat(sprintf("참값: beta0=%.2f, beta1=%.2f, sigma2=%.2f\n\n", 
        tp["beta0"], tp["beta1"], tp["sigma2"]))
    agg
}

summarize_multistart_bias <- 
function (sim_out) 
{
    tp <- attr(sim_out, "true_par")
    agg <- aggregate(cbind(ms_beta0, ms_beta1, ms_sigma2) ~ m + 
        n, data = sim_out, FUN = function(x) mean(x, na.rm = TRUE))
    agg$ms_bias_beta0 <- agg$ms_beta0 - tp["beta0"]
    agg$ms_bias_sigma2 <- agg$ms_sigma2 - tp["sigma2"]
    agg[order(agg$n, agg$m), ]
}

verify_with_multistart <- 
function (m_grid, n_grid, reps, beta0 = 1, beta1 = 0.5, sigma2 = 1, 
    family = c("poisson", "gamma"), nu = 2, n_starts_seed_true = TRUE) 
{
    family <- match.arg(family)
    results <- list()
    idx <- 1
    for (n in n_grid) for (m in m_grid) {
        t0 <- Sys.time()
        for (r in 1:reps) {
            set.seed(9e+05 + m * 1000 + n * 100 + r)
            d <- simulate_data(m, n, beta0, beta1, sigma2, family, 
                nu)
            ms <- tryCatch(fit_gva_multistart(d$X, d$Y, family,
                nu, seed_true = if (n_starts_seed_true)
                  c(beta0, beta1, sigma2)
                else NULL), error = function(e) {
                  ## do not swallow: report which rep failed and why
                  message(sprintf("[multistart FAILED m=%d n=%d rep=%d] %s",
                                  m, n, r, conditionMessage(e)))
                  NULL
                })
            if (is.null(ms)) {
                best <- list(beta0 = NA, beta1 = NA, sigma2 = NA)
            }
            else {
                best <- list(beta0 = ms$beta0[1], beta1 = ms$beta1[1], 
                  sigma2 = ms$sigma2[1])
            }
            results[[idx]] <- data.frame(m = m, n = n, rep = r, 
                ms_beta0 = best$beta0, ms_beta1 = best$beta1, 
                ms_sigma2 = best$sigma2)
            idx <- idx + 1
        }
        cat(sprintf("[%s multistart] m=%d n=%d 완료 (%.1fs)\n", 
            family, m, n, as.numeric(Sys.time() - t0, units = "secs")))
    }
    out <- do.call(rbind, results)
    attr(out, "true_par") <- c(beta0 = beta0, beta1 = beta1, 
        sigma2 = sigma2)
    out
}


## --- internal helpers (dot-prefixed; were missed by ls() without all.names) ---

.newton_backtrack <- 
function (mu, lam, family, stat_i, c1, sigma2, nu = NULL, n_iter = 50) 
{
    obj_fun <- function(mu, lam) {
        if (family == "poisson") {
            E <- stat_i * exp(mu + lam/2)
            c1 * mu - E - (mu^2 + lam)/(2 * sigma2) + 0.5 * log(lam)
        }
        else {
            E <- nu * stat_i * exp(-mu + lam/2)
            c1 * mu - E - (mu^2 + lam)/(2 * sigma2) + 0.5 * log(lam)
        }
    }
    cur_obj <- obj_fun(mu, lam)
    for (it in 1:n_iter) {
        gh <- .obj_grad_hess(mu, lam, family, stat_i, c1, sigma2, 
            nu)
        gmu <- gh$gmu
        glam <- gh$glam
        H11 <- gh$H11
        H12 <- gh$H12
        H22 <- gh$H22
        det <- H11 * H22 - H12^2
        use_newton <- is.finite(det) & (det > 1e-12) & (H11 < 
            0)
        dmu <- numeric(length(mu))
        dlam <- numeric(length(mu))
        dmu[use_newton] <- -(H22[use_newton] * gmu[use_newton] - 
            H12[use_newton] * glam[use_newton])/det[use_newton]
        dlam[use_newton] <- -(-H12[use_newton] * gmu[use_newton] + 
            H11[use_newton] * glam[use_newton])/det[use_newton]
        step0 <- 0.05
        dmu[!use_newton] <- step0 * gmu[!use_newton]
        dlam[!use_newton] <- step0 * glam[!use_newton]
        step <- rep(1, length(mu))
        for (bt in 1:30) {
            mu_try <- mu + step * dmu
            lam_try <- lam + step * dlam
            lam_try[lam_try <= 1e-08] <- NA
            new_obj <- obj_fun(mu_try, lam_try)
            improved <- is.finite(new_obj) & (new_obj >= cur_obj - 
                1e-10)
            if (all(improved)) 
                break
            step[!improved] <- step[!improved]/2
            if (max(step) < 1e-10) 
                break
        }
        mu_try[is.na(mu_try)] <- mu[is.na(mu_try)]
        lam_try[is.na(lam_try) | lam_try <= 1e-08] <- lam[is.na(lam_try) | 
            lam_try <= 1e-08]
        new_obj <- obj_fun(mu_try, lam_try)
        keep <- is.finite(new_obj) & (new_obj >= cur_obj - 1e-10)
        mu[keep] <- mu_try[keep]
        lam[keep] <- lam_try[keep]
        cur_obj[keep] <- new_obj[keep]
        if (max(abs(step * dmu)) < 1e-09 && max(abs(step * dlam)) < 
            1e-09) 
            break
    }
    list(mu = mu, lam = lam)
}

.obj_grad_hess <- 
function (mu, lam, family, Bi_or_Di, c1, sigma2, nu = NULL) 
{
    if (family == "poisson") {
        E <- Bi_or_Di * exp(mu + lam/2)
        gmu <- c1 - E - mu/sigma2
        glam <- -0.5 * E - 1/(2 * sigma2) + 1/(2 * lam)
        H11 <- -E - 1/sigma2
        H12 <- -0.5 * E
        H22 <- -0.25 * E - 1/(2 * lam^2)
    }
    else {
        E <- nu * Bi_or_Di * exp(-mu + lam/2)
        gmu <- c1 + E - mu/sigma2
        glam <- -0.5 * E - 1/(2 * sigma2) + 1/(2 * lam)
        H11 <- -E - 1/sigma2
        H12 <- 0.5 * E
        H22 <- -0.25 * E - 1/(2 * lam^2)
    }
    list(gmu = gmu, glam = glam, H11 = H11, H12 = H12, H22 = H22, 
        E = E)
}

