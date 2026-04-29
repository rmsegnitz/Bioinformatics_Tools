#####################################################################
# parallelized parameter tuning for regularized SEM models in regsem
#####################################################################

# Author: Max Segnitz, msegnitz@uw.edu
# Started February 2026
#
# © Richard M Segnitz 2026
# License: This software is licensed under GNU General Public License, and may
# be modified and/or redistributed under GNU GPL version 3 or later. License details
# can be found in the accompanying this script, or at  (http://www.gnu.org/licenses/).
#
# DESCRIPTION:
# A fairly straightforward parallel implimention of the parameter tuning (eg for lasso) in regsem.
# This is designed to behave almost identically to the cv_regsem function, but to receive specified
# number of cores over which to distribute the parameter search. The function also returns a limited
# number of additional outputs which may be useful, eg some canned plots exploring model fits. 

################################################################################################################

par_cv_regsem <- function(
    model,
    type          = "lasso",
    pars_pen      = "regressions",
    metric        = "BIC",
    mult.start    = FALSE,
    multi.iter    = 10,
    n.lambda      = 40,
    jump          = 0.01,
    lambda.start  = 0,
    alpha         = 0,
    gamma         = 3.7,
    fit.ret       = c("rmsea", "BIC", "chisq"),
    optMethod     = "rsolnp",
    gradFun       = "ram",
    hessFun       = "none",
    Start         = "lavaan",
    subOpt        = "nlminb",
    diff_par      = NULL,
    LB            = -Inf,
    UB            = Inf,
    par.lim       = c(-Inf, Inf),
    block         = TRUE,
    full          = TRUE,
    calc          = "normal",
    max.iter      = 2000,
    tol           = 1e-5,
    round         = 3,
    solver        = FALSE,
    quasi         = FALSE,
    solver.maxit  = 5,
    alpha.inc     = FALSE,
    step          = 0.1,
    momentum      = FALSE,
    step.ratio    = FALSE,
    line.search   = FALSE,
    nlminb.control = list(),
    missing       = "listwise",
    random.alpha  = 0.5,
    prerun        = FALSE,
    par_threshold = 20L,
    n.cores       = parallel::detectCores() - 1
) {
  
  require(regsem)
  require(parallel)
  
  # ---- build full lambda grid -----------------------------------------------
  lambda_grid <- round(
    seq(lambda.start, lambda.start + (n.lambda - 1) * jump, by = jump), 10
  )
  
  # ---- bundle all regsem/multi_optim args for workers -----------------------
  worker_args <- list(
    model          = model,
    type           = type,
    pars_pen       = pars_pen,
    alpha          = alpha,
    gamma          = gamma,
    optMethod      = optMethod,
    gradFun        = gradFun,
    hessFun        = hessFun,
    Start          = Start,
    subOpt         = subOpt,
    diff_par       = diff_par,
    LB             = LB,
    UB             = UB,
    par.lim        = par.lim,
    block          = block,
    full           = full,
    calc           = calc,
    max.iter       = max.iter,
    tol            = tol,
    round          = round,
    solver         = solver,
    quasi          = quasi,
    solver.maxit   = solver.maxit,
    alpha.inc      = alpha.inc,
    step           = step,
    momentum       = momentum,
    step.ratio     = step.ratio,
    line.search    = line.search,
    nlminb.control = nlminb.control,
    missing        = missing,
    random.alpha   = random.alpha,
    prerun         = prerun,
    mult.start     = mult.start,
    multi.iter     = multi.iter,
    fit.ret        = fit.ret
  )
  
  # ---- single lambda worker -------------------------------------------------
  # mirrors the cv_parallel() inner function from cv_regsem source
  .worker_single <- function(lam, args) {
    tryCatch(
      withCallingHandlers({
        if (!args$mult.start) {
          out <- regsem::regsem(
            model          = args$model,
            lambda         = lam,
            type           = args$type,
            pars_pen       = args$pars_pen,
            alpha          = args$alpha,
            gamma          = args$gamma,
            optMethod      = args$optMethod,
            gradFun        = args$gradFun,
            hessFun        = args$hessFun,
            Start          = args$Start,
            subOpt         = args$subOpt,
            diff_par       = args$diff_par,
            LB             = args$LB,
            UB             = args$UB,
            par.lim        = args$par.lim,
            block          = args$block,
            full           = args$full,
            calc           = args$calc,
            max.iter       = args$max.iter,
            tol            = args$tol,
            round          = args$round,
            solver         = args$solver,
            quasi          = args$quasi,
            solver.maxit   = args$solver.maxit,
            alpha.inc      = args$alpha.inc,
            step           = args$step,
            momentum       = args$momentum,
            step.ratio     = args$step.ratio,
            line.search    = args$line.search,
            nlminb.control = args$nlminb.control,
            missing        = args$missing,
            random.alpha   = args$random.alpha,
            prerun         = args$prerun
          )
        } else {
          out <- regsem::multi_optim(
            model          = args$model,
            lambda         = lam,
            max.try        = args$multi.iter,
            type           = args$type,
            pars_pen       = args$pars_pen,
            alpha          = args$alpha,
            gamma          = args$gamma,
            optMethod      = args$optMethod,
            gradFun        = args$gradFun,
            hessFun        = args$hessFun,
            LB             = args$LB,
            UB             = args$UB,
            par.lim        = args$par.lim,
            block          = args$block,
            full           = args$full,
            max.iter       = args$max.iter,
            tol            = args$tol,
            round          = args$round,
            solver         = args$solver,
            quasi          = args$quasi,
            solver.maxit   = args$solver.maxit,
            alpha.inc      = args$alpha.inc,
            step           = args$step,
            momentum       = args$momentum,
            step.ratio     = args$step.ratio,
            line.search    = args$line.search,
            nlminb.control = args$nlminb.control,
            random.alpha   = args$random.alpha,
            prerun         = args$prerun,
            diff_par       = args$diff_par
          )
        }
        
        # compute fit indices the same way cv_regsem does internally
        fitt <- tryCatch(
          regsem::fit_indices(out, CV = FALSE)$fits[args$fit.ret],
          error = function(e) setNames(rep(NA_real_, length(args$fit.ret)),
                                       args$fit.ret)
        )
        
        list(
          lambda      = lam,
          conv        = out$convergence,
          fits        = fitt,
          coefficients = out$coefficients,
          df          = out$df
        )
      },
      warning = function(w) invokeRestart("muffleWarning")
      ),
      error = function(e) {
        message(sprintf("\nFailed at lambda=%.4f: %s", lam, conditionMessage(e)))
        NULL
      }
    )
  }
  
  # ---- serial fallback for small grids --------------------------------------
  .run_serial <- function(lambda_grid, worker_args, fit.ret, metric,
                          n.lambda, par_threshold, serial_msg) {
    message(serial_msg)
    pb  <- txtProgressBar(min = 0, max = length(lambda_grid), style = 3)
    res <- vector("list", length(lambda_grid))
    for (i in seq_along(lambda_grid)) {
      res[[i]] <- .worker_single(lambda_grid[[i]], worker_args)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    res
  }
  
  # ---- progress file setup --------------------------------------------------
  prog_file <- tempfile("par_cv_regsem_progress_")
  writeLines("0", prog_file)
  on.exit(unlink(prog_file), add = TRUE)
  
  .poll_progress <- function(prog_file, n_total, start_time) {
    width  <- 40
    last_n <- -1L
    while (TRUE) {
      Sys.sleep(0.5)
      n_done <- tryCatch(as.integer(readLines(prog_file, warn = FALSE)),
                         error = function(e) 0L)
      if (is.na(n_done)) n_done <- 0L
      if (n_done != last_n) {
        pct     <- n_done / n_total
        filled  <- round(pct * width)
        bar     <- paste0(strrep("=", filled),
                          if (filled < width) ">" else "",
                          strrep(" ", max(0, width - filled - 1)))
        elapsed <- round(difftime(Sys.time(), start_time, units = "secs"))
        eta     <- if (pct > 0) round(as.numeric(elapsed) / pct - as.numeric(elapsed)) else "?"
        cat(sprintf("\r  [%s] %d/%d  elapsed: %ds  ETA: %ss   ",
                    bar, n_done, n_total, elapsed, eta))
        flush.console()
        last_n <- n_done
      }
      if (n_done >= n_total) break
    }
    cat("\n")
  }
  
  # ---- parallel worker wrapper (adds progress tick) -------------------------
  .worker_par <- function(lam, args, prog_file) {
    result <- .worker_single(lam, args)
    tryCatch({
      n <- as.integer(readLines(prog_file, warn = FALSE)) + 1L
      writeLines(as.character(n), prog_file)
    }, error = function(e) NULL)
    result
  }
  
  # ---- dispatch: serial or parallel -----------------------------------------
  start_time <- Sys.time()
  
  if (n.lambda < par_threshold) {
    raw_results <- .run_serial(
      lambda_grid, worker_args, fit.ret, metric, n.lambda, par_threshold,
      sprintf("par_cv_regsem: n.lambda=%d is below threshold (%d), running serially.",
              n.lambda, par_threshold)
    )
  } else {
    n.cores <- min(n.cores, length(lambda_grid))
    message(sprintf(
      "par_cv_regsem: parallelising %d lambda values across %d cores...",
      n.lambda, n.cores
    ))
    
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n.cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl,
                              varlist = c("worker_args", ".worker_single", ".worker_par", "prog_file"),
                              envir   = environment()
      )
      parallel::clusterEvalQ(cl, library(regsem))
      raw_results <- parallel::parLapply(
        cl, lambda_grid,
        function(lam) .worker_par(lam, worker_args, prog_file)
      )
      .poll_progress(prog_file, n.lambda, start_time)
    } else {
      poll_pid <- parallel::mcparallel(
        .poll_progress(prog_file, n.lambda, start_time)
      )
      raw_results <- parallel::mclapply(
        lambda_grid,
        function(lam) .worker_par(lam, worker_args, prog_file),
        mc.cores       = n.cores,
        mc.preschedule = FALSE
      )
      parallel::mccollect(poll_pid)
    }
  }
  
  elapsed_total <- round(difftime(Sys.time(), start_time, units = "secs"))
  message(sprintf("par_cv_regsem: completed in %ds.", elapsed_total))
  
  # ---- filter failed lambda values ------------------------------------------
  results_ok <- Filter(Negate(is.null), raw_results)
  n_converged <- length(results_ok)
  
  if (n_converged == 0)
    stop("par_cv_regsem: no lambda values returned results. ",
         "Check model specification and lambda range.")
  
  message(sprintf("par_cv_regsem: %d / %d lambda values returned results.",
                  n_converged, n.lambda))
  
  # ---- assemble output matrices matching cv_regsem format -------------------
  fits_mat <- do.call(rbind, lapply(results_ok, function(r) {
    c(lambda = r$lambda, conv = r$conv, r$fits)
  }))
  
  pars_mat <- do.call(rbind, lapply(results_ok, function(r) r$coefficients))
  rownames(pars_mat) <- NULL
  
  # sort by lambda
  ord      <- order(fits_mat[, "lambda"])
  fits_mat <- fits_mat[ord, , drop = FALSE]
  pars_mat <- pars_mat[ord, , drop = FALSE]
  
  # ---- identify best converged solution -------------------------------------
  conv     <- fits_mat[, "conv"]
  fit_vals <- fits_mat[, metric]
  best_idx <- which.min(ifelse(conv == 0 & !is.nan(fit_vals) & !is.na(fit_vals),
                               fit_vals, Inf))
  best_lambda <- fits_mat[best_idx, "lambda"]
  best_pars   <- pars_mat[best_idx, ]
  
  # ---- return ---------------------------------------------------------------
  structure(
    list(
      fits           = fits_mat,
      parameters     = pars_mat,
      final_pars     = best_pars,
      best_lambda    = best_lambda,
      best_idx       = best_idx,
      pars_pen       = pars_pen,
      metric         = metric,
      n_total_lambda = n.lambda,
      n_converged    = n_converged,
      call           = match.call()
    ),
    class = c("par_cv_regsem", "list")
  )
}

# ---- S3 methods --------------------------------------------------------------

print.par_cv_regsem <- function(x, ...) {
  cat("par_cv_regsem result\n")
  cat(sprintf("  Lambda range:   %.4f -- %.4f\n",
              min(x$fits[, "lambda"]), max(x$fits[, "lambda"])))
  cat(sprintf("  Converged:      %d / %d lambda values\n",
              x$n_converged, x$n_total_lambda))
  cat(sprintf("  Best lambda:    %.4f  (metric: %s)\n", x$best_lambda, x$metric))
  cat(sprintf("  Best %s:     %.2f\n", x$metric, x$fits[x$best_idx, x$metric]))
  cat("\nBest-fit parameters:\n")
  print(round(x$final_pars, 4))
  invisible(x)
}

summary.par_cv_regsem <- function(object, ...) {
  fits  <- as.data.frame(object$fits)
  conv  <- fits$conv
  valid <- conv == 0 & !is.na(fits[[object$metric]])
  delta <- fits[[object$metric]] - min(fits[[object$metric]][valid], na.rm = TRUE)
  
  cat("par_cv_regsem summary\n\n")
  cat(sprintf("Lambda range:      %.4f to %.4f\n",
              min(fits$lambda), max(fits$lambda)))
  cat(sprintf("Converged:         %d / %d lambda values\n",
              object$n_converged, object$n_total_lambda))
  cat(sprintf("Metric:            %s\n", object$metric))
  cat(sprintf("Best lambda:       %.4f\n", object$best_lambda))
  cat(sprintf("Best %s:        %.2f\n", object$metric,
              fits[[object$metric]][object$best_idx]))
  lams_in_window <- fits$lambda[valid & delta <= 6]
  cat(sprintf("Lambda within d%s <= 6:  %.4f to %.4f\n",
              object$metric,
              min(lams_in_window, na.rm = TRUE),
              max(lams_in_window, na.rm = TRUE)))
  cat("\nNon-zero parameters at best lambda:\n")
  nz <- object$final_pars[abs(object$final_pars) > 1e-6]
  print(round(nz, 4))
  invisible(object)
}

plot.par_cv_regsem <- function(x, bic_thresh = 6, ...) {
  fits  <- as.data.frame(x$fits)
  
  # filter non-converged
  n_total   <- nrow(fits)
  fits      <- fits[!is.na(fits[[x$metric]]) & fits$conv == 0, ]
  n_removed <- n_total - nrow(fits)
  
  delta <- fits[[x$metric]] - min(fits[[x$metric]], na.rm = TRUE)
  col   <- ifelse(delta <= 2,          "#1B9E77",
                  ifelse(delta <= bic_thresh, "#D95F02", "#BDBDBD"))
  
  op <- par(mar = c(5.1, 4.1, 5.5, 2.1))
  on.exit(par(op), add = TRUE)
  
  plot(fits$lambda, fits[[x$metric]],
       col  = col, pch = 19, cex = 0.8,
       xlab = expression(lambda),
       ylab = x$metric,
       main = sprintf("%s across lambda (par_cv_regsem)", x$metric), ...)
  abline(v = x$best_lambda,
         lty = 2, col = "black")
  abline(h = min(fits[[x$metric]], na.rm = TRUE) + bic_thresh,
         lty = 3, col = "grey50")
  
  legend("topleft",
         inset  = c(0, -0.18),
         xpd    = TRUE,
         horiz  = TRUE,
         bty    = "n",
         cex    = 0.8,
         legend = c(expression(Delta*"BIC 0-2"),
                    expression(Delta*"BIC 2-6"),
                    expression(Delta*"BIC 6+")),
         col    = c("#1B9E77", "#D95F02", "#BDBDBD"),
         pch    = 19)
  
  if (n_removed > 0)
    mtext(sprintf("%d non-converged lambda value%s omitted",
                  n_removed, if (n_removed > 1) "s" else ""),
          side = 1, line = 4, cex = 0.75, col = "grey40", adj = 1)
  
  invisible(x)
}