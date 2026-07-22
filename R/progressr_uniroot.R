uniroot_with_progress <- function(
  fun,
  target,
  bounds,
  abs_tol = 1e-6,
  maxiter = 100
) {
  last_val <- NA
  first_interval <- NA
  last_reported_progress <- 0

  left_val_dist <- fun(bounds[1]) - target
  right_val_dist <- fun(bounds[2]) - target

  if (left_val_dist * right_val_dist > 0) {
    stop("The bounds do not bracket a fitting solution for the target value.")
  }

  p <- progressr::progressor(100)

  optimization_goal <- function(target_val, target) {
    logger::log_debug("Evaluating at ", target_val)
    est <- fun(target_val)
    error <- est - target
    if (!is.na(last_val)) {
      if (est - last_val == 0) {
        progress <- 1
      } else {
        conv_crit <- .Machine$double.eps * abs(target_val) + abs_tol / 2
        interval <- abs(last_val - target_val)
        if (is.na(first_interval)) {
          first_interval <<- interval
        }
        progress <- (log(interval) - log(first_interval)) /
          (log(conv_crit) - log(first_interval))
      }
      logger::log_debug(
        "Current estimate: ", est,
        ", error: ", error,
        ", progress: ", progress
      )
      progress_percent <- min(ceiling(100 * progress), 100)
      add_progress <- progress_percent - last_reported_progress
      if (add_progress > 0) {
        message <- paste(
          "Est.: ", round(est, 4), "(", round(target_val, 4),
          ")"
        )
        p(message, amount = add_progress)
      }
    }
    last_val <<- target_val
    error
  }

  results <- stats::uniroot(optimization_goal,
    interval = bounds,
    target = target,
    tol = abs_tol,
    maxiter = maxiter
  )

  if (is.na(last_reported_progress) || last_reported_progress < 100) {
    p(step = 100)
  }

  results
}
