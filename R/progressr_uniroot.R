uniroot_with_progress <- function(
  fun,
  target,
  bounds,
  abs_tol = 1e-6,
  maxiter = 100
) {
  state <- new.env(parent = emptyenv())
  state$last_val <- NA
  state$first_interval <- NA
  state$last_reported_progress <- 0

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
    if (!is.na(state$last_val)) {
      if (est - state$last_val == 0) {
        progress <- 1
      } else {
        conv_crit <- .Machine$double.eps * abs(target_val) + abs_tol / 2
        interval <- abs(state$last_val - target_val)
        if (is.na(state$first_interval)) {
          state$first_interval <- interval
        }
        progress <- (log(interval) - log(state$first_interval)) /
          (log(conv_crit) - log(state$first_interval))
      }
      logger::log_debug(
        "Current estimate: ", est,
        ", error: ", error,
        ", progress: ", progress
      )
      progress_percent <- min(ceiling(100 * progress), 100)
      add_progress <- progress_percent - state$last_reported_progress
      if (add_progress > 0) {
        message <- paste(
          "Est.: ", round(est, 4), "(", round(target_val, 4),
          ")"
        )
        p(message, amount = add_progress)
        state$last_reported_progress <- progress_percent +
          state$last_reported_progress
      }
    }
    state$last_val <- target_val
    error
  }

  results <- stats::uniroot(optimization_goal,
    interval = bounds,
    target = target,
    tol = abs_tol,
    maxiter = maxiter
  )

  last_p <- state$last_reported_progress
  if (is.na(last_p) || last_p < 100) {
    p(step = 100)
  }

  results
}
