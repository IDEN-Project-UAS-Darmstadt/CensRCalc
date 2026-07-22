partial <- function(f, fixed_args) {
  force(f)
  force(fixed_args)
  new_fun <- function(...) {
    args <- list(...)
    do.call(f, c(args, fixed_args))
  }
  class(new_fun) <- c("censrcalc_partial", class(new_fun))
  new_fun
}

#' @importFrom methods formalArgs
partial_formalargs <- function(x) {
  if (inherits(x, "censrcalc_partial")) {
    env <- environment(x)
    f <- get("f", envir = env, inherits = TRUE)
    fixed_args <- get("fixed_args", envir = env, inherits = TRUE)
    setdiff(
      setdiff(formalArgs(f), names(fixed_args)),
      names(formals(x))
    )
  } else {
    formalArgs(x)
  }
}
