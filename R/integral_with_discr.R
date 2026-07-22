integral_with_discr <- function(
  fun,
  bounds,
  abs_tol = 1e-6,
  rel_tol = abs_tol, # We assume values around 1 bec of probabilities/densities
  vectorize = FALSE
) {
  coll <- checkmate::makeAssertCollection()
  checkmate::assert_list(
    bounds,
    types = c("numeric", "list"),
    any.missing = FALSE,
    min.len = 1,
    names = "unique",
    add = coll
  )

  args <- partial_formalargs(fun)

  checkmate::assert_function(fun, add = coll)
  if (!"..." %in% args) {
    checkmate::assert_names(
      names(bounds),
      type = "unique",
      permutation.of = args,
      add = coll
    )
    checkmate::assert_names(
      args,
      type = "unique",
      permutation.of = names(bounds),
      add = coll,
      .var.name = "function arguments"
    )
  }
  checkmate::assert_numeric(
    abs_tol,
    len = 1,
    any.missing = FALSE,
    finite = TRUE,
    lower = 0,
    add = coll
  )
  checkmate::assert_numeric(
    rel_tol,
    len = 1,
    any.missing = FALSE,
    finite = TRUE,
    lower = 0,
    add = coll
  )
  checkmate::assert_flag(vectorize, add = coll)
  checkmate::reportAssertions(coll)

  is_cont <- vapply(bounds, is.numeric, logical(1))
  continuous <- bounds[is_cont]
  checkmate::qassertr(continuous, "N2", .var.name = "bounds (continuous part)")
  discrete <- bounds[!is_cont]
  checkmate::qassertr(discrete, "L+", .var.name = "bounds (discrete part)")

  int_res <- 0
  int_objects <- list()

  integrate_once <- function(fun, bounds) {
    int <- calculus::integral(
      fun,
      bounds = bounds,
      absTol = abs_tol,
      relTol = rel_tol,
      vectorize = vectorize
    )
    if (checkmate::test_numeric(
      int,
      len = 1,
      any.missing = FALSE,
      finite = TRUE
    )) {
      int <- list(
        value = int,
        error = 0,
        message = "No integration",
        bounds = bounds
      )
    }
    logger::log_trace(
      "Integral result: ",
      int$value,
      " with error ",
      int$error
    )
    int
  }

  if (length(discrete) == 0 && length(continuous) > 0) {
    int <- integrate_once(fun, continuous)
    int_res <- int$value
    int_objects <- list(int)
  } else {
    discrete_args_transformed <- list()
    for (discrete_arg_name in names(discrete)) {
      to_transform <- as.numeric(unlist(discrete[[discrete_arg_name]]))
      checkmate::assert_numeric(
        to_transform,
        any.missing = FALSE,
        len = length(to_transform),
        .var.name = paste0(
          "bounds for discrete variable ",
          discrete_arg_name
        )
      )
      discrete_args_transformed[[discrete_arg_name]] <- to_transform
    }
    # This is a little bit inefficient, but for now ok
    cart_product <- expand.grid(
      discrete_args_transformed,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    n_ints <- nrow(cart_product)

    int_objects <- vector("list", n_ints)
    for (i in seq_len(n_ints)) {
      args <- stats::setNames(as.list(cart_product[i, ]), names(cart_product))
      fun_part <- partial(fun, args)

      if (length(continuous) == 0) {
        int <- list(
          value = fun_part(),
          error = 0,
          message = "No integration, only discrete"
        )
        int_res <- int_res + int$value
        int_objects[[i]] <- int
      } else {
        int <- integrate_once(fun_part, continuous)
        int_res <- int_res + int$value
        int_objects[[i]] <- int
      }
    }
  }
  errors <- vapply(int_objects, `[[`, numeric(1), "error")
  logger::log_trace(
    "Integral result: ",
    int_res,
    ", errors: ",
    paste(errors, collapse = ", ")
  )
  list(
    value = int_res,
    error = sum(errors),
    int_objects = int_objects,
    errors = errors
  )
}
