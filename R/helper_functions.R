#' @noRd
make_formula <- function(response = "Invasive", predictors = NULL) {

  # intercept-only model if no predictors
  if (is.null(predictors) || length(predictors) == 0) {
    return(stats::reformulate("1", response = response))
  }
  # ensure predictors are character vector
  if (is.data.frame(predictors)) {
    stop("predictors must be a character vector, not a data frame")
  }
  stats::reformulate(predictors, response = response)
}

