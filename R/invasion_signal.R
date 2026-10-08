#' Compute phylogenetic signal and baseline invasion model
#'
#' This function must be run before all prediction/evaluation steps.
#'
#' 1. Computes Fritz & Purvis D statistic
#' 2. Fits phylogenetic logistic regression
#' 3. Extracts phylogenetic dependence parameter (alpha)
#' 4. Returns an "invasible_signal" object for prediction functions
#'
#' @param prepared Output of prepare_invasible()
#' #'
#' @return An object of class \code{invasible_signal}, a list containing:
#' \describe{
#'   \item{phylo_model}{ }
#'   \item{predictors}{ }
#'   \item{prepared}{ }
#'   \item{tree}{ }
#'   \item{D}{ }
#'   \item{p_random}{ },
#'   \item{p_brownian}{ }
#'   \item{alpha}{ }
#' }
#'
#' @examples
#' species_list <- fish_beginning_with_E
#' prep <- prepare_invasible(species_list)
#' signal <- invasion_signal(prep)
#' signal$D
#' signal$p_random
#'
#' @export
invasion_signal <- function(prepared) {

  if (!inherits(prepared, "invasible_prepared")) {
    stop("Input must be output of prepare_invasible().")
  }

  df <- prepared$data
  tree <- prepared$tree
  predictors <- prepared$predictors

  if (!requireNamespace("caper", quietly = TRUE)) stop("caper required.")
  if (!requireNamespace("phylolm", quietly = TRUE)) stop("phylolm required.")

  # ---- 1. Phylogenetic signal (D)
  comp <- caper::comparative.data(tree, df, names.col = "Species")

  D_res <- caper::phylo.d(
    comp,
    binvar = Invasive,
    permut = 10000
  )

  # ---- 2. Phylogenetic logistic regression
  fml <- make_formula("Invasive", predictors)

  phylo_fit <- phylolm::phyloglm(
    formula = fml,
    data = comp$data,
    phy = comp$phy,
    method = "logistic_MPLE"
  )

  structure(
    class = "invasible_signal",
    list(
            phylo_model = phylo_fit,
      predictors = predictors,
     prepared = prepared,
      tree = prepared$tree,
    D = D_res$DEstimate,
      p_random = D_res$Pval1,
      p_brownian = D_res$Pval0,
      alpha = phylo_fit$alpha
          )
  )
}
