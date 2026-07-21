#' Predict invasion risk and evaluate model performance
#'
#' 1. Fits a phylogenetic generalized linear mixed model (PGLMM) from the
#' output of \code{invasion_signal()}
#' 2. Calculates the probability of being invasive for all species,
#' based on the phylogeny (and optionally species traits) given to the
#' original \code{prepare_invasible()} function
#' 3. Computes predictive performance metrics such as ROC AUC
#' 4. Plots the predicted values for observed invasive species (1) vs
#' non-invading species (0) (optional)
#'
#' @param signal Output of \code{invasion_signal()}
#' @param plot Logical. If TRUE, returns diagnostic plot of observed vs predicted
#' @param vcv_alpha Numeric, optional: set a value of alpha for OU model which is
#' different to the value inherited from the \code{invasion_signal()} output
#'
#' @return An object of class \code{pred_output}, a list containing:
#' \describe{
#'   \item{model}{Fitted PGLMM model object.}
#'   \item{predictions}{Data frame with observed and predicted values for all species.}
#'   \item{ranked_predictions}{Candidate invasives ordered by decreasing probability.}
#'   \item{roc}{ROC object}
#'   \item{auc}{Numeric AUC (area under the curve) value.}
#'   \item{plot}{ggplot object (if requested; run pred_output$plot to show).}
#' }
#'
#' @examples
#' species_list <- fish_beginning_with_E
#' prep <- prepare_invasible(species_list,rho=1, predictors=c("Fake_continuous_trait"))
#' signal <- invasion_signal(prep)
#' pred <- predict_invasible(signal,plot=TRUE)
#' pred$ranked_predictions
#' pred$auc
#' pred$plot
#'
#' @export
predict_invasible <- function(signal,vcv_alpha=NULL,plot=FALSE) {

  if (!inherits(signal, "invasible_signal")) {
    stop("Input must be from invasion_signal().")
  }

  if (isTRUE(signal$p_random > 0.05)) {
    warning("Phylogenetic signal (D) is not significant. Interpret results with caution!")
  }

  if (!requireNamespace("phyr", quietly = TRUE)) stop("Package 'phyr' required.")
  if (!requireNamespace("phytools", quietly = TRUE)) stop("Package 'phytools' required.")
  if (!requireNamespace("pROC", quietly = TRUE)) stop("Package 'pROC' required.")

  # ---------------------------
  # core objects
  # ---------------------------
  df <- signal$prepared$data
  tree <- signal$tree
  predictors <- signal$predictors

  df$Species <- as.factor(df$Species)
  df$Invasive <- as.numeric(df$Invasive)

  if (!inherits(tree, "phylo")) {
    stop("Tree must be of class 'phylo'.")
  }

  # ---------------------------
  # phylogenetic covariance
  # ---------------------------
  if (is.null(vcv_alpha)) {
    V <- phytools::vcvPhylo(tree, model = "OU", alpha = signal$alpha)
  } else {
    V <- phytools::vcvPhylo(tree, model = "OU", alpha = vcv_alpha)
    }
  # ---------------------------
  # formula construction
  # ---------------------------
  if (is.null(predictors) || length(predictors) == 0) {
    fixed_formula <- "1"
  } else {
    predictors <- predictors[predictors %in% names(df)]
    fixed_formula <- if (length(predictors) == 0) "1" else paste(predictors, collapse = " + ")
  }

  fml <- stats::as.formula(
    paste("Invasive ~", fixed_formula, "+ (1 | Species__)")
  )

  # ---------------------------
  # PGLMM fit
  # ---------------------------
  fit <- phyr::pglmm(
    formula = fml,
    data = df,
    family = "gaussian",
    cov_ranef = list(Species = V),
    REML = FALSE
  )

  preds <- phyr::pglmm_predicted_values(fit)

  pred_df <- data.frame(
    Species = df$Species,
    Observed = df$Invasive,
    Predicted = preds$Y_hat
  )

  # ---------------------------
  # ROC + AUC
  # ---------------------------
  roc_obj <- pROC::roc(
    response = pred_df$Observed,
    predictor = pred_df$Predicted,
    quiet = TRUE
  )

  auc_value <- as.numeric(pROC::auc(roc_obj))

  # ---------------------------
  # ranking
  # ---------------------------
  df_ranked <- pred_df[order(-pred_df$Predicted), ]
  df_ranked$rank <- seq_len(nrow(df_ranked))
  df_ranked <- subset(df_ranked,Observed==0)
  # ---------------------------
  # optional plot
  # ---------------------------
  plot_obj <- NULL

  if (plot) {

        if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' required for plotting.")
    }

    plot_obj <- ggplot2::ggplot(
      pred_df,
      ggplot2::aes(
        x = as.factor(Observed),
        y = Predicted,
        color = as.factor(Observed)
      )
    ) +
      ggplot2::geom_jitter(alpha = 0.3, size = 0.9, width = 0.35) +
      ggplot2::stat_summary(
        fun = mean,
        geom = "point",
        size = 3
      ) +
      ggplot2::stat_summary(
        fun.data = ggplot2::mean_sdl,
        fun.args = list(mult = 1),
        geom = "errorbar",
        width = 0.75
      ) +
      ggplot2::scale_color_manual(values = c("0" = "deepskyblue3", "1" = "red")) +
      ggplot2::theme_classic()
      plot_obj
  }

  # ---------------------------
  # return object
  # ---------------------------
  structure(
    list(
      model = fit,
      predictions = pred_df,
      ranked_predictions = df_ranked,
      roc = roc_obj,
      auc = auc_value,
      plot = plot_obj
    ),
    class = "pred_output"
  )
}
