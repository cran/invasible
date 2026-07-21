#' Simulated monitoring performance for predicted invaders
#'
#' Evaluates the probability of :
#' 1. Detecting at least one invasion within n given species, under ranked
#' monitoring, i.e. the highest-probability n species, vs a random selection.
#' 2. Given an invasion, detecting it within the n  species, under ranked
#' monitoring, i.e. the highest-probability n species, vs a random selection.
#'
#' @param pred_out Output from \code{predict_invasible()}.
#' @param plot Logical. If TRUE, returns ggplot objects.
#' @param max_n Integer. Highest n value to be plotted (n=25 as default).
#'
#' @return a list containing:
#' \describe{
#'   \item{results}{simulation results for each value of n}
#'   \item{plot1}{1st plot of \code{results} for each value until \code{max_n}}
#'   \item{plot2}{2nd plot of \code{results} for each value until \code{max_n}}
#' }
#'dev
#' @examples
#' species_list <- fish_beginning_with_E
#' prep <- prepare_invasible(species_list,rho=1, predictors=c("Fake_continuous_trait"))
#' signal <- invasion_signal(prep)
#' pred <- predict_invasible(signal)
#' set.seed(10)
#' probs <- monitor_species(pred,plot=TRUE,max_n=16)
#' probs$plot1
#' probs$plot2
#'
#' @export
monitor_species <- function(pred_out, plot = FALSE,max_n=25) {

  if (!inherits(pred_out, "pred_output")) {
    stop("Input must be output of predict_invasible().")
  }

  df <- pred_out$ranked_predictions

  if (!all(c("Observed", "Predicted") %in% names(df))) {
    stop("predictions must contain Observed and Predicted columns.")
  }

  df$Observed <- as.numeric(df$Observed)
  df$Predicted <- as.numeric(df$Predicted)

  # ---------------------------
  # filter candidate invaders
  # ---------------------------
  df2 <- df[df$Observed == 0, ]
  df2 <- df2[order(-df2$Predicted), ]

  p <- df2$Predicted
  N <- length(p)

  if (N == 0) {
    stop("No non-invader species (Observed == 0) found.")
  }

  total_p <- sum(p)

  # random permutation baseline
  df_shuffled <- df2[sample(nrow(df2)), ]
  p2 <- df_shuffled$Predicted

  # ---------------------------
  # compute curves
  # ---------------------------
  results <- data.frame(
    n = seq_len(N),

    p_at_least_one_top_n = sapply(seq_len(N), function(n) {
      1 - prod(1 - p[1:n])
    }),

    p_at_least_one_random_n = sapply(seq_len(N), function(n) {
      1 - prod(1 - p2[1:n])
    }),

    p_invasion_outside_top_n = sapply(seq_len(N), function(n) {
      1 - sum(p[(n + 1):N]) / total_p
    }),

    p_invasion_outside_random_n = sapply(seq_len(N), function(n) {
      1 - sum(p2[(n + 1):N]) / total_p
    })
  )

  # ---------------------------
  # optional plots
  # ---------------------------
  plot1 <- NULL
  plot2 <- NULL
  message <- NULL
  if (plot) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("ggplot2 required for plotting.")
    }

    # ---- Plot 1: outside invasion probability
    plot_df1 <- data.frame(
      n = rep(results$n, 2),
      Probability = c(results$p_invasion_outside_random_n,
                      results$p_invasion_outside_top_n),
      Species_set = rep(c("Random", "Top-ranked"), each = N)
    )

    plot1 <- ggplot2::ggplot(plot_df1,
                             ggplot2::aes(x = n, y = Probability, color = Species_set)) +
      ggplot2::geom_line(linewidth = 0.75) +
      ggplot2::scale_color_manual(values = c("Random" = "deepskyblue",
                                             "Top-ranked" = "red")) +
      ggplot2::labs(
        x = "Number of species in monitored set",
        y = "Given 1 invasion, probability it occurs in species set"
      ) +  ggplot2::theme_classic()+ggplot2::xlim(1,max_n)

    # ---- Plot 2: at least one invasion probability
    plot_df2 <- data.frame(
      n = rep(results$n, 2),
      Probability = c(results$p_at_least_one_random_n,
                      results$p_at_least_one_top_n),
      Species_set = rep(c("Random", "Top-ranked"), each = N)
    )

   message <- message("Note that there are TWO plots")

    plot2 <- ggplot2::ggplot(plot_df2,
                             ggplot2::aes(x = n, y = Probability, color = Species_set)) +
      ggplot2::geom_line(linewidth = 0.75) +
      ggplot2::scale_color_manual(values = c("Random" = "deepskyblue",
                                             "Top-ranked" = "red")) +
      ggplot2::labs(
        x = "Number of species in monitored set",
        y = "Probability of 1 or more invasion in species set"
      ) +
      ggplot2::theme_classic()+ggplot2::xlim(1,max_n)
  }
  # ---------------------------
  # return
  # ---------------------------
  list(
    results = results,
    plot1 = plot1,
    plot2 = plot2,
    message = message
  )
}
