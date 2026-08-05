#' Prepare data and phylogeny for analyses (mandatory step)
#'
#' Aligns species-level invasion data with a phylogeny and stores
#' predictor information for downstream analyses.
#'
#' If no phylogeny is supplied, a tree is retrieved from the Open Tree
#' of Life using \pkg{rotl}; Grafen branch lengths are added if branch
#' lengths are absent. We recommend providing a high-quality tree if available,
#' and also recommend fully matching the species in the data and phylogeny.
#'
#' @param df Data frame containing at least a \code{Species} column formatted as
#' "Genus_species", and an \code{Invasive} column (0: not invasive, 1: invasive).
#' @param tree Optional object of class \code{"phylo"}.
#' @param predictors Optional character vector of column names (additional
#' predictors beyond phylogeny; if NULL, models will be based on the tree alone).
#' @param species_col Name of species column (default "Species" is recommended).
#' @param rho Grafen branch length scaling parameter (default 0.5).
#' @param plot Logical. If TRUE, visualizes the phylogeny in \code{tree}.
#'
#' @return An object of class \code{"invasible_prepared"}.
#'
#' @examples
#' species_list <- fish_beginning_with_E
#' prep <- prepare_invasible(species_list,rho=1,predictors=
#'         c("Fake_categorical_trait","Fake_continuous_trait"))
#' prep$tree
#'
#' @export
prepare_invasible <- function(
    df,
    tree = NULL,
    predictors = NULL,
    species_col = "Species",
    rho = 0.5,
    plot = FALSE
) {

  # --------------------------------------------------
  # Basic validation
  # --------------------------------------------------

  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame.")
  }

  if (!(species_col %in% names(df))) {
    stop("Column '", species_col, "' not found in df.")
  }

  if (!("Invasive" %in% names(df))) {
    stop("Column 'Invasive' not found in df.")
  }

  if (anyDuplicated(df[[species_col]]) > 0) {

    dup <- unique(
      df[[species_col]][duplicated(df[[species_col]])]
    )

    stop(
      "Duplicate species names detected: ",
      paste(dup, collapse = ", ")
    )
  }

  # --------------------------------------------------
  # Retrieve tree from Open Tree of Life
  # --------------------------------------------------

  if (is.null(tree)) {

    if (!requireNamespace("rotl", quietly = TRUE)) {
      stop(
        "Package 'rotl' is required when tree = NULL."
      )
    }

    tol <- rotl::tnrs_match_names(df[[species_col]])

    # Multiple names resolving to same OTT taxon
    dup_ott <- unique(
      tol$ott_id[duplicated(tol$ott_id)]
    )

    if (length(dup_ott) > 0) {

      offenders <- tol[
        tol$ott_id %in% dup_ott,
        c("search_string", "unique_name", "ott_id")
      ]

      msg <- paste(
        apply(
          offenders,
          1,
          function(x) {
            paste(
              x["search_string"],
              "->",
              x["unique_name"],
              "(OTT",
              x["ott_id"],
              ")"
            )
          }
        ),
        collapse = "\n"
      )

      stop(
        "Multiple submitted names resolved to the same OTT taxon:\n",
        msg
      )
    }

    tree <- suppressWarnings(rotl::tol_induced_subtree(
      ott_ids = tol$ott_id
    ))
    tree$tip.label <- sub("_ott[0-9]+$", "", tree$tip.label)
  }
  tree$node.label <- NULL

  # --------------------------------------------------
  # Validate tree
  # --------------------------------------------------

  if (!inherits(tree, "phylo")) {
    stop("'tree' must be an object of class 'phylo'.")
  }

  if (anyDuplicated(tree$tip.label) > 0) {

    dup <- unique(
      tree$tip.label[
        duplicated(tree$tip.label)
      ]
    )

    stop(
      "Tree contains duplicated tip labels: ",
      paste(dup, collapse = ", ")
    )
  }

  # --------------------------------------------------
  # Add branch lengths if missing
  # --------------------------------------------------

  if (is.null(tree$edge.length)) {

    tree <- ape::compute.brlen(
      tree,
      method = "Grafen",
      power = rho
    )
  }

  # --------------------------------------------------
  # Check overlap
  # --------------------------------------------------

  missing_from_data <- setdiff(
    tree$tip.label,
    df[[species_col]]
  )

  if (length(missing_from_data) > 0) {

    stop(
      "Species present in tree but absent from data: ",
      paste(missing_from_data, collapse = ", ")
    )
  }

  missing_from_tree <- setdiff(
    df[[species_col]],
    tree$tip.label
  )

  if (length(missing_from_tree) > 0) {

    warning(
      "Removing species absent from tree: ",
      paste(missing_from_tree, collapse = ", ")
    )

    df <- df[
      df[[species_col]] %in% tree$tip.label,
      ,
      drop = FALSE
    ]
  }

  # --------------------------------------------------
  # Align data to tree
  # --------------------------------------------------

  df <- df[
    match(tree$tip.label, df[[species_col]]),
    ,
    drop = FALSE
  ]

  rownames(df) <- df[[species_col]]

  # --------------------------------------------------
  # Validate predictors
  # --------------------------------------------------

  if (!is.null(predictors)) {

    missing_pred <- setdiff(
      predictors,
      names(df)
    )

    if (length(missing_pred) > 0) {

      stop(
        "Predictor columns not found: ",
        paste(missing_pred, collapse = ", ")
      )
    }
  }
  if (plot) {
    if (!requireNamespace("phytools", quietly = TRUE)) {
      stop("Package 'phytools' required for plotting.")
    }
    tre <- phytools::plotTree(tree,type="fan",lwd=0.5,fsize=0.5)
  }

  # --------------------------------------------------
  # Return object
  # --------------------------------------------------

  structure(
    list(
      data = df,
      tree = tree,
      predictors = predictors
    ),
    class = "invasible_prepared"
  )
}
