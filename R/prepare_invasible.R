#' Prepare data and phylogeny for analyses (mandatory step)
#'
#' Aligns species-level invasion data with a phylogeny and stores
#' predictor information for downstream analyses.
#'
#' If no phylogeny is supplied, a tree is retrieved from the Open Tree
#' of Life using \pkg{rotl}. Unmatched species names are removed with a
#' warning. When multiple submitted names resolve to the same Open Tree
#' taxon, the first occurrence is retained and subsequent occurrences
#' are removed with a warning.
#'
#' If a supplied phylogeny has species that are absent from the data,
#' those tips are removed. Species present in the data but absent from
#' the phylogeny are removed from the data with a warning.
#'
#' If branch lengths are absent, Grafen branch lengths are added.
#'
#' @param df Data frame containing at least a \code{Species} column formatted
#'   as "Genus_species", and an \code{Invasive} column (0: not invasive,
#'   1: invasive).
#' @param tree Optional object of class \code{"phylo"}. If \code{NULL},
#'   a phylogeny is retrieved from the Open Tree of Life.
#' @param predictors Optional character vector of column names containing
#'   additional predictors. If \code{NULL}, no additional predictors are
#'   stored.
#' @param species_col Character string giving the name of the species column.
#'   Defaults to \code{"Species"}.
#' @param rho Numeric value controlling the power parameter used for Grafen
#'   branch lengths. Defaults to \code{0.5}.
#' @param plot Logical. If \code{TRUE}, plots the resulting phylogeny as a
#'   fan tree using \pkg{phytools}.
#'
#' @return An object of class \code{"invasible_prepared"} containing:
#'   \itemize{
#'     \item \code{data}: the aligned data frame;
#'     \item \code{tree}: the aligned phylogeny;
#'     \item \code{predictors}: the supplied predictor names.
#'   }
#'
#' @examples
#' \dontrun{
#' prep <- prepare_invasible(
#'   fish_beginning_with_E,
#'   rho = 1,
#'   predictors = c(
#'     "Fake_categorical_trait",
#'     "Fake_continuous_trait"
#'   )
#' )
#'
#' prep$tree
#' }
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

  if (!is.character(species_col) ||
      length(species_col) != 1L ||
      is.na(species_col)) {
    stop("'species_col' must be a single character string.")
  }

  if (!(species_col %in% names(df))) {
    stop(
      "Column '", species_col, "' not found in df."
    )
  }

  if (!("Invasive" %in% names(df))) {
    stop("Column 'Invasive' not found in df.")
  }

  if (!is.numeric(rho) ||
      length(rho) != 1L ||
      is.na(rho) ||
      !is.finite(rho) ||
      rho <= 0) {
    stop("'rho' must be a single positive finite numeric value.")
  }

  if (!is.logical(plot) ||
      length(plot) != 1L ||
      is.na(plot)) {
    stop("'plot' must be TRUE or FALSE.")
  }

  if (!is.null(predictors)) {

    if (!is.character(predictors)) {
      stop("'predictors' must be a character vector or NULL.")
    }

    if (anyNA(predictors)) {
      stop("'predictors' cannot contain NA values.")
    }

    if (anyDuplicated(predictors) > 0L) {
      predictors <- unique(predictors)
    }
  }

  if (nrow(df) == 0L) {
    stop("'df' contains no rows.")
  }

  # --------------------------------------------------
  # Validate species names
  # --------------------------------------------------

  if (anyNA(df[[species_col]]) ||
      any(!nzchar(trimws(as.character(df[[species_col]]))))) {

    bad <- which(
      is.na(df[[species_col]]) |
        !nzchar(trimws(as.character(df[[species_col]])))
    )

    warning(
      length(bad),
      " rows contain missing or empty species names ",
      "and will be removed."
    )

    df <- df[-bad, , drop = FALSE]
  }

  if (nrow(df) == 0L) {
    stop("No rows remain after removing missing species names.")
  }

  # Convert species names to character explicitly
  df[[species_col]] <- as.character(df[[species_col]])

  # --------------------------------------------------
  # Handle exact duplicate species names
  # --------------------------------------------------

  if (anyDuplicated(df[[species_col]]) > 0L) {

    dup <- unique(
      df[[species_col]][
        duplicated(df[[species_col]])
      ]
    )

    warning(
      "Duplicate species names detected. ",
      "Keeping the first occurrence and removing duplicates:\n",
      paste(dup, collapse = ", ")
    )

    df <- df[
      !duplicated(df[[species_col]]),
      ,
      drop = FALSE
    ]
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

    # OpenTree generally performs better with spaces than underscores.
    submitted_names <- df[[species_col]]

    search_names <- gsub(
      "_",
      " ",
      submitted_names,
      fixed = TRUE
    )

    # Match names to Open Tree of Life
    tol <- suppressWarnings(
      rotl::tnrs_match_names(search_names)
    )

    if (nrow(tol) == 0L) {
      stop(
        "No species could be matched to the Open Tree of Life."
      )
    }

    # ------------------------------------------------
    # Handle unmatched names
    # ------------------------------------------------

    unmatched <- is.na(tol$ott_id)

    if (any(unmatched)) {

      unmatched_names <- submitted_names[unmatched]

      warning(
        length(unmatched_names),
        " species could not be matched to the Open Tree ",
        "of Life and will be removed:\n",
        paste(unmatched_names, collapse = ", ")
      )

      df <- df[!unmatched, , drop = FALSE]

      tol <- tol[!unmatched, , drop = FALSE]

      submitted_names <- submitted_names[!unmatched]
    }

    if (nrow(df) == 0L || nrow(tol) == 0L) {
      stop(
        "No species remained after removing names that could ",
        "not be matched to the Open Tree of Life."
      )
    }

    # ------------------------------------------------
    # Handle multiple names mapping to same OTT taxon
    # ------------------------------------------------

    duplicated_ott <- duplicated(tol$ott_id) |
      duplicated(tol$ott_id, fromLast = TRUE)

    if (any(duplicated_ott)) {

      offenders <- data.frame(
        submitted = submitted_names[duplicated_ott],
        matched = tol$unique_name[duplicated_ott],
        ott_id = tol$ott_id[duplicated_ott],
        stringsAsFactors = FALSE
      )

      warning(
        "Multiple submitted names resolve to the same ",
        "Open Tree taxon. Keeping the first occurrence ",
        "and removing subsequent occurrences:\n",
        paste(
          apply(
            offenders,
            1L,
            function(x) {
              paste(
                x["submitted"],
                "->",
                x["matched"],
                "(OTT",
                x["ott_id"],
                ")"
              )
            }
          ),
          collapse = "\n"
        )
      )

      # Keep only first occurrence of each OTT taxon
      keep <- !duplicated(tol$ott_id)

      df <- df[keep, , drop = FALSE]
      tol <- tol[keep, , drop = FALSE]
    }

    if (nrow(df) == 0L || nrow(tol) == 0L) {
      stop(
        "No species remained after resolving Open Tree ",
        "taxonomic matches."
      )
    }

    # ------------------------------------------------
    # Construct induced subtree
    # ------------------------------------------------

    tree <- suppressWarnings(
      rotl::tol_induced_subtree(
        ott_ids = tol$ott_id
      )
    )

    if (!inherits(tree, "phylo")) {
      stop(
        "Open Tree of Life did not return a valid phylogenetic tree."
      )
    }

    # Remove OpenTree OTT suffixes
    tree$tip.label <- sub(
      "_ott[0-9]+$",
      "",
      tree$tip.label
    )

    # Convert tree labels to Genus_species format
    tree$tip.label <- gsub(
      " ",
      "_",
      tree$tip.label,
      fixed = TRUE
    )
  }

  # --------------------------------------------------
  # Validate supplied/generated tree
  # --------------------------------------------------

  if (!inherits(tree, "phylo")) {
    stop("'tree' must be an object of class 'phylo'.")
  }

  if (length(tree$tip.label) == 0L) {
    stop("The phylogeny contains no tips.")
  }

  if (anyNA(tree$tip.label) ||
      any(!nzchar(trimws(tree$tip.label)))) {
    stop(
      "The phylogeny contains missing or empty tip labels."
    )
  }

  tree$tip.label <- as.character(tree$tip.label)

  # --------------------------------------------------
  # Handle duplicated tree tip labels
  # --------------------------------------------------

  if (anyDuplicated(tree$tip.label) > 0L) {

    dup <- unique(
      tree$tip.label[
        duplicated(tree$tip.label)
      ]
    )

    stop(
      "Tree contains duplicated tip labels: ",
      paste(dup, collapse = ", "),
      ". A phylogeny with duplicated species labels ",
      "cannot be aligned unambiguously."
    )
  }

  # Remove node labels
  tree$node.label <- NULL

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
  # Standardize tree labels
  # --------------------------------------------------

  tree$tip.label <- gsub(
    " ",
    "_",
    tree$tip.label,
    fixed = TRUE
  )

  # --------------------------------------------------
  # Match data to tree
  # --------------------------------------------------

  missing_from_tree <- setdiff(
    df[[species_col]],
    tree$tip.label
  )

  if (length(missing_from_tree) > 0L) {

    warning(
      length(missing_from_tree),
      " species present in the data are absent from ",
      "the tree and will be removed:\n",
      paste(missing_from_tree, collapse = ", ")
    )

    df <- df[
      df[[species_col]] %in% tree$tip.label,
      ,
      drop = FALSE
    ]
  }

  # Remove tree tips not represented in the data
  extra_tree_tips <- setdiff(
    tree$tip.label,
    df[[species_col]]
  )

  if (length(extra_tree_tips) > 0L) {

    tree <- ape::drop.tip(
      tree,
      extra_tree_tips
    )
  }

  # --------------------------------------------------
  # Final overlap check
  # --------------------------------------------------

  if (nrow(df) == 0L) {
    stop(
      "No species remain after matching the data to the tree."
    )
  }

  if (ape::Ntip(tree) == 0L) {
    stop(
      "No species remain in the phylogeny after matching to the data."
    )
  }

  # Check that data and tree now match exactly
  if (!setequal(df[[species_col]], tree$tip.label)) {

    missing_after <- setdiff(
      df[[species_col]],
      tree$tip.label
    )

    extra_after <- setdiff(
      tree$tip.label,
      df[[species_col]]
    )

    msg <- c()

    if (length(missing_after) > 0L) {
      msg <- c(
        msg,
        paste(
          "Data species absent from tree:",
          paste(missing_after, collapse = ", ")
        )
      )
    }

    if (length(extra_after) > 0L) {
      msg <- c(
        msg,
        paste(
          "Tree species absent from data:",
          paste(extra_after, collapse = ", ")
        )
      )
    }

    stop(
      "Data and tree could not be aligned:\n",
      paste(msg, collapse = "\n")
    )
  }

  # --------------------------------------------------
  # Align data order to tree order
  # --------------------------------------------------

  match_idx <- match(
    tree$tip.label,
    df[[species_col]]
  )

  if (anyNA(match_idx)) {
    stop(
      "Internal error while aligning data to the phylogeny."
    )
  }

  df <- df[
    match_idx,
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

    if (length(missing_pred) > 0L) {
      stop(
        "Predictor columns not found: ",
        paste(missing_pred, collapse = ", ")
      )
    }
  }

  # --------------------------------------------------
  # Optional plotting
  # --------------------------------------------------

  if (plot) {

    if (!requireNamespace("phytools", quietly = TRUE)) {
      stop(
        "Package 'phytools' is required for plotting."
      )
    }

    phytools::plotTree(
      tree,
      type = "fan",
      lwd = 0.5,
      fsize = 0.5
    )
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
