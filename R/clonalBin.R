#' Bin Clones by Frequency or Proportion
#'
#' This function adds a clonal grouping variable (`cloneSize`) to the output
#' of [combineTCR()], [combineBCR()], or [combineExpression()]. It calculates
#' the clonal frequency and proportion, then bins clones into categories based
#' on customizable thresholds. This is useful for categorizing clones prior to
#' downstream analysis or visualization.
#'
#' @examples
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list,
#'                         samples = c("P17B", "P17L", "P18B", "P18L",
#'                                     "P19B","P19L", "P20B", "P20L"))
#'
#' # Adding clonal bins with default settings (proportion-based)
#' combined <- clonalBin(combined)
#'
#' # Adding clonal bins based on frequency
#' combined <- clonalBin(combined,
#'                       proportion = FALSE,
#'                       clone.size = c(Rare = 1, Small = 5, Medium = 20,
#'                                      Large = 100, Hyperexpanded = 500))
#'
#' # Using a custom grouping variable
#' combined <- addVariable(combined,
#'                         variable.name = "Type",
#'                         variables = rep(c("B", "L"), 4))
#' combined <- clonalBin(combined, group.by = "Type")
#'
#' @param input.data The product of [combineTCR()], [combineBCR()], or
#' [combineExpression()].
#' @param clone.call Defines the clonal sequence grouping. Accepted values
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL`
#' (for both light chains), `both`.
#' @param group.by A column header in the metadata to group the analysis
#' by (e.g., "sample", "treatment"). If `NULL`, data will be analyzed
#' by list element.
#' @param proportion Whether to use proportion (`TRUE`) or total
#' frequency (`FALSE`) of the clone for binning.
#' @param clone.size The bins for the grouping based on proportion or frequency.
#' If proportion is `FALSE` and the clone.size values are not set high enough
#' based on frequency, the upper limit of clone.size will be automatically
#' updated.
#' @param cloneCall \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `clone.call` instead.
#' @param cloneSize \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `clone.size` instead.
#' @importFrom dplyr bind_rows group_by summarise mutate n
#' @importFrom rlang .data
#' @importFrom stats na.omit
#' @export
#' @concept Clonal_Analysis
#' @author Nick Borcherding
#' @return A list of data frames with clonal frequency, clonal proportion,
#' and cloneSize columns added.
#'
clonalBin <- function(input.data,
                      clone.call = NULL,
                      chain = "both",
                      group.by = NULL,
                      proportion = TRUE,
                      clone.size = NULL,
                      # Deprecated arguments
                      cloneCall = NULL,
                      cloneSize = NULL) {

  # Avoid R CMD check note for dplyr column references
  clonalFrequency <- NULL

  # Handle deprecated arguments
  clone.call <- .deprecate_arg(cloneCall, clone.call, "cloneCall", "clone.call",
                               "clonalBin", default = "strict")
  clone.size <- .deprecate_arg(cloneSize, clone.size, "cloneSize", "clone.size",
                               "clonalBin",
                               default = c(Rare = 1e-4, Small = 0.001,
                                          Medium = 0.01, Large = 0.1,
                                          Hyperexpanded = 1))

  # Suppress summarise messages
  options(dplyr.summarise.inform = FALSE)

 # Validate clone.size for frequency-based binning
  if (!proportion && any(clone.size < 1)) {
    stop("Adjust the clone.size parameter - there are groupings < 1")
  }

  # Add the baseline bin
  clone.size <- c(None = 0, clone.size)

  # Resolve clone.call to actual column name
  clone.call <- .theCall(input.data, clone.call)

  # Handle chain filtering if needed
  if (chain != "both") {
    for (i in seq_along(input.data)) {
      input.data[[i]] <- .offTheChain(input.data[[i]], chain, clone.call, check = FALSE)
    }
  }

  # Ensure input.data is a list
  input.data <- .checkList(input.data)

  # Remove any existing clonal columns to allow recalculation
  cols_to_remove <- c("clonalFrequency", "clonalProportion", "cloneSize")
  for (i in seq_along(input.data)) {
    existing_cols <- intersect(cols_to_remove, colnames(input.data[[i]]))
    if (length(existing_cols) > 0) {
      input.data[[i]] <- input.data[[i]][, !colnames(input.data[[i]]) %in% existing_cols, drop = FALSE]
    }
  }

  # Initialize output list
  output.data <- list()

  # Compute clonal frequency and proportion
  if (is.null(group.by) || group.by == "none") {
    # Calculate per list element
    for (i in seq_along(input.data)) {
      data <- data.frame(input.data[[i]], stringsAsFactors = FALSE)
      data2 <- unique(data[, c("barcode", clone.call)])
      data2 <- na.omit(data2)

      # Calculate frequency and proportion
      data2 <- data2 %>%
        group_by(data2[, clone.call]) %>%
        summarise(clonalProportion = dplyr::n() / nrow(data2),
                  clonalFrequency = dplyr::n())
      colnames(data2)[1] <- clone.call

      # Merge back to original data
      data <- merge(data, data2, by = clone.call, all = TRUE)
      output.data[[i]] <- data
    }
    names(output.data) <- names(input.data)

  } else {
    # Calculate across all elements grouped by specified variable
    data <- data.frame(bind_rows(input.data), stringsAsFactors = FALSE)
    data2 <- na.omit(unique(data[, c("barcode", clone.call, group.by)]))

    data2 <- data2 %>%
      group_by(.data[[clone.call]], .data[[group.by]]) %>%
      summarise(clonalFrequency = n(), .groups = "drop") %>%
      group_by(.data[[group.by]]) %>%
      mutate(clonalProportion = clonalFrequency / sum(clonalFrequency))

    colnames(data2)[c(1, 2)] <- c(clone.call, group.by)
    data <- merge(data, data2, by = c(clone.call, group.by), all = TRUE)

    # Split back into list by original list element names
    if ("element.names" %in% colnames(data)) {
      output.data <- split(data, data$element.names)
    } else {
      # If no element.names, re-split based on input names
      output.data <- list(data)
      names(output.data) <- "combined"
    }
  }

  # Combine all data for bin assignment to get consistent bin thresholds
  Con.df <- bind_rows(output.data)

  # Use shared helper to assign bins (returns both df and formatted clone.size)
  bin_result <- .assignCloneSizeBins(Con.df, clone.size, proportion)
  clone.size <- bin_result$clone.size

  # Apply binning to each element in output.data
  for (i in seq_along(output.data)) {
    element_result <- .assignCloneSizeBins(output.data[[i]], clone.size, proportion)
    output.data[[i]] <- element_result$df

    # Convert to factor with reversed levels (largest first)
    output.data[[i]]$cloneSize <- factor(
      output.data[[i]]$cloneSize,
      levels = rev(names(clone.size))
    )
  }

  return(output.data)
}
