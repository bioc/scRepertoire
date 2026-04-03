#Making lodes to function in alluvial plots
#' @importFrom ggalluvial to_lodes_form
.makingLodes <- function(meta2, color, alpha, facet, set.axes) {
  diffuse_elements <- c()
  if (!is.null(color)) {
    diffuse_elements <- c(diffuse_elements, color)
  }
  if (!is.null(alpha)) {
    diffuse_elements <- c(diffuse_elements, alpha)
  }
  if (!is.null(facet)) {
    diffuse_elements <- c(diffuse_elements, facet)
  }

  lodes <- to_lodes_form(meta2, key = "x", value = "stratum",
                         id = "alluvium", axes = set.axes,
                         diffuse = if(length(diffuse_elements) > 0) diffuse_elements else NULL)
  return(lodes)
}

#' Alluvial Plotting for Single-Cell Object
#'
#' View the proportional contribution of clones by Seurat or SCE object
#' meta data after [combineExpression()]. The visualization
#' is based on the ggalluvial package, which requires the aesthetics
#' to be part of the axes that are visualized. Therefore, alpha, facet,
#' and color should be part of the the axes you wish to view or will
#' add an additional stratum/column to the end of the graph.
#'
#' @examples
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list,
#'                         samples = c("P17B", "P17L", "P18B", "P18L",
#'                                     "P19B","P19L", "P20B", "P20L"))
#'
#' # Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#'
#' # Using combineExpresion()
#' scRep_example <- combineExpression(combined, scRep_example)
#' scRep_example$Patient <- substring(scRep_example$orig.ident, 1,3)
#'
#' # Using alluvialClones()
#' alluvialClones(scRep_example,
#'                    clone.call = "gene",
#'                    y.axes = c("Patient", "ident"),
#'                    color = "ident")
#'
#' # Show only top 50 most frequent clones
#' alluvialClones(scRep_example,
#'                    clone.call = "aa",
#'                    y.axes = c("Patient", "ident"),
#'                    top.clones = 50)
#'
#' # Highlight specific clones
#' alluvialClones(scRep_example,
#'                    clone.call = "aa",
#'                    y.axes = c("Patient", "ident"),
#'                    highlight.clones = c("CVVSDNTGGFKTIF_CASSVRRERANTGELFF"))
#'
#' @param sc.data The product of [combineExpression()].
#' @param clone.call Defines the clonal sequence grouping. Accepted values
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL`,
#' `IGK`, `Light` (for both light chains), or `both` (for TRA/B and Heavy/Light).
#' @param y.axes The columns that will separate the proportional visualizations.
#' @param color The column header or clone(s) to be highlighted.
#' @param facet The column label to separate.
#' @param alpha The column header to have gradated opacity.
#' @param top.clones Show only the top N clones by frequency. If `NULL` (default),
#' show all clones.
#' @param min.freq Minimum frequency threshold for displaying flows. Clones
#' appearing fewer than this many times are filtered out.
#' @param highlight.clones Character vector of specific clone sequences to highlight.
#' These clones will be colored distinctly while others are shown in gray.
#' @param highlight.color Color to use for highlighted clones (default: "red").
#' @param stratum.width Width of the stratum bars (default: 0.2).
#' @param flow.alpha Transparency of the flows (default: 0.5). Highlighted clones
#' use full opacity.
#' @param show.labels If `TRUE` (default), display stratum labels.
#' @param label.size Text size for stratum labels (default: 2).
#' @param order.strata Named list specifying the order of levels within each stratum.
#' Names should match column names in y.axes.
#' @param export.table If `TRUE`, returns a data frame of the results
#' instead of a plot.
#' @param palette Colors to use in visualization - input any
#' [hcl.pals][grDevices::hcl.pals].
#' @param cloneCall \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `clone.call` instead.
#' @param exportTable \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `export.table` instead.
#' @param ... Additional arguments passed to the ggplot theme
#'
#' @importFrom ggalluvial StatStratum geom_flow geom_stratum to_lodes_form geom_alluvium
#'
#' @export
#' @concept SC_Functions
#' @return A ggplot object visualizing categorical distribution of clones, or a
#' data.frame if `export.table = TRUE`.
alluvialClones <- function(sc.data,
                           clone.call = NULL,
                           chain = "both",
                           y.axes = NULL,
                           color = NULL,
                           facet = NULL,
                           alpha = NULL,
                           top.clones = NULL,
                           min.freq = 0,
                           highlight.clones = NULL,
                           highlight.color = "red",
                           stratum.width = 0.2,
                           flow.alpha = 0.5,
                           show.labels = TRUE,
                           label.size = 2,
                           order.strata = NULL,
                           export.table = NULL,
                           palette = "inferno",
                           # Deprecated arguments
                           cloneCall = NULL,
                           exportTable = NULL,
                           ...) {

  # Handle deprecated arguments
  clone.call <- .deprecate_arg(cloneCall, clone.call, "cloneCall", "clone.call",
                               "alluvialClones", default = "strict")
  export.table <- .deprecate_arg(exportTable, export.table, "exportTable", "export.table",
                                 "alluvialClones", default = FALSE)

  x <- alluvium <- stratum <- NULL
  .checkSingleObject(sc.data)
  clone.call <- .theCall(.grabMeta(sc.data), clone.call)
  if (length(y.axes) == 0) {
    stop("Make sure you have selected the variable(s) to visualize")
  }
  meta <- .grabMeta(sc.data)
  if (chain != "both") {
    meta <- .offTheChain(meta, chain, clone.call)
  }
  meta$barcodes <- rownames(meta)
  meta <- meta[!is.na(meta[,clone.call]),]

  # Filter by minimum frequency
  if(min.freq > 0) {
    clone.counts <- table(meta[[clone.call]])
    keep.clones <- names(clone.counts)[clone.counts >= min.freq]
    meta <- meta[meta[[clone.call]] %in% keep.clones, ]
  }

  # Filter to top clones by frequency
  if(!is.null(top.clones)) {
    clone.counts <- sort(table(meta[[clone.call]]), decreasing = TRUE)
    keep.clones <- names(clone.counts)[seq_len(min(top.clones, length(clone.counts)))]
    meta <- meta[meta[[clone.call]] %in% keep.clones, ]
  }

  # Handle highlight.clones - creates a special highlighting column
  highlight.mode <- FALSE
  if(!is.null(highlight.clones)) {
    highlight.mode <- TRUE
    meta$`.highlight` <- ifelse(meta[[clone.call]] %in% highlight.clones,
                                "Highlighted", "Other")
    # Ensure highlighted clones are plotted last (on top)
    meta <- meta[order(meta$`.highlight` == "Highlighted"), ]
  }

  # Handle color parameter - check if it's clone sequences or a column
  check <- colnames(meta) == color
  if (length(unique(check)) == 1 & unique(check)[1] == FALSE &
      !is.null(color)) {
    meta <- meta %>% mutate("clone(s)" = ifelse(meta[,clone.call] %in%
                                                      color, "Selected", "Other"))
    color <- "clone(s)"
  }

  # Prepping the data for calculating lodes
  y.axes.plot <- unique(c(y.axes, color, alpha, facet))
  set.axes <- seq_along(y.axes.plot)
  cols.to.keep <- c(y.axes.plot, clone.call, "barcodes")
  if(highlight.mode) {
    cols.to.keep <- c(cols.to.keep, ".highlight")
  }
  meta2 <- meta[, cols.to.keep]
  meta2 <- unique(na.omit(meta2[!duplicated(as.list(meta2))]))

  # Apply stratum ordering if specified
  if(!is.null(order.strata)) {
    for(axis.name in names(order.strata)) {
      if(axis.name %in% colnames(meta2)) {
        meta2[[axis.name]] <- factor(meta2[[axis.name]],
                                     levels = order.strata[[axis.name]])
      }
    }
  }

  lodes <- .makingLodes(meta2, color, alpha, facet, set.axes)

  # Filtering the lodes
  if(any(lodes[,clone.call] != "")) {
    lodes <- lodes[lodes[,clone.call] != "",]
  }
  if(any(is.na(lodes[,clone.call]))) {
    lodes <- lodes[!is.na(lodes[,clone.call]),]
  }

  # Add frequency and proportion info for export
  if(export.table) {
    # Calculate frequency per clone
    clone.freq <- as.data.frame(table(meta[[clone.call]]))
    colnames(clone.freq) <- c(clone.call, "freq")
    lodes <- merge(lodes, clone.freq, by = clone.call, all.x = TRUE)

    # Calculate proportion
    total.cells <- nrow(meta)
    lodes$prop <- lodes$freq / total.cells

    # Add rank
    lodes$rank <- match(lodes[[clone.call]],
                        names(sort(table(meta[[clone.call]]), decreasing = TRUE)))
    return(lodes)
  }

  # Plotting
  plot <- ggplot(data = lodes, aes(x = x,
                                   stratum = stratum,
                                   alluvium = alluvium,
                                   label = stratum)) +
                geom_stratum(width = stratum.width, color = "black")

  # Handle highlighting mode
  if(highlight.mode) {
    # Use highlight column for coloring
    plot <- plot + geom_flow(aes(fill = lodes[[".highlight"]],
                                 alpha = lodes[[".highlight"]]),
                             stat = "alluvium",
                             lode.guidance = "forward",
                             width = stratum.width) +
                   scale_fill_manual(values = c("Highlighted" = highlight.color,
                                               "Other" = "grey70"),
                                    name = "Clone") +
                   scale_alpha_manual(values = c("Highlighted" = 1,
                                                "Other" = flow.alpha),
                                     guide = "none")
  } else if (is.null(color) & is.null(alpha)) {
    plot <- plot + geom_alluvium(width = stratum.width, alpha = flow.alpha)
  } else if (!is.null(color) & is.null(alpha)) {
    plot <- plot + geom_flow(aes(fill = lodes[,color]),
                             stat = "alluvium",
                             lode.guidance = "forward",
                             width = stratum.width,
                             alpha = flow.alpha) +
                   labs(fill = color)
  } else if (is.null(color) & !is.null(alpha)) {
    plot <- plot + geom_flow(aes(alpha = lodes[,alpha]),
                             stat = "alluvium",
                             lode.guidance = "forward",
                             width = stratum.width) +
                   labs(alpha = alpha)
  } else {
    plot <- plot + geom_flow(aes(alpha = lodes[,alpha],
                                 fill = lodes[,color]),
                             stat = "alluvium",
                             lode.guidance = "forward",
                             width = stratum.width) +
      labs(fill = color, alpha = alpha)
  }

  # Add faceting
  if (length(facet) == 1 & length(facet) < 2) {
    plot <- plot +
              facet_wrap(.~lodes[,facet], scales="free_y")
  }

  # Add labels
  if(show.labels) {
    plot <- plot +
            geom_text(stat = ggalluvial::StatStratum,
                      infer.label = FALSE,
                      reverse = TRUE,
                      size = label.size)
  }

  # Apply colors if not in highlight mode
  if(!highlight.mode && !is.null(color)) {
    n.colors <- length(unique(lodes[,color]))
    plot <- plot +
            scale_fill_manual(values = .colorizer(palette, n.colors))
  }

  plot <- plot +
            scale_x_discrete(expand = c(0.025, 0.025)) +
            .themeRepertoire(...) +
            theme(axis.title.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  line = element_blank())

  return(plot)
}
