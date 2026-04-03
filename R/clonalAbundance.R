#' Plot the Relative Abundance of Clones
#'
#' Displays the number of clones at specific frequencies by sample
#' or group. Visualization can either be a line graph (
#' `scale` = FALSE) using calculated numbers or density
#' plot (`scale` = TRUE). Multiple sequencing runs can
#' be group together using the group parameter. If a matrix
#' output for the data is preferred, set
#' `export.table` = TRUE.
#'
#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list,
#'                         samples = c("P17B", "P17L", "P18B", "P18L",
#'                                     "P19B","P19L", "P20B", "P20L"))
#'
#' # Using clonalAbundance()
#' clonalAbundance(combined,
#'                 clone.call = "gene",
#'                 scale = FALSE)
#'
#' @param input.data The product of [combineTCR()],
#' [combineBCR()], or [combineExpression()].
#' @param clone.call Defines the clonal sequence grouping. Accepted values
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL`,
#' `IGK`, `Light` (for both light chains), or `both` (for TRA/B and Heavy/Light).
#' @param group.by A column header in the metadata or lists to group the analysis
#' by (e.g., "sample", "treatment"). If `NULL`, data will be analyzed
#' by list element or active identity in the case of single-cell objects.
#' @param order.by A character vector defining the desired order of elements
#' of the `group.by` variable. Alternatively, use `alphanumeric` to sort groups
#' automatically.
#' @param scale Converts the graphs into density plots in order to show
#' relative distributions.
#' @param export.table If `TRUE`, returns a data frame or matrix of the results
#' instead of a plot.
#' @param palette Colors to use in visualization - input any
#' [hcl.pals][grDevices::hcl.pals].
#' @param cloneCall \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `clone.call` instead.
#' @param exportTable \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `export.table` instead.
#' @param ... Additional arguments passed to the ggplot theme
#'
#' @importFrom ggplot2 ggplot
#' @importFrom dplyr distinct
#' @export
#' @concept Visualizing_Clones
#' @author Nick Borcherding, Justin Reimertz
#' @return A ggplot object visualizing clonal abundance by group, or a
#' data.frame if `export.table = TRUE`.
clonalAbundance <- function(input.data,
                            clone.call = NULL,
                            chain = "both",
                            scale = FALSE,
                            group.by = NULL,
                            order.by = NULL,
                            export.table = NULL,
                            palette = "inferno",
                            # Deprecated arguments
                            cloneCall = NULL,
                            exportTable = NULL,
                            ...) {

  # Handle deprecated arguments
  clone.call <- .deprecate_arg(cloneCall, clone.call, "cloneCall", "clone.call",
                               "clonalAbundance", default = "strict")
  export.table <- .deprecate_arg(exportTable, export.table, "exportTable", "export.table",
                                 "clonalAbundance", default = FALSE)
  Con.df <- NULL
  xlab <- "Abundance"
  input.data <- .dataWrangle(input.data,
                             group.by,
                             .theCall(input.data, clone.call,
                                      check.df = FALSE, silent = TRUE),
                             chain)
  clone.call <- .theCall(input.data, clone.call)
  
  names <- names(input.data)
  if (!is.null(group.by)) {
  	Con.df <- bind_rows(lapply(seq_along(input.data), function(x) {
  		input.data %>%
  			.parseContigs(x, names, clone.call) %>%
  			mutate("{group.by}" := input.data[[x]][[group.by]])
  	})) %>% 
  		dplyr::distinct() %>%
  	  as.data.frame()
  	col <- length(unique(Con.df[[group.by]]))
    fill <- group.by
    if(!is.null(order.by)) {
        Con.df <- .orderingFunction(vector = order.by,
                                     group.by = group.by, 
                                     data.frame =  Con.df)
    }
    if (scale == TRUE) { 
        ylab <- "Density of Clones"
        plot <- ggplot(Con.df, aes(x=Abundance, 
                                   fill=Con.df[,group.by])) +
                      geom_density(aes(y=after_stat(scaled)), 
                                   alpha=0.5, 
                                   lwd=0.25, 
                                   color="black", 
                                   bw=0.5)  +
                      scale_fill_manual(values = .colorizer(palette,col)) +
                      labs(fill = fill)
    } else { 
        ylab <- "Number of Clones"
        plot <- ggplot(Con.df, aes(x=Abundance, group.by = values, 
                               color = Con.df[,group.by])) +
                        geom_line(stat="count") +
                        scale_color_manual(values = .colorizer(palette,col)) +
                        labs(color = fill)
    }
  } else {
    for (i in seq_along(input.data)) {
      data1 <- .parseContigs(input.data, i, names, clone.call) 
      Con.df<- rbind.data.frame(Con.df, data1) %>%
                  dplyr::distinct() %>%
                  as.data.frame()
    }
    Con.df <- data.frame(Con.df)
    if(!is.null(order.by)) {
      Con.df <- .orderingFunction(vector = order.by,
                                  group.by = "values", 
                                  data.frame = Con.df)
    }
    
    col <- length(unique(Con.df$values))
    fill <- "Samples"
    if (scale == TRUE) { 
      ylab <- "Density of Clones"
      plot <- ggplot(Con.df, aes(Abundance, fill=values)) +
                      geom_density(aes(y=after_stat(scaled)), 
                                   alpha=0.5, 
                                   lwd=0.25, 
                                   color="black", 
                                   bw=0.5) +
                      scale_fill_manual(values = .colorizer(palette,col)) +
                      labs(fill = fill)
    } else { 
      ylab <- "Number of Clones"
      plot <- ggplot(Con.df, aes(x=Abundance, group = values, 
                               color = values)) +
                      geom_line(stat="count") +
                      scale_color_manual(values = .colorizer(palette,col)) +
                      labs(color = fill)
  } }
  if (export.table == TRUE) { 
    return(Con.df) 
  }
  plot <- plot + scale_x_log10() + ylab(ylab) + xlab(xlab) +
    .themeRepertoire(...)
  return(plot) 
}
