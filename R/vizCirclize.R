#' Visualize Clonal Relationships as a Chord Diagram
#'
#' This function creates a chord diagram visualization of shared clones between
#' groups using the circlize package. It provides a convenient wrapper around
#' [getCirclize()] that handles the circlize plotting code automatically.
#'
#' @details
#' This function requires the circlize package to be installed. If circlize is

#' not available, the function will return the data that would be used for plotting
#' and provide instructions for manual plotting.
#'
#' The chord diagram shows relationships between groups (sectors) where the width
#' of each chord represents the number or proportion of shared clones between
#' the connected groups.
#'
#' @examples
#' \dontrun{
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list,
#'                         samples = c("P17B", "P17L", "P18B", "P18L",
#'                                     "P19B","P19L", "P20B", "P20L"))
#'
#' # Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#' scRep_example <- combineExpression(combined,
#'                                    scRep_example)
#'
#' # Simple chord diagram
#' vizCirclize(scRep_example, group.by = "seurat_clusters")
#'
#' # Directional chord diagram with arrows
#' vizCirclize(scRep_example,
#'             group.by = "seurat_clusters",
#'             directional = TRUE)
#'
#' # Multi-level grouping
#' scRep_example$Patient <- substring(scRep_example$orig.ident, 1, 3)
#' vizCirclize(scRep_example,
#'             group.by = c("Patient", "seurat_clusters"))
#' }
#'
#' @param sc.data The single-cell object after [combineExpression()].
#' @param clone.call Defines the clonal sequence grouping. Accepted values
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param group.by A column header (or vector of column headers for hierarchical
#' grouping) in the metadata to group the analysis by.
#' @param method The method for calculating link values: `"unique"` (default) counts
#' unique shared clones, `"jaccard"` calculates Jaccard similarity,
#' `"overlap"` calculates overlap coefficient.
#' @param proportion Calculate the relationship by unique clones (`FALSE`, default)
#' or normalized by proportion (`TRUE`).
#' @param directional If `TRUE`, show directional arrows on chords. Default is `FALSE`.
#' @param self.link How to handle self-links. `1` = show as loops, `2` = show as
#' parallel lines. Default is `1`.
#' @param include.self Include self-links (clones within a single group). Default is `TRUE`.
#' @param transparency Transparency of the chord links (0-1). Default is `0.5`.
#' @param link.visible If `FALSE`, hide the chord links and show only sectors.
#' @param annotate.sectors If `TRUE` (default), display sector names.
#' @param min.shared Minimum number of shared clones to include a link (default 0).
#' @param palette Colors to use for sectors - input any [hcl.pals][grDevices::hcl.pals].
#' @param sector.colors Named vector of colors for specific sectors. Overrides palette.
#' @param export.table If `TRUE`, returns the data instead of plotting.
#'
#' @export
#' @concept SC_Functions
#' @return Invisibly returns the circlize data list. If circlize is not installed
#' or `export.table = TRUE`, returns the data that would be used for plotting.
#' @seealso [getCirclize()] for generating the underlying data
#' @author Nick Borcherding
vizCirclize <- function(sc.data,
                        clone.call = NULL,
                        group.by = NULL,
                        method = c("unique", "jaccard", "overlap"),
                        proportion = FALSE,
                        directional = FALSE,
                        self.link = 1,
                        include.self = TRUE,
                        transparency = 0.5,
                        link.visible = TRUE,
                        annotate.sectors = TRUE,
                        min.shared = 0,
                        palette = "inferno",
                        sector.colors = NULL,
                        export.table = FALSE) {

  method <- match.arg(method)

  # Check for circlize package
  if(!requireNamespace("circlize", quietly = TRUE)) {
    message("The circlize package is required for visualization.")
    message("Install it with: install.packages('circlize')")
    message("Returning data for manual plotting instead.")
    export.table <- TRUE
  }

  # Get the circlize data
  circ.data <- getCirclize(
    sc.data = sc.data,
    clone.call = clone.call,
    group.by = group.by,
    method = method,
    proportion = proportion,
    symmetric = !directional,
    include.self = include.self,
    include.metadata = TRUE,
    min.shared = min.shared,
    palette = palette
  )

  # Apply custom sector colors if provided
  if(!is.null(sector.colors)) {
    for(s in names(sector.colors)) {
      if(s %in% names(circ.data$colors)) {
        circ.data$colors[s] <- sector.colors[s]
      }
    }
  }

  if(export.table) {
    return(circ.data)
  }

  # Check if there are links to plot
  if(nrow(circ.data$links) == 0) {
    warning("No links to plot after filtering. Try reducing min.shared or check your data.")
    return(invisible(circ.data))
  }

  # Clear any existing circlize plots
  circlize::circos.clear()

  # Set up circlize parameters
  circlize::circos.par(
    gap.after = 2,
    track.margin = c(0.01, 0.01)
  )

  # Determine annotation track setting
  annotation_track <- if(annotate.sectors) "name" else NULL

  # Create the chord diagram
  if(directional) {
    circlize::chordDiagram(
      circ.data$links,
      grid.col = circ.data$colors,
      transparency = transparency,
      self.link = self.link,
      link.visible = link.visible,
      directional = 1,
      direction.type = "arrows",
      link.arr.type = "big.arrow",
      annotationTrack = annotation_track
    )
  } else {
    circlize::chordDiagram(
      circ.data$links,
      grid.col = circ.data$colors,
      transparency = transparency,
      self.link = self.link,
      link.visible = link.visible,
      annotationTrack = annotation_track
    )
  }

  # Clear circos parameters at the end
  on.exit(circlize::circos.clear(), add = TRUE)

  # Return data invisibly
  invisible(circ.data)
}
