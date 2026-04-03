#' Highlighting Specific Clones
#'
#' Use a specific clonal sequence to highlight on top of the dimensional 
#' reduction in single-cell object.
#'
#' @examples
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' # Getting a sample of a Seurat object
#' scRep_example  <- get(data("scRep_example"))
#' 
#' # Using combineExpresion()
#' scRep_example  <- combineExpression(combined, 
#'                                     scRep_example)
#' 
#' # Using highlightClones()
#' scRep_example   <- highlightClones(scRep_example,
#'                                    clone.call= "aa",
#'                                    sequence = c("CVVSDNTGGFKTIF_CASSVRRERANTGELFF"))
#' 
#' @param sc.data The single-cell object to attach after
#' [combineExpression()]
#' @param clone.call Defines the clonal sequence grouping. Accepted values
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param sequence The specific sequence or sequence to highlight
#' @param cloneCall \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `clone.call` instead.
#' @importFrom S4Vectors DataFrame
#' @export
#' @concept SC_Functions
#' @return Single-cell object object with new meta data column
#' for indicated clones
highlightClones <- function(sc.data,
                            clone.call = NULL,
                            sequence = NULL,
                            # Deprecated arguments
                            cloneCall = NULL){

  # Handle deprecated arguments
  clone.call <- .deprecate_arg(cloneCall, clone.call, "cloneCall", "clone.call",
                               "highlightClones", default = "strict")

  if (!.is.seurat.or.se.object(sc.data)) {
    stop("Please select a single-cell object")
  }

  clone.call <- .theCall(.grabMeta(sc.data), clone.call)
  meta <- .grabMeta(sc.data)
  meta$highlight <- NA
  for(i in seq_along(sequence)) {
    meta$highlight <-  ifelse(meta[,clone.call] == sequence[i],
                              sequence[i], meta$highlight)
  }
  
  meta <- meta[,-(which(colnames(meta) == "ident"))]
  
  if(.is.se.object(sc.data)) {
    colData(sc.data) <- DataFrame(meta)
  } else {
    col.name <- names(meta) %||% colnames(meta)
    sc.data[[col.name]] <- meta
  }
  return(sc.data)
}
