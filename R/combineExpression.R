#' Adding Clonal Information to Single-Cell Object
#'
#' This function adds the immune receptor information to the Seurat or 
#' SCE object to the meta data. By default this function also calculates 
#' the frequencies and proportion of the clones by sequencing 
#' run (`group.by` = NULL). To change how the frequencies/proportions
#' are calculated, select a column header for the `group.by` variable. 
#' Importantly, before using [combineExpression()] ensure the 
#' barcodes of the single-cell object object match the barcodes in the output 
#' of the [combineTCR()] or [combineBCR()]. 
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
#' 
#' @param input.data The product of [combineTCR()], [combineBCR()] or a list of
#' both c([combineTCR()], [combineBCR()]).
#' @param sc.data The Seurat or Single-Cell Experiment (SCE) object to attach
#' @param clone.call Defines the clonal sequence grouping. Accepted values
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL`
#' (for both light chains), `both`.
#' @param group.by A column header in lists to group the analysis
#' by (e.g., "sample", "treatment"). If `NULL`, will be based on the list element.
#' @param proportion Whether to proportion (`TRUE`) or total
#' frequency (`FALSE`) of the clone based on the group.by variable.
#' @param clone.size The bins for the grouping based on proportion or frequency.
#' If proportion is `FALSE` and the clone.sizes are not set high enough
#' based on frequency, the upper limit of clone.sizes will be automatically
#' updated.
#' @param filter.na Method to subset Seurat/SCE object of barcodes without
#' clone information
#' @param add.label This will add a label to the frequency header, allowing
#' the user to try multiple group.by variables or recalculate frequencies after
#' subsetting the data.
#' @param cloneCall \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `clone.call` instead.
#' @param cloneSize \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `clone.size` instead.
#' @param filterNA \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `filter.na` instead.
#' @param addLabel \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `add.label` instead.
#' @importFrom dplyr left_join all_of coalesce
#' @importFrom  rlang %||% sym :=
#' @importFrom SummarizedExperiment colData<- colData
#' @importFrom S4Vectors DataFrame
#' @export
#' @concept SC_Functions
#' @return Single-cell object with clone information added to meta data
#' information
#'
combineExpression <- function(input.data,
                              sc.data,
                              clone.call = NULL,
                              chain = "both",
                              group.by = NULL,
                              proportion = TRUE,
                              filter.na = NULL,
                              clone.size = NULL,
                              add.label = NULL,
                              # Deprecated arguments
                              cloneCall = NULL,
                              cloneSize = NULL,
                              filterNA = NULL,
                              addLabel = NULL) {

  # Handle deprecated arguments
  clone.call <- .deprecate_arg(cloneCall, clone.call, "cloneCall", "clone.call",
                               "combineExpression", default = "strict")
  clone.size <- .deprecate_arg(cloneSize, clone.size, "cloneSize", "clone.size",
                               "combineExpression", default = c(Rare = 1e-4, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = 1))
  filter.na <- .deprecate_arg(filterNA, filter.na, "filterNA", "filter.na",
                              "combineExpression", default = FALSE)
  add.label <- .deprecate_arg(addLabel, add.label, "addLabel", "add.label",
                              "combineExpression", default = FALSE)

  clonalFrequency <- NULL
  call_time <- Sys.time()
    options( dplyr.summarise.inform = FALSE )
    if (!proportion && any(clone.size < 1)) {
        stop("Adjust the clone.size parameter - there are groupings < 1")
    }
    clone.size <- c(None = 0, clone.size)

    clone.call <- .theCall(input.data, clone.call)
    if (chain != "both") {
      #Retain the full clone information
      full.clone <- lapply(input.data, function(x) {
        x[, c("barcode", clone.call)]
      })
      full.clone <- bind_rows(full.clone)
      for(i in seq_along(input.data)) {
        input.data[[i]] <- .offTheChain(input.data[[i]], chain, clone.call, check = FALSE)
      }
    }
    input.data <- .checkList(input.data)

    # Remove any existing clonal columns to allow recalculation
    cols_to_remove <- c("clonalFrequency", "clonalProportion", "cloneSize")
    for (i in seq_along(input.data)) {
      existing_cols <- intersect(cols_to_remove, colnames(input.data[[i]]))
      if (length(existing_cols) > 0) {
        input.data[[i]] <- input.data[[i]][, !colnames(input.data[[i]]) %in% existing_cols, drop = FALSE]
      }
    }

    #Getting Summaries of clones from combineTCR() or combineBCR()
    Con.df <- NULL
    meta <- .grabMeta(sc.data)
    cell.names <- rownames(meta)

    conDfColnamesNoCloneSize <- unique(c(
        "barcode", CT_lines, clone.call, "clonalProportion", "clonalFrequency"
    ))

    # Carry any retained full-length sequence columns (from combineTCR/BCR with
    # retain.sequences) into the single-cell metadata. Additive and present-only:
    # when none are retained this is a no-op and the default path is unchanged.
    seq_candidates <- paste0(rep(.seqColMap, each = 2), c("1", "2"))
    retained_present <- Reduce(intersect,
        c(list(seq_candidates), lapply(input.data, colnames)))
    if (length(retained_present) > 0) {
        conDfColnamesNoCloneSize <- unique(c(conDfColnamesNoCloneSize, retained_present))
    }

    # Computes the clonalProportion and clonalFrequency for each clone
    if (is.null(group.by) || group.by == "none") {

        for (i in seq_along(input.data)) {

            data <- data.frame(input.data[[i]], stringsAsFactors = FALSE)
            data2 <- unique(data[,c("barcode", clone.call)])
            #This ensures all calculations are based on the cells in the SCO
            data2 <- na.omit(data2[data2[,"barcode"] %in% cell.names,])
            data2 <- data2 %>%
                        group_by(data2[,clone.call]) %>%
                        summarise(clonalProportion = dplyr::n()/nrow(data2),
                                  clonalFrequency = dplyr::n())
            colnames(data2)[1] <- clone.call
            data <- merge(data, data2, by = clone.call, all = TRUE)
            data <- data[, conDfColnamesNoCloneSize]
            Con.df <- rbind.data.frame(Con.df, data)
        }

    } else {
        data <- data.frame(bind_rows(input.data), stringsAsFactors = FALSE)
        data2 <- na.omit(unique(data[,c("barcode", clone.call, group.by)]))
        #This ensures all calculations are based on the cells in the SCO
        data2 <- data2[data2[,"barcode"] %in% cell.names, ]
        data2 <- data2 %>%
          group_by(.data[[clone.call]], .data[[group.by]]) %>%
          summarise(clonalFrequency = n(), .groups = "drop") %>%
          group_by(.data[[group.by]]) %>%
          mutate(clonalProportion = clonalFrequency / sum(clonalFrequency))

        colnames(data2)[c(1,2)] <- c(clone.call, group.by)
        data <- merge(data, data2, by = c(clone.call, group.by), all = TRUE)
        Con.df <- data[, conDfColnamesNoCloneSize]
    }

    #Use shared helper to assign cloneSize bins
    bin_result <- .assignCloneSizeBins(Con.df, clone.size, proportion)
    Con.df <- bin_result$df
    clone.size <- bin_result$clone.size

    #Formating the meta data to add and removing any duplicate barcodes
    PreMeta <- unique(Con.df[, c(conDfColnamesNoCloneSize, "cloneSize")])
    dup <- PreMeta$barcode[which(duplicated(PreMeta$barcode))]
    PreMeta <- PreMeta[!PreMeta$barcode %in% dup,]

    #Re-adding full clones
    if (chain != "both") {
      clone_sym <- sym(clone.call)
      PreMeta <- PreMeta %>%
        left_join(full.clone, by = "barcode", suffix = c("", ".from_full_clones")) %>%
        mutate(!!clone_sym := coalesce(!!sym(paste0(clone.call, ".from_full_clones")), !!clone_sym)) %>%
        dplyr::select(-all_of(paste0(clone.call, ".from_full_clones")))
    }
    barcodes <- PreMeta$barcode
    PreMeta <- PreMeta[,-1]
    rownames(PreMeta) <- barcodes
    if (!is.null(group.by) && group.by != "none" && add.label) {
      location <- which(colnames(PreMeta) %in% c("clonalProportion",
                          "clonalFrequency"))
      colnames(PreMeta)[location] <- paste0(c("clonalProportion",
                                            "clonalFrequency"), group.by)
    }
    
    if (.is.seurat.object(sc.data)) { 
        if (length(which(rownames(PreMeta) %in% 
                         rownames(sc.data[[]])))/length(rownames(sc.data[[]])) < 0.01) {
          getHighBarcodeMismatchError()
        }
        col.name <- names(PreMeta) %||% colnames(PreMeta)
        sc.data[[col.name]] <- PreMeta
    } else {
      rownames <- rownames(colData(sc.data))
      if (length(which(rownames(PreMeta) %in% 
                       rownames))/length(rownames) < 0.01) {
        getHighBarcodeMismatchError() }
      
      combined_col_names <- unique(c(colnames(colData(sc.data)), colnames(PreMeta)))
      full_data <- merge(colData(sc.data), PreMeta[rownames, , drop = FALSE], by = "row.names", all.x = TRUE)
      # at this point, the rows in full_data are shuffled. match back with the original colData
      full_data <- full_data[match(rownames, full_data[,1]), ]
      rownames(full_data) <- full_data[, 1]
      full_data  <- full_data[, -1]
      colData(sc.data) <- DataFrame(full_data[, combined_col_names])  
    }
    if (filter.na) {
      sc.data <- .filteringNA(sc.data)
    }
    sc.data$cloneSize <- factor(sc.data$cloneSize, levels = rev(names(clone.size)))
    
    if(.is.seurat.object(sc.data)) {
        sc.data@commands[["combineExpression"]] <- .makeScrepSeurat(
              call_time, sc.data@active.assay)
    }
    return(sc.data)
}

getHighBarcodeMismatchError <- function() {
  stop(
    paste0(
      "< 1% of barcodes match.\n",
      "Ensure the barcodes in the single-cell object match the barcodes from scRepertoire.\n",
      "For help, see: https://www.borch.dev/uploads/screpertoire/articles/faq"
    ),
    call. = FALSE 
  )
}
