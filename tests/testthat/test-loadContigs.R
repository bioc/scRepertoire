# test script for loadContigs.R - testcases are NOT comprehensive!

check_loadContigs_output <- function(loaded_data) {
  # Check if the output is a list containing a single data frame
  expect_type(loaded_data, "list")
  expect_length(loaded_data, 1)
  df <- loaded_data[[1]]
  expect_s3_class(df, "data.frame")
  
  # Check if the data frame is not empty
  expect_gt(nrow(df), 0)
  
  # Check for the presence of essential standardized columns
  expected_cols <- c("barcode", "chain", "reads", "v_gene", 
                     "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3")
  expect_true(all(expected_cols %in% names(df)))
  
  # Check data types of key columns
  expect_type(df$chain, "character")
  expect_type(df$cdr3, "character")
  # Reads should be numeric/integer after parsing
  expect_true(is.numeric(df$reads) || is.integer(df$reads))
}


# Minimal Cell Ranger-style 10X contig table that includes the framework/CDR
# region columns recent Cell Ranger versions emit (older versions lack these).
make_10x_with_regions <- function() {
  data.frame(
    barcode    = c("cellA", "cellB"),
    is_cell    = TRUE,
    contig_id  = c("cellA_contig_1", "cellB_contig_1"),
    high_confidence = TRUE,
    length     = 500L,
    chain      = c("TRA", "TRB"),
    v_gene     = c("TRAV1", "TRBV2"),
    d_gene     = c("None", "TRBD1"),
    j_gene     = c("TRAJ1", "TRBJ1"),
    c_gene     = c("TRAC", "TRBC1"),
    full_length = TRUE,
    productive = TRUE,
    fwr1       = c("FAA", "FBB"),       fwr1_nt = c("TTTGCAGCA", "TTTGCTGCT"),
    cdr1       = c("CA", "CB"),         cdr1_nt = c("TGTGCA", "TGTGCT"),
    fwr2       = c("FA2", "FB2"),       fwr2_nt = c("TTTGCAGCT", "TTTGCTGCC"),
    cdr2       = c("DA", "DB"),         cdr2_nt = c("GATGCA", "GATGCT"),
    fwr3       = c("FA3", "FB3"),       fwr3_nt = c("TTTGCAGCG", "TTTGCTGCG"),
    cdr3       = c("CASSA", "CASSB"),   cdr3_nt = c("TGTGCAAGCAGCGCA", "TGTGCAAGCAGCGCT"),
    fwr4       = c("FA4", "FB4"),       fwr4_nt = c("TTTGGAGGA", "TTTGGTGGT"),
    reads      = c(100L, 120L),
    umis       = c(10L, 12L),
    raw_clonotype_id  = c("clone1", "clone2"),
    raw_consensus_id  = c("consensus1", "consensus2"),
    stringsAsFactors = FALSE
  )
}

test_that("loadContigs reconstructs full-length sequence from 10X region columns", {
  df <- make_10x_with_regions()
  loaded <- loadContigs(df, format = "10X")[[1]]

  expect_true(all(c("sequence", "sequence_aa") %in% colnames(loaded)))

  # Re-derive expected full-length sequences for the row matching each barcode.
  exp_nt <- paste0(df$fwr1_nt, df$cdr1_nt, df$fwr2_nt, df$cdr2_nt,
                   df$fwr3_nt, df$cdr3_nt, df$fwr4_nt)
  exp_aa <- paste0(df$fwr1, df$cdr1, df$fwr2, df$cdr2,
                   df$fwr3, df$cdr3, df$fwr4)
  names(exp_nt) <- df$barcode
  names(exp_aa) <- df$barcode

  expect_identical(loaded$sequence,    unname(exp_nt[loaded$barcode]))
  expect_identical(loaded$sequence_aa, unname(exp_aa[loaded$barcode]))
})

test_that("loadContigs adds NA sequence columns for formats without full-length data", {
  WAT3R <- read.csv("https://www.borch.dev/uploads/contigs/WAT3R_contigs.csv")
  loaded <- loadContigs(WAT3R, format = "WAT3R")[[1]]
  expect_true(all(c("sequence", "sequence_aa") %in% colnames(loaded)))
  expect_true(all(is.na(loaded$sequence)))
  expect_true(all(is.na(loaded$sequence_aa)))
})

test_that("loadContigs retains native sequence columns from AIRR-family input", {
  airr <- data.frame(
    cell_id = c("cellA", "cellB"),
    locus   = c("TRA", "TRB"),
    consensus_count = c(5L, 7L),
    v_call  = c("TRAV1", "TRBV2"),
    d_call  = c("TRAD1", "TRBD1"),
    j_call  = c("TRAJ1", "TRBJ1"),
    c_call  = c("TRAC", "TRBC1"),
    junction = c("ATGCGT", "ATGCGA"),
    junction_aa = c("ML", "MR"),
    sequence = c("AAAATGCGTAAA", "CCCATGCGACCC"),
    sequence_aa = c("KMRK", "PMRP"),
    sequence_alignment = c("AAA...ATGCGTAAA", "CCC...ATGCGACCC"),
    germline_alignment = c("AAA...ATGCGTAAA", "CCC...ATGCGGCCC"),
    stringsAsFactors = FALSE
  )
  loaded <- loadContigs(airr, format = "AIRR")[[1]]
  expect_true(all(c("sequence", "sequence_aa",
                    "sequence_alignment", "germline_alignment") %in% colnames(loaded)))
  # Match by barcode since .order_df may reorder rows.
  idx <- match(loaded$barcode, airr$cell_id)
  expect_identical(loaded$sequence, airr$sequence[idx])
  expect_identical(loaded$germline_alignment, airr$germline_alignment[idx])
})

test_that("loadContigs correctly processes various formats from URL", {
  #TRUST4 format
  TRUST4 <- read.csv("https://www.borch.dev/uploads/contigs/TRUST4_contigs.csv")
  trial_trust4 <- loadContigs(TRUST4, format = "TRUST4")
  check_loadContigs_output(trial_trust4)
  
  # BD format
  BD <- read.csv("https://www.borch.dev/uploads/contigs/BD_contigs.csv")
  trial_bd <- loadContigs(BD, format = "BD")
  check_loadContigs_output(trial_bd)
  
  # WAT3R format
  WAT3R <- read.csv("https://www.borch.dev/uploads/contigs/WAT3R_contigs.csv")
  trial_wat3r <- loadContigs(WAT3R, format = "WAT3R")
  check_loadContigs_output(trial_wat3r)
  
  # 10X format (using pre-loaded contig_list data)
  data("contig_list")
  trial_10x <- loadContigs(contig_list[[1]], format = "10X")
  check_loadContigs_output(trial_10x)
  
  # MiXCR format
  MIXCR <- read.csv("https://www.borch.dev/uploads/contigs/MIXCR_contigs.csv")
  trial_mixcr <- loadContigs(MIXCR, format = "MiXCR")
  check_loadContigs_output(trial_mixcr)
  
  # Immcantation format
  Immcantation <- read.csv("https://www.borch.dev/uploads/contigs/Immcantation_contigs.csv")
  trial_immcantation <- loadContigs(Immcantation, format = "Immcantation")
  check_loadContigs_output(trial_immcantation)
  
  # ParseBio format
  Parse <- read.csv("https://www.borch.dev/uploads/contigs/Parse_contigs.csv")
  trial_parse <- loadContigs(Parse, format = "ParseBio")
  check_loadContigs_output(trial_parse)
  
  # Dandelion format
  Dandelion <- read.csv("https://www.borch.dev/uploads/contigs/Dandelion_contigs.csv")
  trial_dandelion <- loadContigs(Dandelion, format = "Dandelion")
  check_loadContigs_output(trial_dandelion)
})

test_that("loadContigs correctly auto-detects and processes various formats", {
  
  BD <- read.csv("https://www.borch.dev/uploads/contigs/BD_contigs.csv")
  trial_bd <- expect_message(
    loadContigs(BD, format = "auto"), 
    "Automatically detected format from data: AIRR"
  )
  
  WAT3R <- read.csv("https://www.borch.dev/uploads/contigs/WAT3R_contigs.csv")
  trial_wat3r <- expect_message(
    loadContigs(WAT3R, format = "auto"), 
    "Automatically detected format from data: WAT3R"
  )
  
  trial_10x <- expect_message(
      loadContigs(contig_list[[1]], format = "auto"), 
      "Automatically detected format from data: 10X"
  )
  
  MIXCR <- read.csv("https://www.borch.dev/uploads/contigs/MIXCR_contigs.csv")
  trial_mixcr <- expect_message(
    loadContigs(MIXCR, format = "auto"), 
    "Automatically detected format from data: MiXCR"
  )
  
  Immcantation <- read.csv("https://www.borch.dev/uploads/contigs/Immcantation_contigs.csv")
  trial_immcantation <- expect_message(
    loadContigs(Immcantation, format = "auto"), 
    "Automatically detected format from data: Immcantation"
  )
  
  Parse <- read.csv("https://www.borch.dev/uploads/contigs/Parse_contigs.csv")
  trial_parse <- expect_message(
    loadContigs(Parse, format = "auto"), 
    "Automatically detected format from data: ParseBio"
  )
})

test_that("loadContigs works with AIRR input (directory mode)", {
  # This test remains unchanged as it is self-contained and does not use internal data.
  # Create a temporary TSV file for AIRR format
  tmp_tsv <- tempfile(fileext = ".tsv")
  airr_data <- data.frame(
    cell_id = "cellA",
    locus = "TRB",
    consensus_count = 5,
    v_call = "TRBV1",
    d_call = "TRBD1",
    j_call = "TRBJ1",
    c_call = "TRBC1",
    junction = "ATGCGT",
    junction_aa = "ML",
    stringsAsFactors = FALSE
  )
  write.table(airr_data, tmp_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Copy to a temporary directory with proper filename
  tmp_dir_airr <- file.path(tempdir(), "airr_test")
  dir.create(tmp_dir_airr, showWarnings = FALSE)
  file.copy(tmp_tsv, file.path(tmp_dir_airr, "airr_rearrangement.tsv"), overwrite = TRUE)
  
  # Explicit format test
  result_airr <- loadContigs(tmp_dir_airr, format = "AIRR")
  expected_airr <- list(
    data.frame(
      barcode = "cellA",
      chain = "TRB",
      reads = as.integer(5),
      v_gene = "TRBV1",
      d_gene = "TRBD1",
      j_gene = "TRBJ1",
      c_gene = "TRBC1",
      cdr3_nt = "ATGCGT",
      cdr3 = "ML",
      sequence = NA_character_,
      sequence_aa = NA_character_,
      stringsAsFactors = FALSE
    )
  )
  expect_identical(result_airr, expected_airr)
})

test_that("loadContigs returns empty list when no matching files are found in a directory", {
  empty_dir <- file.path(tempdir(), "empty_test_dir")
  dir.create(empty_dir, showWarnings = FALSE)
  file.remove(list.files(empty_dir, full.names = TRUE))
  expect_warning(res <- loadContigs(empty_dir, format = "10X"))
  expect_identical(res, list())
})
