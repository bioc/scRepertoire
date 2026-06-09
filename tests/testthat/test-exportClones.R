# test script for eportClones.R - testcases are NOT comprehensive!

# Setup: Create a combined object for testing using internal package data
data("contig_list", package = "scRepertoire")
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))

test_that("exportClones handles invalid format", {
  expect_error(
    exportClones(combined, format = "invalid_format"),
    "'arg' should be one of"
  )
})

test_that("'paired' format export works correctly", {
  # Test returning a data frame
  paired_df <- exportClones(combined, format = "paired", write.file = FALSE)
  expect_s3_class(paired_df, "data.frame")
  expect_named(paired_df, c("chain1_aa", "chain1_nt", "chain1_genes",
                            "chain2_aa", "chain2_nt", "chain2_genes", "group"))
  expect_true(nrow(paired_df) > 0)
  
  # Test writing to file
  temp_dir <- tempdir()
  file_path <- file.path(temp_dir, "paired_test.csv")
  exportClones(combined, format = "paired", dir = temp_dir, file.name = "paired_test.csv")
  expect_true(file.exists(file_path))
  
  # Verify file content
  read_df <- read.csv(file_path, row.names = 1) # row.names are barcodes
  expect_s3_class(read_df, "data.frame")
  expect_equal(ncol(read_df), ncol(paired_df))
  unlink(file_path) # Clean up
})

test_that("'TCRMatch' format export works correctly", {
  tcr_match_df <- exportClones(combined, format = "TCRMatch", write.file = FALSE)
  expect_s3_class(tcr_match_df, "data.frame")
  expect_named(tcr_match_df, c("chain2_aa", "group", "clonalFrequency"))
  expect_true(nrow(tcr_match_df) > 0)
})

test_that("'tcrpheno' format export works correctly", {
  tcr_pheno_df <- exportClones(combined, format = "tcrpheno", write.file = FALSE)
  expect_s3_class(tcr_pheno_df, "data.frame")
  expect_named(tcr_pheno_df, c("cell", "TCRA_cdr3aa", "TCRA_vgene", "TCRA_jgene",
                               "TCRA_cdr3nt", "TCRB_cdr3aa", "TCRB_vgene",
                               "TCRB_jgene", "TCRB_cdr3nt"))
  expect_true(nrow(tcr_pheno_df) > 0)
})

test_that("'airr' format export works correctly", {
  airr_df <- exportClones(combined, format = "airr", write.file = FALSE)
  expect_s3_class(airr_df, "data.frame")
  expect_named(airr_df, c("cell_id", "locus", "v_call", "d_call", "j_call",
                          "c_call", "junction", "junction_aa"))
  expect_true(nrow(airr_df) > 0)
})

test_that("'airr' format export includes retained full-length sequences", {
  seq_contigs <- list(data.frame(
    barcode = c("bc1", "bc1"),
    chain   = c("TRA", "TRB"),
    v_gene  = c("TRAV1", "TRBV1"),
    d_gene  = c(NA, "TRBD1"),
    j_gene  = c("TRAJ1", "TRBJ1"),
    c_gene  = c("TRAC", "TRBC1"),
    cdr3    = c("CA1", "CB1"),
    cdr3_nt = c("AAA", "TTT"),
    reads   = c(10, 12),
    sequence    = c("FULLA", "FULLB"),
    sequence_aa = c("FA", "FB"),
    stringsAsFactors = FALSE
  ))
  comb <- combineTCR(seq_contigs, retain.sequences = TRUE)
  airr_df <- exportClones(comb, format = "airr", write.file = FALSE)

  expect_true(all(c("sequence", "sequence_aa") %in% colnames(airr_df)))
  tra_row <- airr_df[airr_df$locus == "TRA1", ]
  expect_identical(tra_row$sequence, "FULLA")
  expect_identical(tra_row$sequence_aa, "FA")
})

test_that("'immunarch' format export works correctly", {
  skip_if_not_installed("dplyr")

  # Test returning a list object
  immunarch_list <- exportClones(combined, format = "immunarch", write.file = FALSE)
  expect_type(immunarch_list, "list")
  expect_named(immunarch_list, c("data", "meta"))
  expect_s3_class(immunarch_list$meta, "data.frame")
  expect_type(immunarch_list$data, "list")
  expect_true(length(immunarch_list$data) > 0)
  expect_named(immunarch_list$data[[1]], c("Clones", "Proportion", "CDR3.nt",
                                           "CDR3.aa", "V.name", "D.name",
                                           "J.name", "C.name", "Barcode"))
  
  # Test writing to file (with the bug fix)
  temp_dir <- tempdir()
  file_path <- file.path(temp_dir, "immunarch_test.csv")
  exportClones(combined, format = "immunarch", dir = temp_dir, file.name = "immunarch_test.csv")
  expect_true(file.exists(file_path))
  
  # Verify file content
  read_df <- read.csv(file_path)
  expect_s3_class(read_df, "data.frame")
  expect_true(nrow(read_df) > 0)
  expect_true("Sample" %in% names(read_df))
  unlink(file_path) # Clean up
})

# --- dowser bridge ------------------------------------------------------------

make_bcr_aligned_contigs <- function() {
  list(data.frame(
    barcode = c("bc1", "bc1", "bc2", "bc2"),
    chain   = c("IGH", "IGK", "IGH", "IGL"),
    v_gene  = c("IGHV1", "IGKV1", "IGHV1", "IGLV1"),
    d_gene  = c("IGHD1", NA, "IGHD1", NA),
    j_gene  = c("IGHJ1", "IGKJ1", "IGHJ1", "IGLJ1"),
    c_gene  = c("IGHG1", "IGKC", "IGHG1", "IGLC1"),
    cdr3    = c("CARH1", "CQKL1", "CARH1", "CQLL2"),
    cdr3_nt = c("TGTGCA", "TGTCAA", "TGTGCA", "TGTCAG"),
    reads   = c(20, 18, 22, 16),
    sequence_alignment = c("ACGT...H1", "ACGT...K1", "ACGT...H2", "ACGT...L2"),
    germline_alignment = c("ACGT...G1", "ACGT...G2", "ACGT...G3", "ACGT...G4"),
    stringsAsFactors = FALSE
  ))
}

test_that("exportDowser errors when alignment/germline columns are absent", {
  expect_error(
    exportDowser(combined, write.file = FALSE),
    "sequence_alignment"
  )
})

test_that("exportDowser returns an AIRR data frame with clone_id when sequences retained", {
  comb <- combineBCR(make_bcr_aligned_contigs(),
                     call.related.clones = FALSE,
                     retain.sequences = c("sequence_alignment", "germline_alignment"))
  dres <- exportDowser(comb, write.file = FALSE)
  expect_s3_class(dres, "data.frame")
  expect_true(all(c("sequence_id", "clone_id", "locus",
                    "sequence_alignment", "germline_alignment") %in% colnames(dres)))
  expect_true(nrow(dres) > 0)
  expect_false(any(duplicated(dres$sequence_id)))
})

test_that("exportClones dispatches to the dowser format", {
  comb <- combineBCR(make_bcr_aligned_contigs(),
                     call.related.clones = FALSE,
                     retain.sequences = c("sequence_alignment", "germline_alignment"))
  dres <- exportClones(comb, format = "dowser", write.file = FALSE)
  expect_s3_class(dres, "data.frame")
  expect_true("clone_id" %in% colnames(dres))
})
