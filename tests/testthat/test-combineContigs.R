# test script for combineContigs.R - testcases are NOT comprehensive!

test_that("combineTCR works with default parameters", {
  combined <- combineTCR(contig_list[1:2], samples = c("P17B", "P17L"))
  expect_type(combined, "list")
  expect_length(combined, 2)
  expect_s3_class(combined[[1]], "data.frame")
  # Check if barcodes are prefixed
  expect_true(startsWith(combined[[1]]$barcode[1], "P17B_"))
})

test_that("combineTCR `samples` and `ID` parameters work", {
  combined <- combineTCR(contig_list[1], samples = "S1", ID = "A")
  expect_equal(names(combined)[1], "S1_A")
  expect_true(startsWith(combined[[1]]$barcode[1], "S1_A_"))
})

test_that("combineTCR `filterNonproductive = FALSE` keeps non-productive chains", {
  contig_mock <- contig_list[[1]]
  contig_mock$productive[1:50] <- "False"
  combined_filtered <- combineTCR(list(contig_mock), samples="S1")
  combined_unfiltered <- combineTCR(list(contig_mock), samples="S1", filterNonproductive = FALSE)
  expect_lt(nrow(combined_filtered[[1]]), nrow(combined_unfiltered[[1]]))
})

test_that("combineTCR `removeNA` and `removeMulti` work", {
  contig_mock <- contig_list[[1]]
  combined_removeNA <- combineTCR(list(contig_mock), samples="S1", removeNA = TRUE)[[1]]
  expect_true(all(!grepl("NA_", combined_removeNA$CTaa)))
  expect_true(all(!grepl("_NA", combined_removeNA$CTnt)))
  
  combined_removeMulti <- combineTCR(list(contig_mock), samples="S1", removeMulti = TRUE)
  expect_true(all(!grepl(";", combined_removeMulti$CTaa)))
  expect_true(all(!grepl(";", combined_removeMulti$CTnt)))
})

# --- combineBCR testing -------------------------------------------------------

BCR_SOURCE <- read.csv("https://www.borch.dev/uploads/contigs/b_contigs.csv")


BCR_LIST <- list(P1 = BCR_SOURCE,
                 P2 = BCR_SOURCE)

BCR_LIST$P2$barcode <- paste0(BCR_LIST$P2$barcode, "_2")

test_that("Standard combineBCR functionality (Legacy & Basic)", {
  combined_bcr <- combineBCR(BCR_SOURCE, samples = "Patient1")
  
  expect_true(any(grepl("cluster", combined_bcr[[1]]$CTstrict)))
  expect_type(combined_bcr, "list")
  expect_length(combined_bcr, 1)
  expect_s3_class(combined_bcr[[1]], "data.frame")
  expect_true(startsWith(combined_bcr[[1]]$barcode[1], "Patient1_"))
  expect_true(all(c("cdr3_aa1", "cdr3_nt1", "CTgene", "CTnt") %in% colnames(combined_bcr[[1]])))
})

test_that("combineBCR with Alignment Metrics", {
  combined_nw <- combineBCR(BCR_LIST[1], 
                            samples = "Patient1",
                            dist_type = "nw", 
                            dist_mat = "BLOSUM62",
                            threshold = 0.85, # Normalized score
                            normalize = "length")
  
  expect_true("CTstrict" %in% colnames(combined_nw[[1]]))
  expect_true(any(grepl("^cluster", combined_nw[[1]]$CTstrict)))
  
  # Test Smith-Waterman (Local Alignment)
  combined_sw <- combineBCR(BCR_LIST[1], 
                            samples = "Patient1",
                            dist_type = "sw", 
                            dist_mat = "PAM30",
                            threshold = 2, # Raw score threshold
                            normalize = "none") # Raw score usually requires normalize='none'
  
  expect_true("CTstrict" %in% colnames(combined_sw[[1]]))
})

test_that("Clustering Logic: call.related.clones = FALSE", {
  combined_exact <- combineBCR(BCR_LIST[1], 
                               samples = "Patient1", 
                               call.related.clones = FALSE)
  
  sample_ct <- combined_exact[[1]]$CTstrict[1]
  expect_false(grepl("cluster", sample_ct))
  expect_true(grepl("_", sample_ct)) 
})

test_that("Output Structure and samples/ID handling", {
  # Test with Sample + ID
  combined_id <- combineBCR(BCR_LIST[1], 
                            samples = "P1", 
                            ID = "Timepoint1")
  
  first_barcode <- combined_id[[1]]$barcode[1]
  # Format should be: Sample_ID_Barcode
  expect_true(startsWith(first_barcode, "P1_Timepoint1_"))
  
  # Check column existence for specific BCR chains
  cols <- colnames(combined_id[[1]])
  expect_true(all(c("IGH", "IGLC") %in% cols))
})
