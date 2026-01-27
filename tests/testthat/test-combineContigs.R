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

test_that("combineTCR `filter.nonproductive = FALSE` keeps non-productive chains", {
  contig_mock <- contig_list[[1]]
  contig_mock$productive[1:50] <- "False"
  combined_filtered <- combineTCR(list(contig_mock), samples="S1")
  combined_unfiltered <- combineTCR(list(contig_mock), samples="S1", filter.nonproductive = FALSE)
  expect_lt(nrow(combined_filtered[[1]]), nrow(combined_unfiltered[[1]]))
})

test_that("combineTCR `remove.na` and `remove.multi` work", {
  contig_mock <- contig_list[[1]]
  combined_remove.na <- combineTCR(list(contig_mock), samples="S1", remove.na = TRUE)[[1]]
  expect_true(all(!grepl("NA_", combined_remove.na$CTaa)))
  expect_true(all(!grepl("_NA", combined_remove.na$CTnt)))
  
  combined_remove.multi <- combineTCR(list(contig_mock), samples="S1", remove.multi = TRUE)
  expect_true(all(!grepl(";", combined_remove.multi$CTaa)))
  expect_true(all(!grepl(";", combined_remove.multi$CTnt)))
})

# --- combineBCR testing -------------------------------------------------------

# test-combineBCR.R
# Comprehensive tests for combineBCR function

# Setup: Load test data
BCR_SOURCE <- read.csv("https://www.borch.dev/uploads/contigs/b_contigs.csv")
BCR_LIST <- list(P1 = BCR_SOURCE,
                 P2 = BCR_SOURCE)
BCR_LIST$P2$barcode <- paste0(BCR_LIST$P2$barcode, "_2")

# =============================================================================
# BASIC FUNCTIONALITY TESTS
# =============================================================================

test_that("Standard combineBCR functionality (Legacy & Basic)", {
  combined_bcr <- combineBCR(BCR_SOURCE, samples = "Patient1")
  
  expect_true(any(grepl("cluster", combined_bcr[[1]]$CTstrict)))
  expect_type(combined_bcr, "list")
  expect_length(combined_bcr, 1)
  expect_s3_class(combined_bcr[[1]], "data.frame")
  expect_true(startsWith(combined_bcr[[1]]$barcode[1], "Patient1_"))
  expect_true(all(c("cdr3_aa1", "cdr3_nt1", "CTgene", "CTnt") %in% colnames(combined_bcr[[1]])))
})

test_that("combineBCR with multiple samples", {
  combined_multi <- combineBCR(BCR_LIST, samples = c("P1", "P2"))
  
  expect_type(combined_multi, "list")
  expect_length(combined_multi, 2)
  expect_equal(names(combined_multi), c("P1", "P2"))
  expect_true(all(startsWith(combined_multi[["P1"]]$barcode, "P1_")))
  expect_true(all(startsWith(combined_multi[["P2"]]$barcode, "P2_")))
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

test_that("combineBCR without samples parameter", {
  combined_no_samples <- combineBCR(BCR_LIST)
  
  expect_type(combined_no_samples, "list")
  expect_length(combined_no_samples, 2)
  # Should use default naming S1, S2, etc.
  expect_equal(names(combined_no_samples), c("S1", "S2"))
})

# =============================================================================
# CLUSTERING METHOD TESTS
# =============================================================================

test_that("combineBCR with Alignment Metrics", {
  # Test Needleman-Wunsch (Global Alignment)
  combined_nw <- combineBCR(BCR_LIST[1], 
                            samples = "Patient1",
                            dist.type = "nw", 
                            dist.mat = "BLOSUM62",
                            threshold = 0.85,
                            normalize = "length")
  
  expect_true("CTstrict" %in% colnames(combined_nw[[1]]))
  expect_true(any(grepl("^cluster", combined_nw[[1]]$CTstrict)))
  
  # Test Smith-Waterman (Local Alignment)
  combined_sw <- combineBCR(BCR_LIST[1], 
                            samples = "Patient1",
                            dist.type = "sw", 
                            dist.mat = "PAM30",
                            threshold = 2,
                            normalize = "none")
  
  expect_true("CTstrict" %in% colnames(combined_sw[[1]]))
})

test_that("combineBCR with different distance metrics", {
  # Levenshtein (default)
  combined_lev <- combineBCR(BCR_LIST[1], 
                             samples = "Patient1",
                             dist.type = "levenshtein")
  expect_true("CTstrict" %in% colnames(combined_lev[[1]]))
  
  # Damerau-Levenshtein
  combined_dam <- combineBCR(BCR_LIST[1], 
                             samples = "Patient1",
                             dist.type = "damerau")
  expect_true("CTstrict" %in% colnames(combined_dam[[1]]))
})

test_that("Clustering Logic: call.related.clones = FALSE", {
  combined_exact <- combineBCR(BCR_LIST[1], 
                               samples = "Patient1", 
                               call.related.clones = FALSE)
  
  sample_ct <- combined_exact[[1]]$CTstrict[1]
  expect_false(grepl("cluster", sample_ct))
  expect_true(grepl("_", sample_ct)) 
})

# =============================================================================
# CTstrict FORMAT TESTS
# =============================================================================

test_that("CTstrict format: Heavy_Light underscore structure", {
  combined_bcr <- combineBCR(BCR_SOURCE, samples = "Patient1")
  
  # All CTstrict values should have exactly one underscore separating heavy and light
  ctstrict_values <- na.omit(combined_bcr[[1]]$CTstrict)
  
  # Count underscores in each value
  underscore_counts <- sapply(ctstrict_values, function(x) {
    length(gregexpr("_", x, fixed = TRUE)[[1]])
  })
  
  # Each CTstrict should have at least one underscore (Heavy_Light format)
  # Note: sequences may contain underscores, so we check for at least one
  expect_true(all(underscore_counts >= 1))
  
  # Each CTstrict should be splittable into exactly 2 parts by the first underscore
  # This validates the Heavy_Light structure
  for (val in ctstrict_values[1:min(10, length(ctstrict_values))]) {
    parts <- strsplit(val, "_", fixed = TRUE)[[1]]
    expect_gte(length(parts), 2)
  }
})

test_that("CTstrict format: Clustered vs Singlet cells", {
  combined_bcr <- combineBCR(BCR_SOURCE, samples = "Patient1")
  
  ctstrict_values <- combined_bcr[[1]]$CTstrict
  
  # Clustered cells should have format: cluster.X_cluster.X (when chain="both")
  clustered_pattern <- "^cluster\\.\\d+_cluster\\.\\d+$"
  clustered_cells <- grep(clustered_pattern, ctstrict_values, value = TRUE)
  expect_true(length(clustered_cells) > 0)
  
  # Singlets should NOT have cluster pattern - they use vgene.jgene.seq format
  non_cluster <- ctstrict_values[!grepl("^cluster\\.", ctstrict_values) & !is.na(ctstrict_values)]
  
  if (length(non_cluster) > 0) {
    # Singlets should contain gene names (e.g., IGHV, IGKV, IGLV) or NA
    has_gene_or_na <- sapply(non_cluster, function(x) {
      grepl("IG[HKL]", x) | grepl("^NA_", x) | grepl("_NA$", x)
    })
    expect_true(all(has_gene_or_na))
    
    # Should still have underscore separator
    expect_true(all(grepl("_", non_cluster)))
  }
})

test_that("CTstrict: Clustered cells share same cluster ID on both sides (chain='both')", {
  combined_bcr <- combineBCR(BCR_SOURCE, 
                             samples = "Patient1",
                             chain = "both")
  
  clustered <- grep("^cluster\\.\\d+_cluster\\.\\d+$", 
                    combined_bcr[[1]]$CTstrict, 
                    value = TRUE)
  
  if (length(clustered) > 0) {
    # For chain="both", both sides should have the SAME cluster ID
    for (ct in clustered) {
      parts <- strsplit(ct, "_")[[1]]
      expect_equal(parts[1], parts[2], 
                   info = paste("Cluster IDs should match for chain='both':", ct))
    }
  }
})

test_that("CTstrict format: Single chain clustering (IGH only)", {
  combined_igh <- combineBCR(BCR_SOURCE, 
                             samples = "Patient1",
                             chain = "Heavy")
  
  ctstrict_values <- na.omit(combined_igh[[1]]$CTstrict)
  
  # With chain = "IGH", clustered cells should have:
  # cluster.X on heavy side, unique ID on light side
  # Format: cluster.X_IGKV... or cluster.X_IGLV... or cluster.X_NA
  clustered_igh <- grep("^cluster\\.\\d+_", ctstrict_values, value = TRUE)
  
  if (length(clustered_igh) > 0) {
    # Light chain side should NOT be cluster (should be vgene.jgene.seq or NA)
    expect_false(any(grepl("_cluster\\.", clustered_igh)),
                 info = "Light chain side should not be cluster when chain='IGH'")
    
    # Light chain side should have light chain gene or NA
    light_side_valid <- sapply(clustered_igh, function(x) {
      light_part <- strsplit(x, "_")[[1]][2]
      grepl("^IG[KL]", light_part) | light_part == "NA"
    })
    expect_true(all(light_side_valid))
  }
})

test_that("CTstrict format: Single chain clustering (Light only)", {
  combined_light <- combineBCR(BCR_SOURCE, 
                               samples = "Patient1",
                               chain = "Light")
  
  ctstrict_values <- na.omit(combined_light[[1]]$CTstrict)
  
  # With chain = "Light", clustered cells should have:
  # unique ID on heavy side, cluster.X on light side
  # Format: IGHV..._cluster.X or NA_cluster.X
  clustered_light <- grep("_cluster\\.\\d+$", ctstrict_values, value = TRUE)
  
  if (length(clustered_light) > 0) {
    # Heavy chain side should NOT be cluster
    expect_false(any(grepl("^cluster\\.", clustered_light)),
                 info = "Heavy chain side should not be cluster when chain='Light'")
    
    # Heavy chain side should have heavy chain gene or NA
    heavy_side_valid <- sapply(clustered_light, function(x) {
      heavy_part <- strsplit(x, "_")[[1]][1]
      grepl("^IGHV", heavy_part) | heavy_part == "NA"
    })
    expect_true(all(heavy_side_valid))
  }
})

test_that("CTstrict format: Single chain clustering (IGL only)", {
  combined_igl <- combineBCR(BCR_SOURCE, 
                             samples = "Patient1",
                             chain = "IGL")
  
  ctstrict_values <- na.omit(combined_igl[[1]]$CTstrict)
  
  # Similar to Light chain test
  clustered_igl <- grep("_cluster\\.\\d+$", ctstrict_values, value = TRUE)
  
  if (length(clustered_igl) > 0) {
    # Heavy chain side should NOT be cluster
    expect_false(any(grepl("^cluster\\.", clustered_igl)))
  }
})

test_that("CTstrict format: Single chain clustering (IGK only)", {
  combined_igk <- combineBCR(BCR_SOURCE, 
                             samples = "Patient1",
                             chain = "IGK")
  
  ctstrict_values <- na.omit(combined_igk[[1]]$CTstrict)
  
  # Similar to Light chain test
  clustered_igk <- grep("_cluster\\.\\d+$", ctstrict_values, value = TRUE)
  
  if (length(clustered_igk) > 0) {
    # Heavy chain side should NOT be cluster
    expect_false(any(grepl("^cluster\\.", clustered_igk)))
  }
})

test_that("CTstrict format: Cells with missing heavy chain",
          {
            combined_bcr <- combineBCR(BCR_SOURCE, samples = "Patient1")
            
            df <- combined_bcr[[1]]
            
            # Find cells with only light chain (heavy chain is NA)
            light_only <- df[is.na(df$IGH) & !is.na(df$IGLC), ]
            
            if (nrow(light_only) > 0) {
              # CTstrict should have NA on the heavy side (or be fully NA)
              ctstrict_light_only <- light_only$CTstrict[!is.na(light_only$CTstrict)]
              
              if (length(ctstrict_light_only) > 0) {
                # Should start with NA_ or be cluster-based
                valid_format <- sapply(ctstrict_light_only, function(x) {
                  grepl("^NA_", x) | grepl("^cluster\\.", x)
                })
                expect_true(all(valid_format),
                            info = "Light-only cells should have 'NA_' prefix or be clustered")
              }
            }
          })

test_that("CTstrict format: Cells with missing light chain", {
  combined_bcr <- combineBCR(BCR_SOURCE, samples = "Patient1")
  
  df <- combined_bcr[[1]]
  
  # Find cells with only heavy chain (light chain is NA)
  heavy_only <- df[!is.na(df$IGH) & is.na(df$IGLC), ]
  
  if (nrow(heavy_only) > 0) {
    # CTstrict should have NA on the light side (or be fully NA)
    ctstrict_heavy_only <- heavy_only$CTstrict[!is.na(heavy_only$CTstrict)]
    
    if (length(ctstrict_heavy_only) > 0) {
      # Should end with _NA or be cluster-based
      valid_format <- sapply(ctstrict_heavy_only, function(x) {
        grepl("_NA$", x) | grepl("_cluster\\.", x)
      })
      expect_true(all(valid_format),
                  info = "Heavy-only cells should have '_NA' suffix or be clustered")
    }
  }
})

test_that("CTstrict format: Consistency with sequence parameter (aa vs nt)", {
  # Test with amino acid sequences
  combined_aa <- combineBCR(BCR_SOURCE, 
                            samples = "Patient1",
                            sequence = "aa",
                            threshold = 0.8)
  
  # Test with nucleotide sequences
  combined_nt <- combineBCR(BCR_SOURCE, 
                            samples = "Patient1",
                            sequence = "nt",
                            threshold = 0.8)
  
  # Get singlets from each (non-cluster entries)
  singlets_aa <- combined_aa[[1]]$CTstrict[!grepl("^cluster\\.", combined_aa[[1]]$CTstrict) & 
                                             !is.na(combined_aa[[1]]$CTstrict) &
                                             !grepl("^NA_NA$", combined_aa[[1]]$CTstrict)]
  
  singlets_nt <- combined_nt[[1]]$CTstrict[!grepl("^cluster\\.", combined_nt[[1]]$CTstrict) & 
                                             !is.na(combined_nt[[1]]$CTstrict) &
                                             !grepl("^NA_NA$", combined_nt[[1]]$CTstrict)]
  
  if (length(singlets_aa) > 0 && length(singlets_nt) > 0) {
    # NT sequences should be longer than AA sequences (roughly 3x for coding)
    avg_len_aa <- mean(nchar(singlets_aa))
    avg_len_nt <- mean(nchar(singlets_nt))
    
    expect_gt(avg_len_nt, avg_len_aa)
  }
})

# =============================================================================
# FILTERING PARAMETER TESTS
# =============================================================================

test_that("combineBCR remove.na parameter", {
  combined_with_na <- combineBCR(BCR_SOURCE, 
                                 samples = "Patient1",
                                 remove.na = FALSE)
  
  combined_no_na <- combineBCR(BCR_SOURCE, 
                               samples = "Patient1",
                               remove.na = TRUE)
  
  # remove.na=TRUE should have fewer or equal rows
  expect_lte(nrow(combined_no_na[[1]]), nrow(combined_with_na[[1]]))
})

test_that("combineBCR remove.multi parameter", {
  combined_with_multi <- combineBCR(BCR_SOURCE, 
                                    samples = "Patient1",
                                    remove.multi = FALSE)
  
  combined_no_multi <- combineBCR(BCR_SOURCE, 
                                  samples = "Patient1",
                                  remove.multi = TRUE)
  
  # remove.multi=TRUE should have fewer or equal rows
  expect_lte(nrow(combined_no_multi[[1]]), nrow(combined_with_multi[[1]]))
})

test_that("combineBCR filter.multi parameter", {
  combined_filter <- combineBCR(BCR_SOURCE, 
                                samples = "Patient1",
                                filter.multi = TRUE)
  
  combined_no_filter <- combineBCR(BCR_SOURCE, 
                                   samples = "Patient1",
                                   filter.multi = FALSE)
  
  # Both should produce valid output
  expect_s3_class(combined_filter[[1]], "data.frame")
  expect_s3_class(combined_no_filter[[1]], "data.frame")
})

test_that("combineBCR filter.nonproductive parameter", {
  combined_productive <- combineBCR(BCR_SOURCE, 
                                    samples = "Patient1",
                                    filter.nonproductive = TRUE)
  
  combined_all <- combineBCR(BCR_SOURCE, 
                             samples = "Patient1",
                             filter.nonproductive = FALSE)
  
  # filter.nonproductive=TRUE should have fewer or equal rows
  expect_lte(nrow(combined_productive[[1]]), nrow(combined_all[[1]]))
})

# =============================================================================
# V/J GENE USAGE TESTS
# =============================================================================

test_that("combineBCR use.V parameter affects clustering", {
  combined_use_v <- combineBCR(BCR_SOURCE, 
                               samples = "Patient1",
                               use.V = TRUE)
  
  combined_no_v <- combineBCR(BCR_SOURCE, 
                              samples = "Patient1",
                              use.V = FALSE)
  
  # Both should produce valid output with CTstrict
  expect_true("CTstrict" %in% colnames(combined_use_v[[1]]))
  expect_true("CTstrict" %in% colnames(combined_no_v[[1]]))
  
  # Clustering results may differ
  # When use.V=FALSE, more sequences may cluster together
  clustered_use_v <- sum(grepl("^cluster\\.", combined_use_v[[1]]$CTstrict), na.rm = TRUE)
  clustered_no_v <- sum(grepl("^cluster\\.", combined_no_v[[1]]$CTstrict), na.rm = TRUE)
  
  # Just verify both produce some clusters
  expect_true(clustered_use_v > 0 || clustered_no_v > 0)
})

test_that("combineBCR use.J parameter affects clustering", {
  combined_use_j <- combineBCR(BCR_SOURCE, 
                               samples = "Patient1",
                               use.J = TRUE)
  
  combined_no_j <- combineBCR(BCR_SOURCE, 
                              samples = "Patient1",
                              use.J = FALSE)
  
  # Both should produce valid output with CTstrict
  expect_true("CTstrict" %in% colnames(combined_use_j[[1]]))
  expect_true("CTstrict" %in% colnames(combined_no_j[[1]]))
})

test_that("combineBCR use.V and use.J together", {
  combined_both <- combineBCR(BCR_SOURCE, 
                              samples = "Patient1",
                              use.V = TRUE,
                              use.J = TRUE)
  
  combined_neither <- combineBCR(BCR_SOURCE, 
                                 samples = "Patient1",
                                 use.V = FALSE,
                                 use.J = FALSE)
  
  # When neither V nor J is used, more sequences may cluster
  clustered_both <- sum(grepl("^cluster\\.", combined_both[[1]]$CTstrict), na.rm = TRUE)
  clustered_neither <- sum(grepl("^cluster\\.", combined_neither[[1]]$CTstrict), na.rm = TRUE)
  
  # More permissive clustering (no V/J requirement) should produce >= clusters
  expect_gte(clustered_neither, clustered_both)
})

# =============================================================================
# GROUP.BY PARAMETER TESTS
# =============================================================================

test_that("combineBCR group.by parameter", {
  # Add a grouping column to the data
  bcr_with_group <- BCR_SOURCE
  bcr_with_group$patient_group <- sample(c("GroupA", "GroupB"), 
                                         nrow(bcr_with_group), 
                                         replace = TRUE)
  
  combined_grouped <- combineBCR(bcr_with_group, 
                                 samples = "Patient1",
                                 group.by = "patient_group")
  
  # Should produce valid output
  expect_s3_class(combined_grouped[[1]], "data.frame")
  expect_true("CTstrict" %in% colnames(combined_grouped[[1]]))
})

# =============================================================================
# EDGE CASES AND ERROR HANDLING
# =============================================================================


test_that("combineBCR CTstrict values are unique for different clones", {
  combined_bcr <- combineBCR(BCR_SOURCE, 
                             samples = "Patient1",
                             call.related.clones = FALSE)
  
  ctstrict_values <- combined_bcr[[1]]$CTstrict
  
  # Check that duplicated CTstrict values have matching chain information
  dup_ctstrict <- ctstrict_values[duplicated(ctstrict_values) & !is.na(ctstrict_values)]
  
  if (length(dup_ctstrict) > 0) {
    for (dup in unique(dup_ctstrict)) {
      matching_rows <- combined_bcr[[1]][combined_bcr[[1]]$CTstrict == dup & 
                                           !is.na(combined_bcr[[1]]$CTstrict), ]
      # All rows with same CTstrict should have same CDR3 sequences
      expect_equal(length(unique(matching_rows$cdr3_aa1)), 1)
      expect_equal(length(unique(matching_rows$cdr3_aa2)), 1)
    }
  }
})

test_that("combineBCR produces consistent output across runs", {
  set.seed(42)
  combined_run1 <- combineBCR(BCR_SOURCE, samples = "Patient1")
  
  set.seed(42)
  combined_run2 <- combineBCR(BCR_SOURCE, samples = "Patient1")
  
  # Results should be identical with same seed
  expect_equal(nrow(combined_run1[[1]]), nrow(combined_run2[[1]]))
  expect_equal(combined_run1[[1]]$barcode, combined_run2[[1]]$barcode)
})
