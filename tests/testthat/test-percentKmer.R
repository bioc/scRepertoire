# test script for percentKmer.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list, 
                        samples = c("P17B", "P17L", "P18B", "P18L", 
                                     "P19B","P19L", "P20B", "P20L"))
                                    
test_that("percentKmer: output type is correct for aa", {
  
  # Test for ggplot output by default
  p <- percentKmer(combined, 
                   clone.call = "aa", 
                   motif.length = 3, 
                   min.depth = 1)
  expect_s3_class(p, "ggplot")
  
  # Test for matrix output when export.table = TRUE
  mat <- percentKmer(combined, 
                     clone.call = "aa", 
                     motif.length = 3, 
                     min.depth = 1, 
                     export.table = TRUE)
  expect_true(is.matrix(mat))
})

test_that("percentKmer: output type is correct for nt", {
  
  # Test for ggplot output by default
  p <- percentKmer(combined, 
                   clone.call = "nt", 
                   motif.length = 3, 
                   min.depth = 1)
  expect_s3_class(p, "ggplot")
  
  # Test for matrix output when export.table = TRUE
  mat <- percentKmer(combined, 
                     clone.call = "nt", 
                     motif.length = 3, 
                     min.depth = 1, 
                     export.table = TRUE)
  expect_true(is.matrix(mat))
})


test_that("percentKmer: matrix calculations are accurate", {
  
  mat <- percentKmer(combined, 
                     clone.call = "aa", 
                     motif.length = 3, 
                     min.depth = 1, 
                     top.motifs = NULL, 
                     export.table = TRUE)
  
  # Check dimensions
  expect_equal(nrow(mat), 8)
  expect_equal(ncol(mat), 5500)
  expect_equal(rownames(mat), c("P17B", "P17L", "P18B", "P18L", "P19B", 
                                "P19L", "P20B", "P20L"))
  
  # Check if all row sums equal 1 (for proportions)
  expect_true(all(abs(rowSums(mat) - 1) < 1e-9))
})


test_that("percentKmer: parameter handling works as expected", {
  
  mat_top1 <- percentKmer(combined, 
                          clone.call = "aa", 
                          motif.length = 3, 
                          min.depth = 1, 
                          top.motifs = 1,
                          export.table = TRUE)
  expect_true(ncol(mat_top1) == 1)
  expect_equal(colnames(mat_top1), "ASS")
  
  # Test `motif.length` parameter
  mat_len2 <- percentKmer(combined, 
                          clone.call = "aa", 
                          motif.length = 2, 
                          top.motifs = NULL, 
                          export.table = TRUE)
  expect_true(all(nchar(colnames(mat_len2)) == 2))
})


test_that("percentKmer: input validation and error handling", {
  
  # Test for error on invalid `clone.call`
  expect_error(
    percentKmer(combined, clone.call = "invalid_option"),
    "Please select either nucleotide (nt) or amino acid (aa) sequences for clone.call",
    fixed = TRUE
  )
})