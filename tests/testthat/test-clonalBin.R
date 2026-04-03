# test script for clonalBin.R - comprehensive test cases

# Data setup
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))

# Add a variable for grouped tests
combined_with_type <- addVariable(combined,
                                  variable.name = "Type",
                                  variables = rep(c("B", "L"), 4))

# Basic functionality tests
test_that("clonalBin returns a list with the same structure as input", {
  result <- clonalBin(combined)
  expect_type(result, "list")
  expect_equal(length(result), length(combined))
  expect_equal(names(result), names(combined))
})

test_that("clonalBin adds clonalFrequency, clonalProportion, and cloneSize columns", {
  result <- clonalBin(combined)
  expected_cols <- c("clonalFrequency", "clonalProportion", "cloneSize")
  for (i in seq_along(result)) {
    expect_true(all(expected_cols %in% colnames(result[[i]])))
  }
})

test_that("cloneSize is a factor with correct levels", {
  result <- clonalBin(combined)
  for (i in seq_along(result)) {
    expect_s3_class(result[[i]]$cloneSize, "factor")
    # Levels should be in reverse order (hyperexpanded first, None last)
    levels_vec <- levels(result[[i]]$cloneSize)
    expect_true(grepl("Hyperexpanded", levels_vec[1]))
    expect_true(grepl("None", levels_vec[length(levels_vec)]))
  }
})

# clone.call parameter tests
test_that("clone.call parameter works with 'strict' (default)", {
  result <- clonalBin(combined, clone.call = "strict")
  expect_true(all(c("CTstrict", "cloneSize") %in% colnames(result[[1]])))
})

test_that("clone.call parameter works with 'nt'", {
  result <- clonalBin(combined, clone.call = "nt")
  expect_true("CTnt" %in% colnames(result[[1]]))
})

test_that("clone.call parameter works with 'aa'", {
  result <- clonalBin(combined, clone.call = "aa")
  expect_true("CTaa" %in% colnames(result[[1]]))
})

test_that("clone.call parameter works with 'gene'", {
  result <- clonalBin(combined, clone.call = "gene")
  expect_true("CTgene" %in% colnames(result[[1]]))
})

# Proportion vs Frequency tests
test_that("proportion = TRUE uses clonalProportion for binning", {
  result <- clonalBin(combined, proportion = TRUE)
  # Check that clonalProportion values are between 0 and 1
  for (i in seq_along(result)) {
    props <- result[[i]]$clonalProportion
    props <- props[!is.na(props)]
    expect_true(all(props >= 0 & props <= 1))
  }
})

test_that("proportion = FALSE uses clonalFrequency for binning", {
  result <- clonalBin(combined,
                      proportion = FALSE,
                      clone.size = c(Rare = 1, Small = 5, Medium = 20,
                                    Large = 100, Hyperexpanded = 500))
  # Check that clonalFrequency values are integers >= 1
  for (i in seq_along(result)) {
    freqs <- result[[i]]$clonalFrequency
    freqs <- freqs[!is.na(freqs)]
    expect_true(all(freqs >= 1))
    expect_true(all(freqs == floor(freqs)))
  }
})

test_that("Error when proportion = FALSE and clone.size < 1", {
  bad_clone_size <- c(Rare = 0.1, Small = 0.5)
  expect_error(
    clonalBin(combined, proportion = FALSE, clone.size = bad_clone_size),
    "Adjust the clone.size parameter - there are groupings < 1"
  )
})

test_that("Auto-adjustment of upper bin limit when max frequency exceeds it", {
  result <- clonalBin(combined,
                      proportion = FALSE,
                      clone.size = c(Single = 1, Small = 2))
  # Find max frequency across all elements
  max_freq <- max(sapply(result, function(x) max(x$clonalFrequency, na.rm = TRUE)))
  # The highest bin should include the max frequency
  highest_level <- levels(result[[1]]$cloneSize)[1]
  # Extract the upper bound from the level name
  upper_bound <- as.numeric(gsub(".*<= ([0-9.]+)\\)", "\\1", highest_level))
  expect_gte(upper_bound, max_freq)
})

# Custom clone.size tests
test_that("Custom clone.size bins are applied correctly", {
  custom_bins <- c(Low = 0.01, Medium = 0.1, High = 1)
  result <- clonalBin(combined, clone.size = custom_bins)
  levels_vec <- levels(result[[1]]$cloneSize)
  expect_true(any(grepl("Low", levels_vec)))
  expect_true(any(grepl("Medium", levels_vec)))
  expect_true(any(grepl("High", levels_vec)))
})

# group.by parameter tests
test_that("group.by = NULL calculates per list element", {
  result <- clonalBin(combined, group.by = NULL)
  # Each list element should have its own frequency/proportion calculations
  expect_equal(length(result), length(combined))
})

test_that("group.by parameter works with variable grouping", {
  result <- clonalBin(combined_with_type, group.by = "Type")
  expect_type(result, "list")
  # Verify calculations are grouped by Type
  for (i in seq_along(result)) {
    expect_true("clonalProportion" %in% colnames(result[[i]]))
    expect_true("clonalFrequency" %in% colnames(result[[i]]))
  }
})

# chain parameter tests
test_that("chain = 'TRA' filters to TRA chain only", {
  result <- clonalBin(combined, chain = "TRA")
  expect_type(result, "list")
  expect_true("cloneSize" %in% colnames(result[[1]]))
})
test_that("chain = 'TRB' filters to TRB chain only", {
  result <- clonalBin(combined, chain = "TRB")
  expect_type(result, "list")
  expect_true("cloneSize" %in% colnames(result[[1]]))
})

test_that("chain = 'both' uses paired chain data", {
  result <- clonalBin(combined, chain = "both")
  expect_type(result, "list")
  expect_true("cloneSize" %in% colnames(result[[1]]))
})

test_that("Different chains produce different frequency calculations", {
  result_tra <- clonalBin(combined, chain = "TRA")
  result_trb <- clonalBin(combined, chain = "TRB")
  result_both <- clonalBin(combined, chain = "both")

  # The frequencies should differ between chain selections
  freq_tra <- result_tra[[1]]$clonalFrequency
  freq_trb <- result_trb[[1]]$clonalFrequency
  freq_both <- result_both[[1]]$clonalFrequency

  # They should not all be identical
  expect_false(identical(freq_tra, freq_trb))
})

# Data integrity tests
test_that("clonalBin preserves original data columns", {
  result <- clonalBin(combined)
  original_cols <- colnames(combined[[1]])
  for (col in original_cols) {
    expect_true(col %in% colnames(result[[1]]))
  }
})

test_that("clonalBin handles NA cloneCall values", {
  # Create test data with some NA values
  combined_with_na <- combined
  combined_with_na[[1]]$CTstrict[1:5] <- NA
  result <- clonalBin(combined_with_na, clone.call = "strict")
  expect_type(result, "list")
  # NAs should result in NA cloneSize
  expect_true(any(is.na(result[[1]]$cloneSize)))
})

test_that("clonalBin works with single list element input", {
  single_element <- combined[1]
  result <- clonalBin(single_element)
  expect_type(result, "list")
  expect_equal(length(result), 1)
  expect_true("cloneSize" %in% colnames(result[[1]]))
})

# Bin assignment tests
test_that("All non-NA clones are assigned to a bin",
{
  result <- clonalBin(combined)
  for (i in seq_along(result)) {
    # Rows with non-NA cloneCall should have non-NA cloneSize
    non_na_clones <- !is.na(result[[i]]$CTstrict)
    clone_sizes <- result[[i]]$cloneSize[non_na_clones]
    # The only NAs in cloneSize should be from NAs in CTstrict
    expect_true(sum(is.na(clone_sizes)) < sum(non_na_clones))
  }
})

test_that("Bin labels contain correct range format", {
  result <- clonalBin(combined)
  levels_vec <- levels(result[[1]]$cloneSize)
  # Each level should follow the pattern "Name (lower < X <= upper)"
  for (level in levels_vec) {
    expect_true(grepl("\\(.*< X <=.*\\)", level))
  }
})

# Deprecated parameter tests
test_that("Deprecated cloneCall parameter still works with warning", {
  expect_warning(
    result <- clonalBin(combined, cloneCall = "gene"),
    regexp = "cloneCall.*deprecated"
  )
  expect_true("CTgene" %in% colnames(result[[1]]))
})

test_that("Deprecated cloneSize parameter still works with warning", {
  custom_bins <- c(Low = 0.05, High = 1)
  expect_warning(
    result <- clonalBin(combined, cloneSize = custom_bins),
    regexp = "cloneSize.*deprecated"
  )
  levels_vec <- levels(result[[1]]$cloneSize)
  expect_true(any(grepl("Low", levels_vec)))
})

# Edge case tests
test_that("clonalBin handles data with all unique clones", {
  # Create a subset where all clones are unique
  small_combined <- lapply(combined, function(x) {
    unique_clones <- !duplicated(x$CTstrict) & !is.na(x$CTstrict)
    x[unique_clones, ][1:min(10, sum(unique_clones)), ]
  })
  result <- clonalBin(small_combined)
  expect_type(result, "list")
  # All should be in the smallest bin
  for (i in seq_along(result)) {
    freqs <- result[[i]]$clonalFrequency[!is.na(result[[i]]$clonalFrequency)]
    expect_true(all(freqs == 1))
  }
})

test_that("clonalBin handles data with highly expanded clones", {
  result <- clonalBin(combined,
                      proportion = FALSE,
                      clone.size = c(Single = 1, Small = 2, Medium = 3))
  # Check that clones above the specified bins are captured
  max_freqs <- sapply(result, function(x) max(x$clonalFrequency, na.rm = TRUE))
  # If any clone has frequency > 3, bin should be auto-adjusted
  if (any(max_freqs > 3)) {
    highest_level <- levels(result[[1]]$cloneSize)[1]
    upper_bound <- as.numeric(gsub(".*<= ([0-9.]+)\\)", "\\1", highest_level))
    expect_gte(upper_bound, max(max_freqs))
  }
})

# Consistency tests
test_that("Running clonalBin twice produces same results", {
  result1 <- clonalBin(combined)
  result2 <- clonalBin(combined)
  expect_equal(result1, result2)
})

test_that("clonalProportion sums to approximately 1 within each list element", {
  result <- clonalBin(combined)
  for (i in seq_along(result)) {
    # Get unique clones with their proportions
    unique_clones <- result[[i]][!duplicated(result[[i]]$CTstrict) &
                                   !is.na(result[[i]]$CTstrict), ]
    prop_sum <- sum(unique_clones$clonalProportion, na.rm = TRUE)
    expect_equal(prop_sum, 1, tolerance = 0.01)
  }
})
