# test script for clonalCluster.R 

library(testthat)
library(igraph)
library(SingleCellExperiment) 
library(Matrix)

# --- Global Setup ---
# Using the user-provided setup for general structure tests
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))

set.seed(42)
combined <- lapply(combined, function(x) {
  x <- x[sample(nrow(x), nrow(x) * .25),]
  x
})

# --- Existing Tests (Preserved) ---

test_that("Basic functionality and default output structure", {
  clustered_list <- clonalCluster(combined[1:2])
  expect_length(clustered_list, length(combined[1:2]))
  expect_true(all(sapply(clustered_list, is.data.frame)))
  expect_true("TRB.Cluster" %in% names(clustered_list[[1]]))
  
  # Check for data integrity
  expect_equal(nrow(clustered_list[[1]]), nrow(combined[[1]]))
  expect_false(all(is.na(clustered_list[[1]]$TRB.Cluster)))
})

test_that("exportGraph = TRUE returns a valid igraph object", {
  graph_obj <- clonalCluster(combined[1:2], exportGraph = TRUE)
  expect_s3_class(graph_obj, "igraph")
  expect_gt(vcount(graph_obj), 0)
  expect_gt(ecount(graph_obj), 0)
  expect_true("cluster" %in% igraph::vertex_attr_names(graph_obj))
  expect_true("weight" %in% igraph::edge_attr_names(graph_obj))
})

test_that("exportAdjMatrix = TRUE returns a valid sparse matrix", {
  adj_matrix <- clonalCluster(combined[3:4], exportAdjMatrix = TRUE)
  all_barcodes <- unique(do.call(rbind, combined[3:4])[["barcode"]])
  num_barcodes <- length(all_barcodes)
  expect_s4_class(adj_matrix, "dgCMatrix")
  expect_equal(dim(adj_matrix), c(num_barcodes, num_barcodes))
  expect_equal(rownames(adj_matrix), all_barcodes)
})

test_that("chain parameter works correctly", {
  clustered_tra <- clonalCluster(combined[5:6], chain = "TRA")
  expect_true("TRA.Cluster" %in% names(clustered_tra[[1]]))
  expect_false("TRB.Cluster" %in% names(clustered_tra[[1]]))
  clustered_both <- clonalCluster(combined[5:6], chain = "both")
  expect_true("Multi.Cluster" %in% names(clustered_both[[1]]))
})

test_that("group.by parameter functions without error", {
  clustered_grouped <- clonalCluster(combined[1:2], 
                                     group.by = "sample")
  expect_true("TRB.Cluster" %in% names(clustered_grouped[[1]]))
  expect_false(all(is.na(clustered_grouped[[1]]$TRB.Cluster)))
})

test_that("Different `cluster.method` options work", {
  louvain_graph <- clonalCluster(combined[5:6], 
                                 cluster.method = "louvain", 
                                 exportGraph = TRUE)
  expect_s3_class(louvain_graph, "igraph")
  expect_true("cluster" %in% igraph::vertex_attr_names(louvain_graph))
})

test_that("Input validation and error handling", {
  expect_error(
    clonalCluster(combined, exportGraph = TRUE, exportAdjMatrix = TRUE),
    "Please set only one of `exportGraph` or `exportAdjMatrix` to TRUE."
  )
  expect_error(
    clonalCluster(combined, cluster.method = "invalid_method"),
    "Unsupported clustering.method"
  )
})


test_that("Alignment (NW/SW) and Matrix selection", {
  # Create a small controlled dataset to ensure alignment logic runs
  toy_data <- combined[1]
  
  # Test Needleman-Wunsch with BLOSUM62
  res_nw <- clonalCluster(toy_data, 
                          dist_type = "nw", 
                          dist_mat = "BLOSUM62", 
                          threshold = 0.8) 
  expect_true("TRB.Cluster" %in% names(res_nw[[1]]))
  
  # Test Damerau (Transposition)
  res_dam <- clonalCluster(toy_data, dist_type = "damerau")
  expect_true("TRB.Cluster" %in% names(res_dam[[1]]))
})

test_that("Normalization parameters function correctly", {
  # normalize = "maxlen"
  res_max <- clonalCluster(combined[1], normalize = "maxlen", threshold = 0.1)
  expect_true("TRB.Cluster" %in% names(res_max[[1]]))
  
  # normalize = "length" (mean length)
  res_len <- clonalCluster(combined[1], normalize = "length", threshold = 0.1)
  expect_true("TRB.Cluster" %in% names(res_len[[1]]))
})