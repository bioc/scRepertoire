# base R type check functions

isListOfNonEmptyDataFrames <- function(obj) {
    is.list(obj) &&
        all(sapply(obj, function(x) is.data.frame(x) && sum(dim(x)) > 0))
}
assertthat::on_failure(isListOfNonEmptyDataFrames) <- function(call, env) {
    paste0(deparse(call$obj), " is not a list of non-empty `data.frame`s")
}

# bio objects

is_seurat_object <- function(obj) inherits(obj, "Seurat")
assertthat::on_failure(is_seurat_object) <- function(call, env) {
    paste0(deparse(call$obj), " is not a Seurat object")
}

is_se_object <- function(obj) inherits(obj, "SummarizedExperiment")
assertthat::on_failure(is_se_object) <- function(call, env) {
    paste0(deparse(call$obj), " is not a SummarizedExperiment object")
}

is_seurat_or_se_object <- function(obj) {
    is_seurat_object(obj) || is_se_object(obj)
}
assertthat::on_failure(is_seurat_or_se_object) <- function(call, env) {
    paste0(deparse(call$obj), " is not a Seurat or SummarizedExperiment object")
}

# general objects

is_named_numeric <- function(obj) {
    is.numeric(obj) && !is.null(names(obj))
}
assertthat::on_failure(is_named_numeric) <- function(call, env) {
    paste0(deparse(call$obj), " is not a named numeric vector")
}

# functions

assertthat::on_failure(`%in%`) <- function(call, env) {
    paste0(deparse(call$x), " is not in ", deparse(call$table))
}
