# Seurat v5 RPCA integration, called by seurat_rpca.py as
#   Rscript --vanilla seurat_rpca.R <params.json>
# Splits the RNA assay into one layer per batch, scales the data and computes a PCA on all
# exchanged genes (the integration features), then runs
# IntegrateLayers(method = RPCAIntegration). Returns the first `dims` dimensions of the
# integrated reduction: the dimensions used to find anchors, which Seurat's integration
# vignette also uses downstream.
params <- jsonlite::fromJSON(commandArgs(trailingOnly = TRUE)[1])
source(params$utils)
set.seed(params$seed)

obj <- read_seurat_input(params$input, scale_factor = params$scale_factor)
cat(ncol(obj), "cells,", nrow(obj), "genes,", length(unique(obj$batch)), "batches\n")

obj[["RNA"]] <- split(obj[["RNA"]], f = obj[["batch"]][, 1])
features <- rownames(obj)
VariableFeatures(obj) <- features
obj <- ScaleData(obj, features = features, verbose = FALSE)
obj <- RunPCA(obj, features = features, npcs = params$npcs, verbose = FALSE)

reference <- NULL
if (length(params$reference) > 0) {
  layers <- Layers(obj[["RNA"]], search = "data")
  reference <- match(paste0("data.", params$reference), layers)
  if (anyNA(reference)) {
    stop("reference batches not found: ", paste(params$reference[is.na(reference)], collapse = ", "))
  }
  cat("reference layers:", paste(layers[reference], collapse = ", "), "\n")
}

# a wrapper keeps the object out of the call that do.call() builds
integrate <- function(...) {
  IntegrateLayers(obj, method = RPCAIntegration, orig.reduction = "pca",
                  new.reduction = "integrated.rpca", dims = seq_len(params$dims),
                  reference = reference, verbose = TRUE, ...)
}
obj <- do.call(integrate, as.list(params$integrate_args))

embedding <- Embeddings(obj, "integrated.rpca")[colnames(obj), seq_len(params$dims), drop = FALSE]
write_embedding(params$output, embedding)
cat("embedding:", nrow(embedding), "cells x", ncol(embedding), "dimensions\n")
