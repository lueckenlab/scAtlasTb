# Helpers for integration methods that run in R on Seurat objects (seurat_rpca, stacas).
# Sourced by their R scripts; the exchange files are written and read by seurat_utils.py.
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(hdf5r)
  library(jsonlite)
})


# Seurat object (assay "RNA") with the raw counts and the log-normalized data of the exchanged
# genes, and the exchanged obs columns as metadata (empty strings become NA).
# The data are Seurat's LogNormalize, log1p(counts / library size * scale_factor), with each
# cell's library size over all genes (computed in Python before the genes were subset), i.e.
# the values NormalizeData() gives on the full object.
read_seurat_input <- function(file, scale_factor = 1e4) {
  h5 <- H5File$new(file, mode = "r")
  on.exit(h5$close_all())
  read <- function(name) h5[[name]][]

  cells <- read("obs_names")
  # Seurat does not allow underscores in feature names
  genes <- make.unique(gsub("_", "-", read("var_names")))
  shape <- read("counts/shape")
  stopifnot(shape[1] == length(cells), shape[2] == length(genes))
  # cells x genes CSR in the file = genes x cells CSC in R
  counts <- new("dgCMatrix", i = read("counts/indices"), p = read("counts/indptr"),
                x = read("counts/data"), Dim = c(length(genes), length(cells)),
                Dimnames = list(genes, cells))

  meta <- data.frame(row.names = cells)
  obs <- h5[["obs"]]
  for (key in names(obs)) {
    values <- obs[[key]][]
    values[values == ""] <- NA
    meta[[key]] <- values
  }

  data <- counts
  data@x <- log1p(data@x / rep.int(read("library_size"), diff(data@p)) * scale_factor)
  obj <- CreateSeuratObject(counts, meta.data = meta)
  SetAssayData(obj, layer = "data", new.data = data)
}


# Embedding (cells x dimensions, cell names as row names) for read_embedding() in
# seurat_utils.py, stored as a row-major vector with its shape
write_embedding <- function(file, embedding) {
  h5 <- H5File$new(file, mode = "w")
  on.exit(h5$close_all())
  h5[["obs_names"]] <- rownames(embedding)
  h5[["shape"]] <- dim(embedding)
  h5[["X_emb"]] <- as.vector(t(embedding))
}
