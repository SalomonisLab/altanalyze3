#!/usr/bin/env Rscript
# Stream one Seurat assay's count matrix to standard output in the altanalyze3 sparse stream
# layout (sparse_stream.py), for cellHarmony_lite --sparse_stream -.
# R only reads the RDS and writes the dgCMatrix slots unchanged. No analysis happens here.
#
#   Rscript sparse_stream_from_seurat.R --rds obj.rds [options] | python -m ...cellHarmony_lite --sparse_stream - ...
#   Rscript sparse_stream_from_seurat.R --zip archive.zip --member obj.rds [options] | ...
#
# Options: --assay RNA (default: the active assay)  --slot counts  --obs-cols a,b,c
#          --library-col col  (copied to an obs column named Library, the ambient correction unit)
# Diagnostics go to standard error; standard output carries only the stream.

suppressPackageStartupMessages({library(SeuratObject); library(Matrix)})

args <- commandArgs(trailingOnly = TRUE)
opt <- list(rds = NULL, zip = NULL, member = NULL, assay = NULL, slot = "counts",
            `obs-cols` = "", `library-col` = NULL)
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (!key %in% names(opt) || i == length(args)) stop("unknown or valueless option: ", args[i])
  opt[[key]] <- args[i + 1]
  i <- i + 2
}
if (is.null(opt$rds) == is.null(opt$zip)) stop("pass exactly one of --rds or --zip")
if (!is.null(opt$zip) && is.null(opt$member)) stop("--zip needs --member")

t0 <- Sys.time()
if (!is.null(opt$rds)) {
  obj <- readRDS(opt$rds)
} else {
  cmd <- sprintf("unzip -p %s %s", shQuote(opt$zip), shQuote(opt$member))
  magic <- system(sprintf("%s | head -c 2 | od -An -tx1", cmd), intern = TRUE)
  magic <- gsub("[[:space:]]", "", paste(magic, collapse = ""))
  con <- pipe(cmd, "rb")
  if (magic == "1f8b") con <- gzcon(con) else if (!magic %in% c("580a", "410a", "420a"))
    stop("zip member is compressed as '", magic, "'; only gzip or uncompressed RDS can be streamed; extract it first")
  obj <- readRDS(con)
  close(con)
}
message(sprintf("[stream-R] loaded %s in %.1f s", class(obj)[1], as.numeric(Sys.time() - t0, units = "secs")))

assay <- if (is.null(opt$assay)) DefaultAssay(obj) else opt$assay
if (!assay %in% names(obj@assays)) stop("assay '", assay, "' not in object; has ", paste(names(obj@assays), collapse = ", "))
A <- obj@assays[[assay]]
m <- if (inherits(A, "Assay5")) LayerData(A, layer = opt$slot) else methods::slot(A, opt$slot)
if (!inherits(m, "dgCMatrix")) stop("assay ", assay, " slot ", opt$slot, " is ", class(m)[1], ", not dgCMatrix")
cells <- colnames(m); genes <- rownames(m)
if (anyDuplicated(cells)) stop(sum(duplicated(cells)), " duplicate cell names")

md <- obj@meta.data
if (!all(cells %in% rownames(md))) stop(sum(!cells %in% rownames(md)), " cells have no meta.data row")
obs_cols <- strsplit(opt$`obs-cols`, ",", fixed = TRUE)[[1]]
obs_cols <- obs_cols[nzchar(obs_cols)]
lib_col <- opt$`library-col`
missing <- setdiff(c(obs_cols, lib_col), colnames(md))
if (length(missing)) stop("meta.data lacks: ", paste(missing, collapse = ", "))
obs <- data.frame(obs_name = cells, stringsAsFactors = FALSE)
for (col in obs_cols) obs[[col]] <- as.character(md[cells, col])
if (!is.null(lib_col)) obs[["Library"]] <- as.character(md[cells, lib_col])
obs[is.na(obs)] <- ""
if (any(vapply(obs, function(v) any(grepl("[\t\n\r]", v)), logical(1)))) stop("an obs value holds a tab or newline")
if (any(grepl("\n", c(cells, genes)))) stop("a cell or gene name holds a newline")

x_int <- all(m@x == round(m@x))
data_code <- if (x_int && max(abs(m@x)) < 2^24) 3L else 4L    # float32 is exact for these counts
totals <- as.double(Matrix::colSums(m))
message(sprintf("[stream-R] %s %s: %d genes x %d cells, %d nonzero, sum %.0f, integer %s, data %s",
                assay, opt$slot, nrow(m), ncol(m), length(m@x), sum(totals), x_int,
                if (data_code == 3L) "float32" else "float64"))

out <- file("/dev/stdout", "wb")
wr_i64 <- function(v) for (x in v) {
  if (x < 0 || x >= 2^31) stop("value ", x, " does not fit the writer's int64 encoding")
  writeBin(c(as.integer(x), 0L), out, size = 4, endian = "little")
}
wr_text <- function(txt) {
  raw <- charToRaw(enc2utf8(txt))
  wr_i64(length(raw))
  writeBin(raw, out)
}
wr_chunked <- function(v, size) {
  n <- length(v); step <- 5e7
  if (n == 0) return(invisible())
  for (s in seq(1, n, by = step)) writeBin(v[s:min(n, s + step - 1)], out, size = size, endian = "little")
}

writeBin(charToRaw("AH3SPS01"), out)
wr_i64(c(ncol(m), nrow(m), length(m@x)))
writeBin(c(1L, 1L, data_code, 1L), out, size = 4, endian = "little")   # int32 indptr, int32 indices, flags: row totals
wr_text(paste(cells, collapse = "\n"))
wr_text(paste(genes, collapse = "\n"))
tsv <- paste(c(paste(colnames(obs), collapse = "\t"),
               do.call(paste, c(unname(as.list(obs)), sep = "\t"))), collapse = "\n")
wr_text(if (ncol(obs) > 1) paste0(tsv, "\n") else "")
wr_chunked(m@p, 4)
wr_chunked(m@i, 4)
wr_chunked(m@x, if (data_code == 3L) 4 else 8)
wr_chunked(totals, 8)
writeBin(charToRaw("AH3SPEND"), out)
flush(out); close(out)
message(sprintf("[stream-R] stream written in %.1f s", as.numeric(Sys.time() - t0, units = "secs")))
