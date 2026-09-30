cache_key <- function(prefix, parts) {
  raw <- paste(c(prefix, unlist(parts, use.names = TRUE)), collapse = "__")
  safe <- gsub("[^A-Za-z0-9._-]+", "-", raw)
  if (nchar(safe) > 180L) {
    tmp <- tempfile()
    writeBin(charToRaw(raw), tmp)
    on.exit(unlink(tmp), add = TRUE)
    safe <- paste0(substr(safe, 1L, 120L), "-", unname(tools::md5sum(tmp)))
  }
  safe
}

cache_path <- function(cfg, key) {
  if(length(key)!=1L || is.na(key) || !grepl("^[A-Za-z0-9._-]+$",key) || key %in% c(".",".."))stop("Invalid cache key")
  file.path(cfg$cache_dir, paste0(key, ".rds"))
}

cache_get <- function(cfg, key) {
  path <- cache_path(cfg, key)
  lock<-filelock::lock(paste0(path,".lock"),timeout=10000)
  if(is.null(lock))return(NULL)
  on.exit(filelock::unlock(lock),add=TRUE)
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

cache_put <- function(cfg, key, value) {
  path <- cache_path(cfg, key)
  lock<-filelock::lock(paste0(path,".lock"),timeout=10000)
  if(is.null(lock))stop("Cache is busy; retry shortly")
  on.exit(filelock::unlock(lock),add=TRUE)
  tmp <- tempfile("cache-write-",tmpdir=dirname(path),fileext=".tmp")
  on.exit(unlink(tmp),add=TRUE)
  saveRDS(value, tmp, version = 3)
  if(.Platform$OS.type=="windows" && file.exists(path))unlink(path)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Could not write cache file: ", path, call. = FALSE)
  }
  invisible(path)
}

cache_inventory <- function(cfg) {
  files <- list.files(cfg$cache_dir, pattern = "\\.rds$", full.names = TRUE)
  if (!length(files)) return(data.frame())
  rows <- lapply(files, function(path) {
    obj <- tryCatch(readRDS(path), error = function(e) NULL)
    data.frame(
      file = basename(path),
      created_at = obj$provenance$retrieved_at %||% format(file.info(path)$mtime, tz = "UTC"),
      dataset = obj$provenance$dataset_table %||% NA_character_,
      contrast = obj$provenance$contrast %||% NA_character_,
      rows = if (is.data.frame(obj$data)) nrow(obj$data) else NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
