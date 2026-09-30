# Administrator utility: inventory by default; prune derived records only while app is stopped.
.libPaths(c(".R-library",.libPaths()));source("R/load_core.R");core<-load_r2_core(".")
cfg<-core$read_server_config("config/defaults.json")
args<-commandArgs(trailingOnly=TRUE)
files<-list.files(cfg$cache_dir,pattern="[.]rds$",full.names=TRUE)
age<-as.numeric(difftime(Sys.time(),file.info(files)$mtime,units="days"))
# Preserve authoritative source records, rate-limit state and raw reference files.
eligible<-grepl("^(pubmed_(discovery|identity)|depmap_explorer|depmap_partial|biology_(current|result))",basename(files)) & is.finite(age) & age>90
cat("Cache root:",cfg$cache_dir,"\nDerived records older than 90 days:",sum(eligible),"\nBytes:",sum(file.info(files[eligible])$size),"\n")
if("--apply-offline" %in% args){
  root<-paste0(normalizePath(cfg$cache_dir,winslash="/"),"/")
  targets<-normalizePath(files[eligible],winslash="/",mustWork=TRUE)
  stopifnot(all(startsWith(targets,root)))
  # Non-recursive deletion only. Operator must stop all app workers first.
  unlink(targets)
  cat("Expired derived records removed; shared reference sources preserved.\n")
} else cat("Dry run. Stop all app workers before using --apply-offline.\n")
