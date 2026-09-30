# Bounded concurrency across workers on one shared persistent filesystem.
deployment_worker_slot <- function(cfg, timeout_seconds = NULL) {
  n<-as_integer_scalar(Sys.getenv("MB_MAX_WORKERS","1"),1,1,16)
  if(is.null(timeout_seconds)) {
    timeout_seconds<-suppressWarnings(as.numeric(Sys.getenv("MB_WORKER_TIMEOUT_SECONDS","600")))
    if(length(timeout_seconds)!=1L||!is.finite(timeout_seconds)||timeout_seconds<30)timeout_seconds<-600
  }
  stopifnot(length(timeout_seconds)==1L,is.finite(timeout_seconds),timeout_seconds>=0)
  folder<-file.path(cfg$cache_dir,"worker-slots");dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  # A query launches several independent modules. With one slot they must queue
  # for longer than a cold module's retrieval, rather than fail after 30 seconds.
  # The session's existing request deadline still bounds queue + execution time
  # and cancellation/session close still terminates the waiting callr process.
  deadline<-proc.time()[["elapsed"]]+timeout_seconds
  repeat {
    for(i in seq_len(n)) {
      lock<-filelock::lock(file.path(folder,paste0(i,".lock")),timeout=0)
      if(!is.null(lock))return(lock)
    }
    remaining<-deadline-proc.time()[["elapsed"]]
    if(remaining<=0)break
    Sys.sleep(min(.5,remaining))
  }
  stop("Analysis capacity wait timed out; retry shortly")
}

analyze_explorer_request <- function(cfg,gene) {
  retrieval<-get_r2_expression(cfg,gene)
  analysis<-analyze_r2_subgroups(retrieval)
  subtypes<-tryCatch(analyze_r2_subtypes(retrieval),error=function(e)NULL)
  clinical<-tryCatch(analyze_r2_clinical_v3(cfg,retrieval,analysis),error=function(e){message("Clinical analysis unavailable");NULL})
  list(gene=retrieval$provenance$queried_gene,retrieval=retrieval,analysis=analysis,
    subtypes=subtypes,subtype_failure=if(is.null(subtypes))"Molecular subtype analysis is unavailable." else NULL,clinical=clinical)
}
