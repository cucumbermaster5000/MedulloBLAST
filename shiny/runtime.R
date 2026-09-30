# Infrastructure only; no statistical transformations.
explorer_session_config <- function(cfg,session) {
  folder<-tempfile("mb-session-",tmpdir=tempdir())
  dir.create(folder,recursive=TRUE)
  cfg$output_dir<-folder
  session$onSessionEnded(function(){
    # Only this generated session directory; never shared reference caches.
    later::later(function(){unlink(folder,recursive=TRUE)},delay=2)
  })
  cfg
}
explorer_worker_expired <- function(started,cfg) {
  limit<-suppressWarnings(as.numeric(Sys.getenv("MB_WORKER_TIMEOUT_SECONDS","600")))
  if(!is.finite(limit)||limit<30)limit<-600
  !is.null(started) && as.numeric(difftime(Sys.time(),started,units="secs"))>limit
}
explorer_log <- function(event,gene=NULL) {
  # Only controlled event names and validated gene symbols; no URLs, secrets or identifiers.
  message(format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC")," [Gene Explorer] ",event,
    if(!is.null(gene))paste0(" gene=",gene))
}

explorer_log_error <- function(stage,error) {
  detail<-conditionMessage(error)
  for(name in c("NCBI_API_KEY","NCBI_EMAIL")){
    secret<-Sys.getenv(name);if(nzchar(secret))detail<-gsub(secret,"[redacted]",detail,fixed=TRUE)
  }
  detail<-gsub("(?i)(api_key|token|authorization)[=:][^&[:space:]]+","credential=[redacted]",detail,perl=TRUE)
  explorer_log(paste(stage,detail))
}

start_explorer_worker <- function(cfg,gene,kind="main") {
  logfile<-tempfile(paste0(kind,"-"),tmpdir=cfg$output_dir,fileext=".log")
  callr::r_bg(function(cfg,gene,kind){
    source(file.path(cfg$app_root,"R","load_core.R"));core<-load_r2_core(cfg$app_root)
    lock<-core$deployment_worker_slot(cfg);on.exit(filelock::unlock(lock),add=TRUE)
    if(kind=="pfister")core$analyze_pfister_pan_cancer(core$get_pfister_expression(cfg,gene)) else core$analyze_explorer_request(cfg,gene)
  },args=list(cfg=cfg,gene=gene,kind=kind),libpath=.libPaths(),stdout=logfile,stderr=logfile,supervise=TRUE)
}
