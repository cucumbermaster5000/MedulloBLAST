# Connect Cloud sets R_CONFIG_ACTIVE. Explicit Cloud Variables take precedence.
apply_cloud_defaults <- function() {
  if(Sys.getenv('R_CONFIG_ACTIVE') != 'connect_cloud')return(invisible(NULL))
  values <- c(MB_MAX_WORKERS='1', MB_INTEGRATED_REPORT='false', MB_REPORT_PDF='false',
    MB_AUTO_PROVISION_DATA='true', MB_ALLOW_DATA_DOWNLOADS='false',
    MB_WORKER_TIMEOUT_SECONDS='600', MB_PROVISION_TIMEOUT_SECONDS='1800',
    OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
    R2_DEPMAP_CACHE_DIR='/tmp/medulloblast-v1/cache', MB_DATA_ROOT='/tmp/medulloblast-v1/references',
    NCBI_EMAIL='davengis@gmail.com')
  missing <- vapply(names(values),function(k)is.na(Sys.getenv(k,unset=NA_character_)),logical(1))
  if(any(missing))do.call(Sys.setenv,as.list(values[missing]))
  invisible(NULL)
}

# Opt-in diagnostics for Cloud capacity acceptance. Logs contain only resource
# counters and filesystem types, never environment variables or request data.
start_cloud_resource_diagnostics <- function() {
  if(Sys.getenv('MB_RESOURCE_DIAGNOSTICS')!='true'||Sys.info()[['sysname']]!='Linux')return(invisible(NULL))
  read_value <- function(path) if(file.exists(path))readLines(path,warn=FALSE)else character()
  started <- Sys.time()
  sample <- function() {
    tryCatch({
      stats <- read_value('/sys/fs/cgroup/memory.stat')
      stats <- stats[grepl('^(anon|file|shmem|inactive_file|active_file|kernel) ',stats)]
      children <- ps::ps_children(ps::ps_handle(),recursive=TRUE)
      processes <- lapply(c(list(ps::ps_handle()),children),function(p)tryCatch(list(pid=ps::ps_pid(p),name=ps::ps_name(p),rss=unname(ps::ps_memory_info(p)[['rss']])),error=function(e)NULL))
      data <- list(time=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC'),current=read_value('/sys/fs/cgroup/memory.current'),peak=read_value('/sys/fs/cgroup/memory.peak'),limit=read_value('/sys/fs/cgroup/memory.max'),stat=stats,processes=processes)
      cat('[Cloud resources] ',jsonlite::toJSON(data,auto_unbox=TRUE),'\n',sep='')
    },error=function(e)message('[Cloud resources] counters unavailable'))
    if(as.numeric(difftime(Sys.time(),started,units='mins'))<25)later::later(sample,10)
  }
  cat('[Cloud filesystems] ',paste(system2('df',c('-T','/tmp','/cloud/project'),stdout=TRUE),collapse=' | '),'\n',sep='')
  sample()
  invisible(NULL)
}
