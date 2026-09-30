# Connect Cloud sets R_CONFIG_ACTIVE. Explicit Cloud Variables take precedence.
apply_cloud_defaults <- function() {
  if(Sys.getenv('R_CONFIG_ACTIVE') != 'connect_cloud')return(invisible(NULL))
  values <- c(MB_MAX_WORKERS='1', MB_INTEGRATED_REPORT='false', MB_REPORT_PDF='false',
    MB_AUTO_PROVISION_DATA='true', MB_ALLOW_DATA_DOWNLOADS='false',
    MB_WORKER_TIMEOUT_SECONDS='600', MB_PROVISION_TIMEOUT_SECONDS='1800',
    OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
    R2_DEPMAP_CACHE_DIR='/tmp/medulloblast/cache', MB_DATA_ROOT='/tmp/medulloblast/references',
    NCBI_EMAIL='davengis@gmail.com')
  missing <- vapply(names(values),function(k)is.na(Sys.getenv(k,unset=NA_character_)),logical(1))
  if(any(missing))do.call(Sys.setenv,as.list(values[missing]))
  invisible(NULL)
}
