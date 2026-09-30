# Reference recovery is separate from per-gene analysis and never runs preparation.
provisioning_catalog <- function(cfg) {
  path<-cfg$dataset_catalog %||% file.path(cfg$app_root,"config","datasets.json")
  x<-jsonlite::read_json(path,simplifyVector=FALSE)
  if(!identical(x$schema,1L)||!length(x$datasets))stop("Unsupported dataset catalog")
  safe<-function(z)is.character(z)&&length(z)==1L&&grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$",z)
  for(id in names(x$datasets)) {
    d<-x$datasets[[id]]
    if(!safe(id)||!safe(d$version)||!length(d$files))stop("Invalid dataset identity")
    names<-vapply(d$files,`[[`,character(1),"name")
    if(anyDuplicated(tolower(names)))stop("Duplicate dataset filenames")
    for(f in d$files) {
      if(!safe(f$name)||!is.numeric(f$bytes)||length(f$bytes)!=1L||!is.finite(f$bytes)||f$bytes<0||
         !is.character(f$sha256)||length(f$sha256)!=1L||!grepl("^[a-f0-9]{64}$",f$sha256))stop("Invalid dataset file")
      if(!is.null(f$bundled)) {
        if(!grepl("^resources/provenance/[A-Za-z0-9._-]+$",f$bundled))stop("Invalid bundled provenance path")
      } else if(is.null(f$url)||length(f$url)!=1L||
                !(grepl("^https://",f$url)||(isTRUE(cfg$provisioning_test_http)&&grepl("^http://127[.]0[.]0[.]1:[0-9]+/",f$url)))||
                grepl("/latest(/|$)",f$url))stop("Dataset URL must be HTTPS and version-pinned")
    }
  }
  x
}

provisioning_config <- function(cfg) {
  enabled<-tolower(Sys.getenv("MB_AUTO_PROVISION_DATA",if(Sys.getenv("R_CONFIG_ACTIVE")=="connect_cloud")"true"else"false"))=="true"
  cfg$auto_provision<-enabled
  if(!enabled)return(cfg)
  cfg$cache_dir<-Sys.getenv("R2_DEPMAP_CACHE_DIR",file.path(tempdir(),"medulloblast-cache"))
  dir.create(cfg$cache_dir,recursive=TRUE,showWarnings=FALSE)
  cfg$reference_root<-Sys.getenv("MB_DATA_ROOT",file.path(cfg$cache_dir,"references"))
  cfg$dataset_catalog<-file.path(cfg$app_root,"config","datasets.json")
  catalog<-provisioning_catalog(cfg)
  cfg$cerebellum_data_dir<-provisioning_path(cfg,"cerebellum",catalog)
  cfg$depmap_archive_dir<-provisioning_path(cfg,"depmap_archive",catalog)
  cfg
}

provisioning_path <- function(cfg,id,catalog=provisioning_catalog(cfg)) {
  if(!id %in% names(catalog$datasets))stop("Unknown reference dataset")
  file.path(cfg$reference_root,id,catalog$datasets[[id]]$version)
}

provisioning_fingerprint <- function(dataset) digest::digest(dataset,algo="sha256")
provisioning_state_config <- function(cfg) {
  p<-file.path(cfg$reference_root,"state");dir.create(p,recursive=TRUE,showWarnings=FALSE)
  list(cache_dir=p)
}
provisioning_set_state <- function(cfg,id,state,file=NULL,error=NULL) {
  cache_put(provisioning_state_config(cfg),paste0("reference-",id),
    list(state=state,file=file,error=error,updated_at=utc_now()))
}
provisioning_state <- function(cfg,id) cache_get(provisioning_state_config(cfg),paste0("reference-",id)) %||% list(state="queued")
provisioning_file_valid <- function(path,file,hash=TRUE) {
  file.exists(path)&&!dir.exists(path)&&isTRUE(file.info(path)$size==file$bytes)&&
    (!hash||identical(digest::digest(file=path,algo="sha256"),file$sha256))
}
provisioning_ready <- function(cfg,id,catalog=provisioning_catalog(cfg),hash=FALSE) {
  d<-catalog$datasets[[id]];p<-provisioning_path(cfg,id,catalog)
  marker_path<-file.path(p,"READY.rds")
  marker<-if(file.exists(marker_path))tryCatch(readRDS(marker_path),error=function(e)NULL)else NULL
  identical(marker,provisioning_fingerprint(d))&&all(vapply(d$files,function(f)
    provisioning_file_valid(file.path(p,f$name),f,hash),logical(1)))
}

provisioning_fetch <- function(url,path,timeout) {
  h<-curl::new_handle(timeout=timeout,connecttimeout=min(30,timeout),followlocation=TRUE,
    failonerror=TRUE,useragent="MedulloBLAST/0.4 (davengis@gmail.com)")
  response<-curl::curl_fetch_disk(url,path,h)
  if(response$status_code<200L||response$status_code>=300L)stop("Dataset service returned HTTP ",response$status_code)
}

provision_dataset <- function(cfg,id,catalog=provisioning_catalog(cfg),fetch=provisioning_fetch) {
  d<-catalog$datasets[[id]];final<-provisioning_path(cfg,id,catalog)
  parent<-dirname(final);dir.create(parent,recursive=TRUE,showWarnings=FALSE)
  lock<-filelock::lock(file.path(parent,"provision.lock"),timeout=0)
  if(is.null(lock))return(invisible("busy"))
  on.exit(filelock::unlock(lock),add=TRUE)
  budget<-cfg$provisioning_timeout %||% suppressWarnings(as.numeric(Sys.getenv("MB_PROVISION_TIMEOUT_SECONDS","1800")))
  if(length(budget)!=1L||!is.finite(budget)||budget<=0)budget<-1800
  deadline<-proc.time()[["elapsed"]]+budget
  remaining<-function(){n<-deadline-proc.time()[["elapsed"]];if(n<=0)stop("Reference provisioning deadline exceeded");n}
  tryCatch({
    provisioning_set_state(cfg,id,"verifying")
    if(provisioning_ready(cfg,id,catalog,hash=TRUE)) {
      provisioning_set_state(cfg,id,"ready");return(invisible("ready"))
    }
    if(dir.exists(final)) {
      # Both resolved paths stay under this catalog-validated dataset directory.
      base<-normalizePath(parent,winslash="/",mustWork=TRUE)
      from<-normalizePath(final,winslash="/",mustWork=TRUE)
      to<-gsub("\\\\","/",tempfile(paste0(d$version,"-invalid-"),tmpdir=base))
      if(!startsWith(from,paste0(base,"/"))||!startsWith(to,paste0(base,"/")))stop("Unsafe reference quarantine path")
      # A previously marked directory that fails revalidation is not ready,
      # even if Windows prevents renaming it because another reader is open.
      unlink(file.path(from,"READY.rds"))
      if(!file.rename(from,to))stop("Could not quarantine incomplete reference data")
    }
    stage<-paste0(final,".staging");dir.create(stage,recursive=TRUE,showWarnings=FALSE)
    for(f in d$files) {
      remaining();dest<-file.path(stage,f$name)
      if(provisioning_file_valid(dest,f))next
      part<-paste0(dest,".part");last<-NULL
      for(attempt in seq_len(3L)) {
        provisioning_set_state(cfg,id,"downloading",file=f$name)
        last<-tryCatch({
          timeout<-remaining()
          if(!is.null(f$bundled)) {
            if(!file.copy(file.path(cfg$app_root,f$bundled),part,overwrite=TRUE))stop("Bundled provenance missing")
          } else fetch(f$url,part,timeout)
          remaining()
          if(!provisioning_file_valid(part,f))stop("Reference file size or SHA-256 mismatch")
          if(file.exists(dest))unlink(dest)
          if(!file.rename(part,dest))stop("Could not finalize reference file")
          NULL
        },error=identity)
        if(is.null(last))break
        if(file.exists(part))unlink(part)
        if(attempt<3L)Sys.sleep(min(attempt,remaining()))
      }
      if(!is.null(last))stop(last)
    }
    remaining()
    if(!all(vapply(d$files,function(f)provisioning_file_valid(file.path(stage,f$name),f),logical(1))))stop("Reference validation failed")
    saveRDS(provisioning_fingerprint(d),file.path(stage,"READY.rds"))
    if(!file.rename(stage,final))stop("Could not publish verified reference directory")
    provisioning_set_state(cfg,id,"ready");invisible("ready")
  },error=function(e){
    provisioning_set_state(cfg,id,"failed",error=conditionMessage(e))
    message("[Reference provisioning] ",id,": ",conditionMessage(e));invisible("failed")
  })
}

provision_reference_data <- function(cfg) {
  catalog<-provisioning_catalog(cfg)
  setNames(lapply(names(catalog$datasets),function(id)provision_dataset(cfg,id,catalog)),names(catalog$datasets))
}
