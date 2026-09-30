# One lightweight streaming provisioner per app; filesystem locks cover other apps.
explorer_provisioning_service <- function(cfg,core) {
  job<-NULL;attempted<-FALSE
  ensure<-function(retry=FALSE) {
    if(!isTRUE(cfg$auto_provision))return(invisible(NULL))
    if(!is.null(job)&&job$is_alive())return(invisible(NULL))
    if(attempted&&!retry)return(invisible(NULL))
    attempted<<-TRUE
    for(id in names(core$provisioning_catalog(cfg)$datasets))core$provisioning_set_state(cfg,id,"queued")
    job<<-callr::r_bg(function(cfg){source(file.path(cfg$app_root,"R/load_core.R"));core<-load_r2_core(cfg$app_root);core$provision_reference_data(cfg)},
      args=list(cfg=cfg),libpath=.libPaths(),supervise=TRUE)
    invisible(NULL)
  }
  status<-function() {
    if(!isTRUE(cfg$auto_provision))return(list())
    ensure()
    ids<-names(core$provisioning_catalog(cfg)$datasets)
    states<-setNames(lapply(ids,function(id)core$provisioning_state(cfg,id)),ids)
    if(!is.null(job)&&!job$is_alive()) {
      tryCatch(job$get_result(),error=function(e){
        for(id in ids)if(states[[id]]$state %in% c("queued","verifying","downloading"))core$provisioning_set_state(cfg,id,"failed",error="Provisioning process stopped")
      })
      states<-setNames(lapply(ids,function(id)core$provisioning_state(cfg,id)),ids)
    }
    states
  }
  list(ensure=ensure,status=status)
}

explorer_provisioning_server <- function(input,output,session,cfg,core,service) {
  # Poll without invalidating downstream plots when dataset state is unchanged.
  # In particular, disabled provisioning must not rebuild cerebellum controls
  # every second and repeatedly invalidate the integrated report snapshot.
  state<-if(is.null(service)||!isTRUE(cfg$auto_provision))shiny::reactive(list())else
    shiny::reactivePoll(1000,session,
      checkFunc=function()jsonlite::toJSON(service$status(),auto_unbox=TRUE),
      valueFunc=function()service$status())
  session$userData$reference_status<-state
  output$reference_status<-shiny::renderUI({
    states<-state();if(!length(states))return(NULL)
    labels<-c(cerebellum="Developing cerebellum",depmap_archive="Legacy DepMap archive")
    pending<-names(states)[vapply(states,function(s)s$state!="ready",logical(1))]
    if(!length(pending))return(NULL)
    failed<-any(vapply(states,function(s)s$state=="failed",logical(1)))
    shiny::div(class=paste("alert",if(failed)"alert-warning"else"alert-info"),role="status",
      shiny::strong(if(failed)"Reference data need a retry."else"Preparing reference datasets..."),
      lapply(pending,function(id)shiny::div(paste(labels[[id]],":",switch(states[[id]]$state,
        queued="waiting to download",verifying="checking saved files",downloading="downloading and verifying files",failed="download or integrity check failed",states[[id]]$state)))),
      shiny::p("Other analyses remain available. Dataset preparation is separate from gene analysis."),
      if(failed)shiny::actionButton("reference_retry","Retry dataset downloads"))
  })
  shiny::observeEvent(input$reference_retry,{if(!is.null(service))service$ensure(retry=TRUE)},ignoreInit=TRUE)
  invisible(state)
}

explorer_reference_state <- function(session,id) {
  f<-session$userData$reference_status
  if(!is.function(f))return(list(state="ready"))
  value<-f()[[id]]
  if(is.null(value))list(state="ready")else value
}
