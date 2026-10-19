# One lightweight streaming provisioner per app; filesystem locks cover other apps.
explorer_provisioning_service <- function(cfg,core,start_process=callr::r_bg,now=Sys.time) {
  job<-NULL;attempts<-0L;retry_at<-NULL
  ids<-if(isTRUE(cfg$auto_provision))names(core$provisioning_catalog(cfg)$datasets)else character()
  states<-function()setNames(lapply(ids,function(id)core$provisioning_state(cfg,id)),ids)
  ensure<-function(retry=FALSE) {
    if(!isTRUE(cfg$auto_provision)||!is.null(job))return(invisible(NULL))
    if(retry){attempts<<-0L;retry_at<<-NULL}
    if(attempts>=3L)return(invisible(NULL))
    if(!is.null(retry_at)&&now()<retry_at)return(invisible(NULL))
    current<-states()
    selected<-if(attempts==0L&&!retry)ids else names(current)[vapply(current,function(s)s$state!='ready',logical(1))]
    if(!length(selected))return(invisible(NULL))
    attempts<<-attempts+1L;retry_at<<-NULL
    for(id in selected)core$provisioning_set_state(cfg,id,"queued")
    job<<-tryCatch(start_process(function(cfg,ids){
      source(file.path(cfg$app_root,"R/load_core.R"));core<-load_r2_core(cfg$app_root)
      core$provision_reference_data(cfg,ids=ids)
    },args=list(cfg=cfg,ids=selected),libpath=.libPaths(),supervise=TRUE),error=function(e){
      for(id in selected)core$provisioning_set_state(cfg,id,"failed",error="Provisioning process could not start")
      NULL
    })
    invisible(NULL)
  }
  status<-function() {
    if(!isTRUE(cfg$auto_provision))return(list())
    if(!is.null(job)&&!job$is_alive()) {
      finished<-job;job<<-NULL
      tryCatch(finished$get_result(),error=function(e){
        current<-states()
        for(id in ids)if(current[[id]]$state %in% c("queued","verifying","downloading"))
          core$provisioning_set_state(cfg,id,"failed",error="Provisioning process stopped")
      })
    }
    current<-states()
    failed<-names(current)[vapply(current,function(s)s$state=="failed",logical(1))]
    if(is.null(job)&&length(failed)&&attempts>0L&&attempts<3L&&is.null(retry_at)) {
      retry_at<<-now()+c(10,30)[attempts]
      for(id in failed)core$provisioning_set_state(cfg,id,"retrying",error="Temporary download failure; automatic retry scheduled")
    }
    ensure();states()
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
        retrying="temporary download failure; retrying automatically",queued="waiting to download",verifying="checking saved files",downloading="downloading and verifying files",failed="download or integrity check failed",states[[id]]$state)))),
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
