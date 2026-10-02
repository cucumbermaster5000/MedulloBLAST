# Shared in-process queue: waiting requests do not allocate R subprocesses.
explorer_create_worker_queue <- function(start_process=callr::r_bg, limit=1L) {
  pending<-list();running<-list();pumping<-FALSE
  pump<-function() {
    if(pumping)return(invisible(NULL))
    pumping<<-TRUE;on.exit(pumping<<-FALSE)
    completed<-FALSE
    running<<-Filter(function(j){
      if(!j$state %in% c('running','cancelling')){completed<<-TRUE;return(FALSE)}
      if(j$state=='cancelling'&&j$process$is_alive())tryCatch(j$process$kill(),error=function(e)FALSE)
      if(j$process$is_alive())return(TRUE)
      completed<<-TRUE
      j$state<-if(j$state=='cancelling')'cancelled' else 'finished';FALSE
    },running)
    # Reclaim completed requests before another worker overlaps with Shiny.
    if(completed)invisible(gc(full=TRUE))
    pending<<-Filter(function(j)j$state=='queued',pending)
    while(length(running)<limit && length(pending)) {
      j<-pending[[1]];pending<<-pending[-1]
      j$started<-Sys.time()
      j$process<-tryCatch(do.call(start_process,j$args),error=function(e){j$error<-e;NULL})
      j$args<-NULL
      if(is.null(j$process))j$state<-'failed' else {j$state<-'running';running[[length(running)+1L]]<<-j}
    }
    invisible(NULL)
  }
  submit<-function(args,owner=NULL) {
    j<-new.env(parent=emptyenv());j$args<-args;args<-NULL;j$owner<-owner
    j$process<-NULL;j$error<-NULL;j$state<-'queued';j$started<-NULL
    j$is_alive<-function(){pump();j$state %in% c('queued','running','cancelling')}
    j$started_at<-function()j$started
    j$phase<-function(){pump();j$state}
    j$kill<-function(){
      j$args<-NULL
      if(!is.null(j$process)&&j$process$is_alive()) {
        # Signal only our callr worker; kill_tree may also encounter protected
        # supervisor processes on Cloud. Never free a still-running slot.
        j$state<-'cancelling'
        tryCatch(j$process$kill(),error=function(e)FALSE)
        if(!j$process$is_alive())j$state<-'cancelled'
      } else j$state<-'cancelled'
      invisible(TRUE)
    }
    j$get_result<-function(){
      pump()
      if(!is.null(j$error))stop(j$error)
      if(j$state=='cancelled')stop('Analysis cancelled')
      if(j$state %in% c('queued','running','cancelling'))stop('Analysis is not complete')
      process<-j$process
      if(is.null(process))stop('Analysis result was already collected')
      on.exit(j$process<-NULL)
      process$get_result()
    }
    j$get_pid<-function()if(is.null(j$process))NA_integer_ else j$process$get_pid()
    pending[[length(pending)+1L]]<<-j
    pump();j
  }
  counts<-function(owner=NULL){
    pump()
    own<-function(j)is.null(owner)||identical(j$owner,owner)
    list(waiting=sum(vapply(pending,own,logical(1))),running=sum(vapply(running,own,logical(1))))
  }
  cancel_owner<-function(owner) {
    for(j in c(pending,running))if(identical(j$owner,owner))j$kill()
    invisible(NULL)
  }
  list(submit=submit,pump=pump,counts=counts,cancel_owner=cancel_owner)
}

explorer_worker_queue <- explorer_create_worker_queue(limit=1L)
explorer_r_bg <- function(...) {
  session<-shiny::getDefaultReactiveDomain()
  job<-explorer_worker_queue$submit(list(...),if(is.null(session))NULL else session$token)
  job
}

explorer_idle_seconds <- function() {
  x<-suppressWarnings(as.numeric(Sys.getenv('MB_SESSION_IDLE_SECONDS','900')))
  if(length(x)!=1L||!is.finite(x)||x<1)900 else x
}
explorer_capacity_notice <- function() shiny::div(id='hosting-notice',class='alert alert-info',
  shiny::strong('Shared free hosting: '),
  'Limited computing resources allow one heavy analysis task at a time. Other tasks wait in a queue. ',
  sprintf('Your session disconnects after %s minutes without interaction to free space for other users. ',format(explorer_idle_seconds()/60,trim=TRUE)),
  'Download any results you want to keep before leaving.',
  shiny::uiOutput('capacity_status'),
  shiny::div(id='mb-idle-warning',role='alert',style='display:none',
    shiny::span(id='mb-idle-warning-text'), ' ',
    shiny::tags$button(type='button',id='mb-keep-active',class='btn btn-sm btn-primary','Keep session active')))

explorer_session_policy <- function(input,output,session,now=Sys.time,limit=explorer_idle_seconds()) {
  # One owner callback, rather than retaining every completed job for the
  # lifetime of the session. Completed handles release their result process.
  session$onSessionEnded(function()explorer_worker_queue$cancel_owner(session$token))
  last<-now();warned<-FALSE;closing<-FALSE
  shiny::observeEvent(input$mb_user_activity,{
    if(!closing){last<<-now();warned<<-FALSE;session$sendCustomMessage('mb-session-policy',list(type='active'))}
  },ignoreInit=FALSE,ignoreNULL=TRUE,priority=100)
  shiny::observe({
    shiny::invalidateLater(1000,session)
    if(closing)return()
    left<-limit-as.numeric(difftime(now(),last,units='secs'))
    if(left<=0) {
      closing<<-TRUE
      session$allowReconnect(FALSE)
      session$sendCustomMessage('mb-session-policy',list(type='expired'))
      later::later(function(){if(!session$isClosed())session$close()},.25)
    } else if(left<=min(60,limit/3)&&!warned) {
      warned<<-TRUE
      session$sendCustomMessage('mb-session-policy',list(type='warning',seconds=ceiling(left)))
    }
  })
  output$capacity_status<-shiny::renderUI({
    shiny::invalidateLater(1000,session)
    n<-explorer_worker_queue$counts(session$token)
    if(n$waiting>0)shiny::tags$small(sprintf('Your analysis: %d task(s) waiting; %d running. Please keep this tab open.',n$waiting,n$running))
    else if(n$running>0)shiny::tags$small('Your analysis is running.')
  })
}
