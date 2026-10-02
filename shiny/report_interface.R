explorer_reports_enabled <- function() tolower(Sys.getenv('MB_INTEGRATED_REPORT','true')) == 'true'

explorer_report_control <- function(id,input=NULL) {
  if(!explorer_reports_enabled()) return(NULL)
  value<-if(is.null(input))TRUE else shiny::isolate(input[[paste0('report_include_',id)]])
  shiny::checkboxInput(paste0('report_include_',id),'Include in Gene Report',value=!identical(value,FALSE))
}

explorer_report_ids <- function() c('pfister_pan_cancer','r2_subgroups','r2_subtypes','depmap_dependency','depmap_expression_dependency','tf_target_enrichment','functional_annotation',
  paste0('hpa_',c('tissue','brain','single_cell','subcellular','cancer','blood','cell_line','structure','interaction')),
  outer(c('metastasis','survival'),c('all','wnt','shh','group3','group4'),paste,sep='_'),
  outer(c('cerebellum_annotation','cerebellum_expression'),c('human','mouse'),paste,sep='_'))

explorer_report_snapshot <- function(cfg, input, result, pfister, biology, hpa, depmap, publications,cerebellum=function()NULL) {
  x <- result(); shiny::req(x)
  or <- function(a,b) if(is.null(a)||!length(a)) b else a
  current <- function(r) tryCatch(r(), error=function(e) NULL)
  pf <- current(pfister); if(!is.null(pf) && !identical(pf$gene,x$gene)) pf <- NULL
  if(!is.null(pf)) pf <- pf$analysis
  b <- current(biology); if(!is.null(b) && !identical(b$gene,x$gene)) b <- NULL
  hp <- current(hpa); if(!is.null(hp) && !identical(hp$identity$canonical_symbol,x$gene)) hp <- NULL
  dm <- current(depmap); if(!is.null(dm) && !identical(dm$canonical_symbol,x$gene)) dm <- NULL
  pub <- current(publications); if(!is.null(pub) && !identical(pub$canonical_symbol,x$gene)) pub <- NULL
  cb <- current(cerebellum); if(!is.null(cb)&&!identical(toupper(cb$query),toupper(x$gene)))cb<-NULL
  version <- tryCatch(as.character(read.dcf(file.path(cfg$app_root,"DESCRIPTION"),fields="Version")[[1]]),error=function(e)"unknown")
  commit <- tryCatch(trimws(system2("git",c("-C",shQuote(cfg$app_root),"rev-parse","HEAD"),stdout=TRUE,stderr=FALSE)[1]),error=function(e)NA_character_)
  sources <- list(
    R2 = x$retrieval$provenance[c("dataset_label","dataset_id","retrieved_at")],
    Pfister = if(!is.null(pf)) pf$provenance else "Unavailable",
    HPA = if(!is.null(hp)) hp$metadata else "Unavailable",
    DepMap = if(!is.null(dm)) dm$provenance else "Unavailable",
    PubMed = if(!is.null(pub)) list(query=pub$query,retrieved_at=pub$retrieved_at,scope=pub$scope,sort=pub$sort) else "Unavailable")
  ids<-explorer_report_ids()
  hpa_keys<-c('tissue','brain','single_cell','subcellular','cancer','blood','cell_line','structure','interaction')
  list(gene=x$gene,main=x,pfister=pf,biology=b,hpa=hp,depmap=dm,publications=pub,cerebellum=cb,
    figure_selection=setNames(lapply(ids,function(id)!identical(input[[paste0('report_include_',id)]],FALSE)),ids),
    settings=list(survival_cutoff=or(input$survival_cutoff,"mean"),
      publication_sort=or(input$publication_sort,"relevance"),
      publication_scope=or(input$publication_scope,"general"),
      biology_target_set=or(input$biology_target_set,"curated"),depmap_view=or(input$depmap_view,'tumor'),
      cerebellum_annotation=or(input$cerebellum_annotation,'broad_lineage'),
      hpa_metrics=setNames(lapply(hpa_keys,function(k)input[[paste0('hpa_metric_',k)]]),hpa_keys)),
    generated_at=format(Sys.time(),tz="UTC",usetz=TRUE),app_version=version,git_commit=commit,sources=sources,
    availability=list(pfister=!is.null(pf),biology=!is.null(b),hpa=!is.null(hp),depmap=!is.null(dm),publications=!is.null(pub),cerebellum=!is.null(cb)))
}

start_integrated_report_worker <- function(cfg, snapshot, out_dir) {
  if(!explorer_reports_enabled()) stop('Integrated reports are disabled for this deployment.')
  snapshot_file <- file.path(out_dir,"snapshot-input.rds"); dir.create(out_dir,recursive=TRUE,showWarnings=FALSE); saveRDS(snapshot,snapshot_file)
  logfile <- file.path(out_dir,"report-build.log")
  explorer_r_bg(function(root,cfg,snapshot_file,out_dir) {
    source(file.path(root,"R","load_core.R")); core <- load_r2_core(root)
    lock<-core$deployment_worker_slot(cfg);on.exit(filelock::unlock(lock),add=TRUE)
    pdf<-tolower(Sys.getenv('MB_REPORT_PDF','true'))=='true'
    core$build_integrated_report(readRDS(snapshot_file),file.path(out_dir,"report"),root,include_pdf=pdf)
  },args=list(root=cfg$app_root,cfg=cfg,snapshot_file=snapshot_file,out_dir=out_dir),libpath=.libPaths(),stdout=logfile,stderr=logfile,supervise=TRUE)
}

explorer_report_server <- function(input,output,session,cfg,core,result,pfister,biology,hpa,depmap,publications,
  worker_start=start_integrated_report_worker,cerebellum=function()NULL) {
  if(!explorer_reports_enabled()) {
    output$report_downloads <- shiny::renderUI(shiny::p('Integrated reports are unavailable on this hosted version. Individual analysis downloads remain available.'))
    output$report_overview <- shiny::renderUI(NULL)
    return(shiny::reactiveVal(list(status='disabled')))
  }
  state <- shiny::reactiveVal(list(status="idle")); tick <- shiny::reactiveVal(0L); job <- NULL; started<-NULL
  current_snapshot<-shiny::reactive(explorer_report_snapshot(cfg,input,result,pfister,biology,hpa,depmap,publications,cerebellum))
  figure_catalog<-shiny::reactive(core$report_figure_catalog(current_snapshot()))
  revision<-shiny::reactiveVal(0L)
  shiny::observeEvent(current_snapshot(),{revision(shiny::isolate(revision())+1L)},ignoreInit=TRUE)
  built_revision<-NULL
  session$onSessionEnded(function(){if(!is.null(job)&&job$is_alive())job$kill()})
  launch <- function() {
    shiny::req(result())
    if(!is.null(job)&&job$is_alive()) return()
    snapshot <- shiny::isolate(current_snapshot());built_revision<<-shiny::isolate(revision())
    parent <- file.path(cfg$output_dir,"integrated-reports"); dir.create(parent,recursive=TRUE,showWarnings=FALSE)
    out_dir <- tempfile(paste0(core$report_safe_name(snapshot$gene),"-"),tmpdir=parent)
    state(list(status="building",gene=snapshot$gene,message="Building the integrated HTML, figures, tables and ZIP package..."))
    started<<-Sys.time()
    job <<- tryCatch(worker_start(cfg,snapshot,out_dir),error=function(e)NULL)
    if(is.null(job)) state(list(status="failed",gene=snapshot$gene,message="The report worker could not start. Please try again."))
    tick(shiny::isolate(tick())+1L)
  }
  shiny::observeEvent(input$report_generate,launch(),ignoreInit=TRUE)
  shiny::observeEvent(input$report_generate_overview,launch(),ignoreInit=TRUE)
  shiny::observeEvent(input$run,{
    if(!is.null(job)&&job$is_alive())job$kill()
    job <<- NULL; state(list(status="idle"))
  },ignoreInit=TRUE)
  shiny::observe({
    tick(); if(is.null(job))return()
    if(exists('explorer_worker_expired',mode='function')&&explorer_worker_expired(started,cfg,job)){job$kill();job<<-NULL;state(list(status='failed',message='Report generation timed out. Please retry.'));return()}
    if(job$is_alive()){shiny::invalidateLater(500,session);return()}
    built <- tryCatch(job$get_result(),error=function(e)NULL); job <<- NULL; gene <- state()$gene
    if(is.null(built)) state(list(status="failed",gene=gene,message="Report generation failed. The dashboard results remain available.")) else state(list(status="ready",gene=gene,artifacts=built,message="Report ready."))
  })
  status_ui <- function(button_id) {
    s <- state(); x <- result()
    if(is.null(x)) return(shiny::p("Run a gene analysis to enable the integrated report."))
    if(s$status=="building") return(shiny::div(class="alert alert-info",role="status",s$message))
    if(s$status=="failed") return(shiny::tagList(shiny::div(class="alert alert-warning",s$message),shiny::actionButton(button_id,"Retry Full Report",class="btn-primary")))
    if(s$status=="ready" && identical(s$gene,x$gene)&&identical(built_revision,revision())) return(shiny::tagList(shiny::div(class="alert alert-success",s$message),
      if(tolower(Sys.getenv('MB_REPORT_PDF','true'))!='true')shiny::p('Download the HTML report and use your browser\'s Print / Save as PDF for a PDF copy.'),
      shiny::div(class="d-flex flex-wrap gap-3",shiny::downloadButton(paste0(button_id,"_html"),"Download Full Report (HTML)"),
        if(!is.null(s$artifacts$pdf))shiny::downloadButton(paste0(button_id,"_pdf"),"Download Full Report (PDF)"),
        shiny::downloadButton(paste0(button_id,"_zip"),"Download Complete Report Package (ZIP)"))))
    shiny::tagList(shiny::p("Creates a frozen report from the current dashboard objects and settings. Modules still loading are recorded as unavailable."),shiny::actionButton(button_id,"Generate Full Report",class="btn-primary"))
  }
  output$report_downloads <- shiny::renderUI({shiny::tagList(shiny::hr(),shiny::h3("Complete integrated report"),status_ui("report_generate"))})
  output$report_overview <- shiny::renderUI({
    if(is.null(result()))return(status_ui('report_generate_overview'))
    s<-current_snapshot();v<-core$report_interpretation(s);o<-v$evidence$observations;figs<-Filter(function(f)f$selected,figure_catalog())
    refs<-function(ids)shiny::span(class='small text-muted',' [',lapply(ids,function(id){z<-o[[id]];shiny::a(href='#',onclick=sprintf("Shiny.setInputValue('report_evidence_click','%s',{priority:'event'});return false;",id),paste(z$source,z$tab,sep=' - '))} ),'] ')
    shiny::tagList(shiny::hr(),shiny::h3('Integrated Interpretation'),
      shiny::p(lapply(v$summary,function(z)shiny::tagList(z$text,refs(z$evidence),' '))),
      shiny::h4('Potentially Interesting Patterns'),if(!is.null(v$fallback))shiny::p(v$fallback),
      lapply(v$patterns,function(p)bslib::card(class='dashboard-card',bslib::card_header(p$title),shiny::p(p$text),shiny::p(class='small text-muted',p$strength),refs(p$evidence))),
      shiny::tags$details(shiny::tags$summary('Evidence and availability'),
        lapply(o,function(z)shiny::p(shiny::strong(z$strength),paste0(' / ',z$type,': '),z$text,refs(z$id))),
        lapply(v$evidence$missing,shiny::p)),
      shiny::h4('Figures'),shiny::p('Selected figures update with the controls beside each plot. Exports use these selections and current plot settings.'),
      if(!length(figs))shiny::p('No figures selected.')else shiny::tags$ol(lapply(figs,function(f)shiny::tags$li(
        shiny::tags$details(shiny::tags$summary(f$caption),if(is.null(f$kind))shiny::div(style='overflow-x:auto',shiny::plotOutput(paste0('report_preview_',f$id),width=if(is.null(attr(f$plot,'hpa_width')))'100%'else paste0(round(attr(f$plot,'hpa_width')*100),'px'),height=if(grepl('^hpa_',f$id))'650px'else'450px'))else plotly::plotlyOutput(paste0('report_preview_',f$id),height='450px'))))),
      status_ui('report_generate_overview'))
  })
  for(id in explorer_report_ids())local({key<-id
    if(grepl('^cerebellum_',key))output[[paste0('report_preview_',key)]]<-plotly::renderPlotly({
      s<-current_snapshot();f<-figure_catalog()[[key]];shiny::req(f$selected)
      sp<-sub('.*_','',key);a<-list(cells=s$cerebellum$cells[[sp]],display_cells=s$cerebellum$display_cells[[sp]])
      shiny::req(a$cells,a$display_cells)
      field<-s$settings$cerebellum_annotation
      assets<-list(species=lapply(s$cerebellum$cells,function(d)list(cells=d)))
      core$cerebellum_plot(a,field,if(f$kind=='expression')s$cerebellum$species[[sp]]else NULL,core$cerebellum_colors(assets,field))
    })else output[[paste0('report_preview_',key)]]<-shiny::renderPlot({f<-figure_catalog()[[key]];shiny::req(f$selected,f$plot);f$plot},res=120)
  })
    shiny::observeEvent(input$report_evidence_click,{
      z<-core$report_interpretation(current_snapshot())$evidence$observations[[input$report_evidence_click]]
      if(!is.null(z))bslib::nav_select('section',selected=z$tab,session=session)
    },ignoreInit=TRUE)
  for(prefix in c("report_generate","report_generate_overview")) local({p<-prefix
    for(kind in c("html","pdf","zip")) local({k<-kind
      output[[paste0(p,"_",k)]] <- shiny::downloadHandler(
        filename=function(){s<-state();shiny::req(s$status=="ready");basename(s$artifacts[[k]])},
        contentType=switch(k,html="text/html",pdf="application/pdf",zip="application/zip"),
        content=function(file){s<-state();shiny::req(s$status=="ready",identical(built_revision,revision()),file.exists(s$artifacts[[k]]));file.copy(s$artifacts[[k]],file,overwrite=TRUE)})
    })
  })
  state
}
