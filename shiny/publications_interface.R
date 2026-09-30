# Independent literature discovery. No expression/clinical/TF services called here.
explorer_publications_ui <- function() bslib::nav_panel("Publications",
  shiny::fluidRow(shiny::column(6,shiny::selectInput("publication_sort","Publication ranking",c("Most relevant"="relevance","Most recent"="recent"),selected="relevance",selectize=FALSE)),
    shiny::column(6,shiny::selectInput("publication_scope","Publication scope",c("Gene - all contexts"="general","Gene + medulloblastoma"="medulloblastoma"),selected="general",selectize=FALSE))),
  shiny::actionButton("publication_refresh","Refresh PubMed"),shiny::uiOutput("publication_status"),
  shiny::uiOutput("publication_query"),shiny::uiOutput("publication_cards"),
  shiny::p(class="text-muted","PubMed ranking and publication counts do not measure evidence quality, biological importance, clinical significance or consensus."),
  shiny::a(href="https://www.ncbi.nlm.nih.gov/About/disclaimer.html",target="_blank",rel="noopener","NCBI disclaimer and copyright information"))

start_publication_worker <- function(cfg,gene,sort,scope,refresh) {
  folder<-file.path(cfg$output_dir,"publication-jobs");dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  logfile<-tempfile("pubmed-",folder,fileext=".log")
  callr::r_bg(function(root,cfg,gene,sort,scope,refresh) {
    source(file.path(root,"R","load_core.R"));core<-load_r2_core(root)
    lock<-core$deployment_worker_slot(cfg);on.exit(filelock::unlock(lock),add=TRUE)
    core$get_pubmed_publications(cfg,gene,n=10L,sort=sort,scope=scope,refresh=refresh)
  },args=list(root=cfg$app_root,cfg=cfg,gene=gene,sort=sort,scope=scope,refresh=refresh),
    libpath=.libPaths(),stdout=logfile,stderr=logfile,supervise=TRUE)
}

explorer_publications_server <- function(input,output,session,cfg,core,worker_start=start_publication_worker) {
  state<-shiny::reactiveVal(list(status="idle"));tick<-shiny::reactiveVal(0L);job<-NULL;started<-NULL;current_query<-NULL
  queried_gene<-shiny::eventReactive(input$run,{toupper(trimws(input$gene))},ignoreInit=TRUE)
  launch<-function(refresh=FALSE) {
    gene<-queried_gene();shiny::req(gene)
    query<-list(gene,input$publication_sort,input$publication_scope)
    if(!is.null(job) && job$is_alive() && identical(query,current_query))return()
    current_query<<-query;started<<-Sys.time()
    if(!is.null(job) && job$is_alive())job$kill()
    job<<-NULL
    valid<-tryCatch({core$biology_symbol(gene);TRUE},error=function(e)FALSE)
    if(!valid) {state(list(status="failed",message="Enter one valid gene symbol to search PubMed."));return()}
    state(list(status="loading",message=paste("Loading PubMed publications for",gene,"in the background.")))
    job<<-tryCatch(worker_start(cfg,gene,input$publication_sort %||% "relevance",input$publication_scope %||% "general",refresh),
      error=function(e) {message("[PubMed worker] Could not start background retrieval");state(list(status="failed",message="PubMed results are currently unavailable."));NULL})
    tick(shiny::isolate(tick())+1L)
  }
  `%||%`<-core$`%||%`
  shiny::observeEvent(list(queried_gene(),input$publication_sort,input$publication_scope),launch())
  shiny::observeEvent(input$publication_refresh,launch(TRUE),ignoreInit=TRUE)
  shiny::observe({
    tick();if(is.null(job))return()
    if(exists("explorer_worker_expired",mode="function") && explorer_worker_expired(started,cfg)){
      job$kill();job<<-NULL;state(list(status="failed",message="PubMed timed out. Please retry later; other analyses remain available."));return()
    }
    if(job$is_alive()){shiny::invalidateLater(500,session);return()}
    x<-tryCatch(job$get_result(),error=function(e)NULL);job<<-NULL
    if(is.null(x))state(list(status="failed",message="PubMed results are currently unavailable. Verify the official gene symbol or try Refresh PubMed later; other analyses remain available.")) else
      state(list(status="complete",data=x))
  })
  session$onSessionEnded(function(){if(!is.null(job) && job$is_alive())job$kill()})
  session$userData$publications_loading_state<-shiny::reactive(state()$status)
  publications<-shiny::reactive({x<-state();if(x$status=="complete")x$data else NULL})
  show_value<-function(x,missing="Not supplied in PubMed")if(is.null(x) || !length(x) || is.na(x[1]) || !nzchar(x[1]))missing else x
  output$publication_status<-shiny::renderUI({
    s<-state();x<-publications()
    if(is.null(x))return(shiny::div(class="alert alert-info",if(s$status=="idle")"Run a gene query to search PubMed." else s$message))
    shiny::tagList(shiny::h3(if(x$sort=="relevance")"Top 10 PubMed results by relevance" else "10 most recent PubMed results"),
      shiny::p(paste(x$canonical_symbol,"-",x$total_hits,"matching PubMed results;",x$returned_count,"shown.")),
      shiny::p(if(x$sort=="relevance")"PubMed relevance / Best Match ordering." else "Publication-date descending (not PubMed indexing date)."),
      if(!x$returned_count)shiny::p("No PubMed results matched this query."),lapply(x$warnings,shiny::p))
  })
  output$publication_query<-shiny::renderUI({
    x<-publications();shiny::req(x)
    shiny::tags$details(shiny::tags$summary("Gene identity, exact query and retrieval details"),
      shiny::p(paste("Original input:",x$queried_gene,"| Canonical human symbol:",x$canonical_symbol,"| Official name:",x$identity$official_name)),
      shiny::a(href=paste0("https://www.ncbi.nlm.nih.gov/gene/",x$identity$gene_id),target="_blank",rel="noopener","NCBI Gene identity"),
      shiny::p(paste("Aliases used:",paste(x$identity$aliases_used,collapse=", "))),
      shiny::pre(style="white-space:pre-wrap",x$query),shiny::p(paste("Scope:",x$scope,"| API sort:",x$api_sort,"| Retrieved:",x$retrieved_at,if(x$cache_hit)"(cached)" else "")),
      shiny::p(x$limitations))
  })
  output$publication_cards<-shiny::renderUI({
    x<-publications();shiny::req(x)
    shiny::tagList(lapply(x$articles,function(a)bslib::card(
      bslib::card_header(shiny::span(paste0(a$rank,". ")),shiny::a(href=a$pubmed_url,target="_blank",rel="noopener",show_value(a$title,"Title not supplied in PubMed"))),
      shiny::p(if(length(a$authors))paste(a$authors,collapse=", ") else "Authors not supplied in PubMed"),
      shiny::p(paste(show_value(a$journal),"|",show_value(a$publication_date),"| PMID:",a$pmid)),
      shiny::p("DOI: ",if(!is.na(a$doi_url))shiny::a(href=a$doi_url,target="_blank",rel="noopener",a$doi) else "Not supplied in PubMed"),
      shiny::p(paste("Publication type:",if(length(a$publication_types))paste(a$publication_types,collapse=", ") else "Not supplied")),
      if(a$retracted)shiny::p(class="alert alert-warning","PubMed flags this as a retracted publication."),
      if(length(a$notices))shiny::p(paste("PubMed notices:",paste(a$notices,collapse=", "))),
      shiny::tags$details(shiny::tags$summary("Abstract from PubMed"),shiny::div(style="white-space:pre-wrap",show_value(a$abstract,"No abstract available in PubMed."))))))
  })
  output$publications_overview<-shiny::renderUI({
    x<-publications();if(is.null(x))return(shiny::p(if(state()$status=="failed")"PubMed results are currently unavailable." else "Publications are pending; see the Publications tab."))
    shiny::tagList(shiny::p(paste(x$canonical_symbol,"|",x$scope,"|",x$total_hits,"matching PubMed publications |",x$sort,"ranking")),
      if(length(x$articles))shiny::p(if(x$sort=="relevance")"Top result by PubMed relevance: " else "First result by publication date: ",
        shiny::a(href=x$articles[[1]]$pubmed_url,target="_blank",rel="noopener",show_value(x$articles[[1]]$title))),
      if(length(x$articles))shiny::p(paste("First result's publication date:",show_value(x$articles[[1]]$publication_date))),
      shiny::p("Counts and ranking do not establish scientific consensus."))
  })
  output$publication_downloads<-shiny::renderUI({shiny::req(publications());shiny::tagList(shiny::h4("Publications (current scope and ranking)"),
    shiny::downloadButton("publication_csv","PubMed CSV"),shiny::downloadButton("publication_json","PubMed JSON"))})
  for(format in c("csv","json"))local({f<-format
    output[[paste0("publication_",f)]]<-shiny::downloadHandler(filename=function(){shiny::req(publications());core$publication_download_name(publications(),f)},
      content=function(file){shiny::req(publications());core$write_publication_download(publications(),file,f)})
  })
  publications
}
