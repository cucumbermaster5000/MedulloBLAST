# Presentation and background orchestration only; scientific functions are in R/.
explorer_biology_ui <- function() {
  categories<-c(GO_BP="Biological processes",GO_MF="Molecular functions",GO_CC="Cellular components",Reactome="Reactome",KEGG="KEGG",WikiPathways="WikiPathways")
  list(functional=bslib::nav_panel("Functional Biology",
    shiny::uiOutput("biology_status"),shiny::actionButton("biology_refresh","Refresh external biology resources"),
    shiny::h3("A. Direct gene annotation / pathway membership"),shiny::uiOutput("functional_summary"),
    shiny::plotOutput("annotation_plot",height=explorer_plot_height("small")),explorer_report_control('functional_annotation'),
    lapply(names(categories),function(cat)shiny::tagList(shiny::h4(categories[[cat]]),
      shiny::div(style="overflow:auto;max-height:350px",shiny::tableOutput(paste0("annotation_",cat))))),
    shiny::h3("B. Target-gene programs"),
    shiny::selectInput("biology_target_set","TF target set",choices=c("Curated human regulatory targets (default)"="curated","Human binding-associated targets"="binding"),selected="curated"),
    shiny::uiOutput("enrichment_summary"),shiny::uiOutput("enrichment_plot_ui"),explorer_report_control('tf_target_enrichment'),
    shiny::numericInput("enrichment_page","Enrichment results page (100 rows per page)",1,min=1,step=1),
    shiny::div(style="overflow:auto;max-height:600px",shiny::tableOutput("target_enrichment_table"))),
    targets=bslib::nav_panel("TF Targets",shiny::uiOutput("tf_classification"),shiny::uiOutput("tf_evidence_summary"),
      shiny::textInput("tf_search","Filter target gene or evidence text",value=""),
      shiny::numericInput("tf_page","Evidence page (100 rows per page)",1,min=1,step=1),
      shiny::uiOutput("tf_page_summary"),
      shiny::div(style="overflow:auto;max-height:650px",shiny::tableOutput("tf_evidence_table"))))
}

explorer_biology_server <- function(input,output,session,cfg,result,core) {
  `%||%` <- core$`%||%`
  state<-shiny::reactiveVal(list(status="idle")); job<-NULL; current_gene<-NULL; tick<-shiny::reactiveVal(0L);started<-NULL
  session$userData$biology_loading_state<-shiny::reactive(state()$status)
  queried_gene<-shiny::eventReactive(input$run,{tryCatch(core$biology_symbol(input$gene),error=function(e)NULL)},ignoreInit=TRUE)
  launch<-function(gene,refresh=FALSE) {
    if(!is.null(job) && job$is_alive() && identical(gene,current_gene))return()
    started<<-Sys.time()
    if(!is.null(job) && job$is_alive())job$kill()
    current_gene<<-gene
    if(!requireNamespace("callr",quietly=TRUE)) {state(list(status="failed",gene=gene,message="Background support is unavailable; install the Shiny dependencies."));return()}
    root<-cfg$app_root %||% normalizePath(".",winslash="/")
    folder<-file.path(cfg$output_dir %||% file.path(root,"outputs"),"biology-jobs");dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    logfile<-tempfile(paste0(gene,"-"),folder,fileext=".log")
    state(list(status="running",gene=gene,message="Loading documented biology in the background. Other tabs remain available."))
    job<<-tryCatch(explorer_r_bg(function(root,cfg,gene,refresh) {
      source(file.path(root,"R","load_core.R"))
      core<-load_r2_core(root)
      lock<-core$deployment_worker_slot(cfg);on.exit(filelock::unlock(lock),add=TRUE)
      core$get_gene_biology(cfg,gene,refresh)
    },args=list(root=root,cfg=cfg,gene=gene,refresh=refresh),libpath=.libPaths(),
      stdout=logfile,stderr=logfile,supervise=TRUE),error=function(e) {
        message("[Biology background] ",conditionMessage(e));state(list(status="failed",gene=gene,message="Functional biology could not be started. Existing results remain available."));NULL
      })
    tick(shiny::isolate(tick())+1L)
  }
  shiny::observeEvent(queried_gene(),{
    gene<-queried_gene()
    if(is.null(gene)) {
      if(!is.null(job) && job$is_alive())job$kill()
      job<<-NULL;current_gene<<-NULL;state(list(status="idle"));return()
    }
    if(!identical(gene,current_gene))launch(gene)
  },ignoreNULL=FALSE)
  shiny::observeEvent(input$biology_refresh,{shiny::req(queried_gene());launch(queried_gene(),TRUE)},ignoreInit=TRUE)
  shiny::observe({
    tick()
    if(is.null(job))return()
    if(exists("explorer_worker_expired",mode="function") && explorer_worker_expired(started,cfg,job)){
      job$kill();job<<-NULL;state(list(status="failed",gene=current_gene,message="Functional biology timed out. Please retry later."));return()
    }
    if(job$is_alive()) {shiny::invalidateLater(500,session);return()}
    x<-tryCatch(job$get_result(),error=function(e) {message("[Biology worker] ",conditionMessage(e));NULL})
    job<<-NULL
    if(is.null(x))state(list(status="failed",gene=current_gene,message="Functional biology is currently unavailable. Existing analysis tabs remain available.")) else state(list(status="complete",gene=current_gene,data=x))
  })
  session$onSessionEnded(function(){if(!is.null(job) && job$is_alive())job$kill()})
  biology<-shiny::reactive({s<-state();if(s$status!="complete" || !identical(s$gene,queried_gene()))NULL else s$data})
  output$biology_status<-shiny::renderUI({
    s<-state()
    shiny::div(class=if(s$status=="failed")"alert alert-warning" else "alert alert-info",
      if(s$status=="idle")"Run a gene analysis to load biological evidence." else if(s$status=="complete")paste("Biological evidence loaded for",s$gene) else s$message)
  })
  output$functional_summary<-shiny::renderUI({
    x<-biology();if(is.null(x))return(shiny::p("Biological evidence has not finished loading."))
    shiny::tagList(shiny::h4(x$gene),shiny::p(x$annotation$summary),lapply(x$annotation$warnings,shiny::p),
      shiny::p("Library memberships do not include per-annotation experimental evidence codes. They must not all be described as experimentally proven."))
  })
  output$annotation_plot<-shiny::renderPlot({shiny::req(biology());explorer_web_plot(biology()$annotation$plot,"annotation_plot",session)},res=120,execOnResize=TRUE)
  for(category in names(core$biology_libraries()))local({
    cat<-category
    output[[paste0("annotation_",cat)]]<-shiny::renderTable({
      shiny::req(biology());d<-biology()$annotation$rows
      d[d$category==cat,c("term","term_identifier","evidence_type","library","retrieved_at"),drop=FALSE]
    },striped=TRUE)
  })
  output$tf_classification<-shiny::renderUI({
    x<-biology();if(is.null(x))return(shiny::p("TF classification is loading or unavailable; check Functional Biology for status."))
    shiny::tagList(shiny::h3(paste(x$gene,"TF classification:",x$tf$status)),shiny::p(paste("Family/domain:",x$tf$family)),
      shiny::p(x$tf$source),shiny::a(href=x$tf$source_url,target="_blank","TF reference catalogue"),shiny::p(x$tf$note),
      if(x$tf$status!="TRUE")shiny::p(x$targets$warnings))
  })
  output$tf_evidence_summary<-shiny::renderUI({
    shiny::req(biology());x<-biology();if(x$tf$status!="TRUE")return(NULL)
    s<-x$targets$summary
    shiny::tagList(shiny::p(sprintf("Unique target identifiers from sources: %d | Binding-associated: %d | Curated regulatory: %d | Multiple sources: %d",
      nrow(s),sum(s$binding),sum(s$curated),sum(s$source_count>1))),
      shiny::p("Source counts precede symbol validation. Invalid identifiers are excluded from enrichment and recorded in the downloadable input audit."),
      shiny::p("Perturbation/prediction evidence sources are not included; their absence here is not evidence that such targets do not exist."),
      lapply(x$targets$warnings,shiny::p),shiny::p("Evidence records retain each experiment/reference. Source counts are descriptive, not an independent confidence score."))
  })
  evidence_filtered<-shiny::reactive({
    shiny::req(biology());d<-biology()$targets$evidence
    if(nrow(d) && nzchar(input$tf_search %||% "")) {
      query<-tolower(input$tf_search)
      d<-d[grepl(query,tolower(paste(d$target_gene,d$source,d$experiment,d$direction)),fixed=TRUE),,drop=FALSE]
    }
    d[order(d$target_gene,d$source,d$experiment),,drop=FALSE]
  })
  page_rows<-function(d,page) {
    page<-as.integer(page %||% 1L)
    if(length(page)!=1L || is.na(page) || page<1L)page<-1L
    head(d[seq_len(nrow(d))>(page-1L)*100L,,drop=FALSE],100)
  }
  output$tf_page_summary<-shiny::renderUI({shiny::p(paste(nrow(evidence_filtered()),"matching evidence records. Complete evidence is downloadable."))})
  output$tf_evidence_table<-shiny::renderTable({
    d<-page_rows(evidence_filtered(),input$tf_page)
    s<-biology()$targets$summary
    if(nrow(d))d$source_support_count<-s$source_count[match(d$target_gene,s$target_gene)]
    d
  },striped=TRUE)
  selected_enrichment<-shiny::reactive({shiny::req(biology());biology()$enrichment[[input$biology_target_set %||% "curated"]]})
  output$enrichment_summary<-shiny::renderUI({
    shiny::req(biology());a<-selected_enrichment()
    if(is.null(a))return(shiny::p("TF-target enrichment is skipped because this gene is not classified as a TF or TF status is unknown."))
    shiny::tagList(shiny::p(a$definition),shiny::p(paste("Target-set size:",length(a$submitted_genes),"| Service-retained genes:",length(a$accepted_genes),"| Returned terms:",nrow(a$rows))),
      shiny::p(a$method),lapply(a$warnings,shiny::p),shiny::p("This is enrichment of a supported target set, not evidence that the queried TF activates these pathways in medulloblastoma."))
  })
  output$enrichment_plot<-shiny::renderPlot({a<-selected_enrichment();shiny::req(a);explorer_web_plot(a$plot,"enrichment_plot",session)},res=120,execOnResize=TRUE)
  output$enrichment_plot_ui<-shiny::renderUI({a<-selected_enrichment();if(is.null(a$plot))shiny::p("No enrichment plot is available for this target set; see the evidence summary.") else shiny::plotOutput("enrichment_plot",height=explorer_plot_height("pathways"))})
  output$target_enrichment_table<-shiny::renderTable({
    a<-selected_enrichment();shiny::req(a);d<-page_rows(a$rows,input$enrichment_page)
    for(f in intersect(c("p_value","adjusted_p_value"),names(d)))d[[f]]<-format(d[[f]],scientific=TRUE,digits=4)
    d
  },striped=TRUE,digits=4)
  annotation_kinds<-c("functional_annotation",paste0(names(core$biology_libraries()),"_annotation"))
  target_kinds<-c("TF_targets","TF_target_evidence","TF_target_summary","TF_target_input_audit","TF_target_library_recognition",
    paste0("TF_target_enrichment_",c("GO_BP","GO_MF","Reactome","KEGG","WikiPathways")),"TF_target_enrichment_plot.pdf","TF_target_enrichment_plot.png")
  all_kinds<-c(annotation_kinds,target_kinds)
  for(kind in all_kinds)local({
    k<-kind;id<-paste0("bio_download_",gsub("[.]","_",k))
    output[[id]]<-shiny::downloadHandler(filename=function() {
      shiny::req(biology());paste0(biology()$gene,"_",if(startsWith(k,"TF_target_enrichment") || k %in% c("TF_target_input_audit","TF_target_library_recognition"))paste0(input$biology_target_set,"_") else "",k,if(grepl("[.](pdf|png)$",k))"" else ".csv")
    },content=function(file){shiny::req(biology());core$write_biology_download(biology(),k,file,input$biology_target_set %||% "curated")})
  })
  output$biology_downloads<-shiny::renderUI({
    shiny::req(biology());x<-biology();a<-selected_enrichment()
    kinds<-annotation_kinds
    if(x$tf$status=="TRUE")kinds<-c(kinds,target_kinds[!grepl("plot[.]",target_kinds)])
    if(!is.null(a$plot))kinds<-c(kinds,target_kinds[grepl("plot[.]",target_kinds)])
    shiny::tagList(shiny::h4("Functional and regulatory evidence"),shiny::p("Target-enrichment downloads use the target set selected in Functional Biology."),
      shiny::div(class="d-flex flex-wrap gap-3",lapply(kinds,function(k)shiny::downloadButton(paste0("bio_download_",gsub("[.]","_",k)),gsub("_"," ",k)))),
      shiny::downloadButton("biology_rds","Complete biological result (RDS)"))
  })
  output$biology_rds<-shiny::downloadHandler(filename=function()paste0(biology()$gene,"_biology.rds"),content=function(file){shiny::req(biology());saveRDS(biology(),file)})
  biology
}
