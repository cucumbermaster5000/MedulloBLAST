explorer_depmap_ui <- function() bslib::nav_panel("DepMap",
  shiny::uiOutput("depmap_status"),shiny::uiOutput("depmap_dataset"),
  shiny::h3("Medulloblastoma and retinoblastoma models"),shiny::uiOutput("depmap_counts"),
  shiny::selectInput("depmap_view","Dependency view",c("Medulloblastoma and retinoblastoma"="tumor","All DepMap models with both highlighted"="all"),selected="tumor",selectize=FALSE),
  shiny::plotOutput("depmap_dependency_plot",height=explorer_plot_height("large")),explorer_report_control('depmap_dependency'),shiny::uiOutput("depmap_rank_definition"),
  shiny::textInput("depmap_search","Search tumor models or ModelIDs"),
  shiny::selectInput("depmap_sort","Sort model table",c("Strongest dependency first"="dependency","Model name"="model","Highest relative percentile first"="percentile"),selectize=FALSE),
  shiny::div(style="overflow:auto",shiny::tableOutput("depmap_mb_table")),
  shiny::h4("Descriptive statistics"),shiny::tableOutput("depmap_statistics"),shiny::uiOutput("depmap_extremes"),
  shiny::h3("Expression vs dependency"),shiny::uiOutput("depmap_correlation"),shiny::plotOutput("depmap_expression_plot",height=explorer_plot_height()),explorer_report_control('depmap_expression_dependency'),
  shiny::h3("Model metadata and subgroup provenance"),shiny::uiOutput("depmap_model_filter"),
  shiny::div(style="overflow:auto",shiny::tableOutput("depmap_metadata")),
  shiny::p("Cell-line dependency is not validated therapeutic efficacy. Relative ranks do not establish MB-specific dependency. Expression is not dependency."))

start_depmap_worker <- function(cfg,gene) {
  folder<-file.path(cfg$output_dir,"depmap-jobs");dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  logfile<-tempfile("depmap-",folder,fileext=".log")
  explorer_r_bg(function(root,cfg,gene){source(file.path(root,"R","load_core.R"));core<-load_r2_core(root);lock<-core$deployment_worker_slot(cfg);on.exit(filelock::unlock(lock),add=TRUE);core$get_depmap_profile(cfg,gene)},
    args=list(root=cfg$app_root,cfg=cfg,gene=gene),libpath=.libPaths(),stdout=logfile,stderr=logfile,supervise=TRUE)
}

explorer_depmap_server <- function(input,output,session,cfg,core,worker_start=start_depmap_worker) {
  `%||%`<-core$`%||%`;state<-shiny::reactiveVal(list(status="idle"));tick<-shiny::reactiveVal(0L);job<-NULL;started<-NULL;current_gene<-NULL
  shiny::observeEvent(input$run,{
    gene<-tryCatch(core$biology_symbol(input$gene),error=function(e)input$gene)
    if(!is.null(job) && job$is_alive() && identical(gene,current_gene))return()
    if(!is.null(job) && job$is_alive())job$kill()
    job<<-NULL
    current_gene<<-gene;started<<-Sys.time()
    if(!tryCatch({core$biology_symbol(gene);TRUE},error=function(e)FALSE)){state(list(status="failed",message="Enter one valid gene symbol for DepMap."));return()}
    state(list(status="loading",message=paste("Loading current same-release DepMap data for",gene,"from the official public service.")))
    job<<-tryCatch(worker_start(cfg,gene),error=function(e){state(list(status="failed",message="DepMap could not start. Other analyses remain available."));NULL})
    tick(shiny::isolate(tick())+1L)
  },ignoreInit=TRUE)
  shiny::observe({
    tick();if(is.null(job))return()
    if(exists("explorer_worker_expired",mode="function") && explorer_worker_expired(started,cfg,job)){
      job$kill();job<<-NULL;state(list(status="failed",message="DepMap timed out. Please retry later; other analyses remain available."));return()
    }
    if(job$is_alive()){shiny::invalidateLater(500,session);return()}
    x<-tryCatch(job$get_result(),error=function(e){message("[DepMap worker] ",conditionMessage(e));NULL});job<<-NULL
    if(is.null(x))state(list(status="failed",message="DepMap is currently unavailable. Please retry later; other analyses remain available.")) else state(list(status="complete",data=x))
  })
  session$onSessionEnded(function(){if(!is.null(job) && job$is_alive())job$kill()})
  session$userData$depmap_loading_state<-shiny::reactive(state()$status)
  data<-shiny::reactive({s<-state();if(s$status=="complete")s$data else NULL})
  output$depmap_status<-shiny::renderUI({s<-state();shiny::div(class="alert alert-info",if(s$status=="idle")"Run a gene query to load DepMap." else if(s$status=="complete")paste("DepMap loaded for",s$data$canonical_symbol) else s$message)})
  output$depmap_dataset<-shiny::renderUI({x<-data();shiny::req(x);p<-x$provenance
    shiny::tagList(shiny::h3(p$release),shiny::p(p$release_note),
      shiny::p(paste("Dependency:",p$files$dependency$name,"| Metric:",p$dependency_metric)),
      shiny::p(paste("Expression:",p$files$expression$name %||% "Unavailable","| Metric:",p$expression_metric)),
      shiny::p(paste("Model metadata:",p$files$model$name,"| Files verified:",p$verified_at,"| Analysis retrieved:",x$retrieved_at)),
      shiny::p(paste("Original query:",x$queried_gene,"| Resolved symbol:",x$canonical_symbol,"| Entrez ID:",x$identity$entrez_id,"| Mapping:",x$identity$source)),
      shiny::p(paste("Dependency feature:",x$mapping$dependency_column,"| Expression feature:",x$mapping$expression_column)),
      shiny::p("More negative Chronos scores indicate stronger dependency. Zero is the nonessential-gene reference; no dependent/nondependent cutoff is applied."),lapply(x$warnings,shiny::p))
  })
  output$depmap_counts<-shiny::renderUI({x<-data();shiny::req(x);c<-x$counts
    shiny::p(sprintf("Medulloblastoma: %d models | MB dependency: %d | Retinoblastoma: %d models | Retinoblastoma dependency: %d | All models with valid dependency: %d",c$mb_models,c$mb_dependency,c$retinoblastoma_models,c$retinoblastoma_dependency,c$all_dependency))})
  output$depmap_dependency_plot<-shiny::renderPlot({x<-data();shiny::req(x);p<-if(identical(input$depmap_view,"all"))x$plots$all_models_ranked_dependency else x$plots$dependency;shiny::validate(shiny::need(!is.null(p),"No valid dependency data for this view."));explorer_web_plot(p,"depmap_dependency_plot",session)},res=120,execOnResize=TRUE)
  output$depmap_rank_definition<-shiny::renderUI({x<-data();shiny::req(x);shiny::tagList(shiny::p(x$rank_definition),shiny::p(x$percentile_definition))})
  output$depmap_mb_table<-shiny::renderTable({x<-data();shiny::req(x);d<-core$depmap_tumor_display(x)
    q<-tolower(input$depmap_search %||% "");if(nzchar(q))d<-d[grepl(q,tolower(paste(d$ModelID,d$CellLineName)),fixed=TRUE),,drop=FALSE]
    order<-input$depmap_sort %||% "dependency";d<-d[if(order=="model")order(d$CellLineName) else if(order=="percentile")order(-d$dependency_percentile,d$ModelID,na.last=TRUE) else order(d$dependency,d$ModelID,na.last=TRUE),,drop=FALSE]
    d
  },digits=4,striped=TRUE,na="Missing")
  output$depmap_statistics<-shiny::renderTable({shiny::req(data());data()$statistics},digits=4)
  output$depmap_extremes<-shiny::renderUI({
    x<-data();shiny::req(x)
    shiny::tagList(lapply(c("strongest","weakest"),function(k){
      d<-x[[k]]
      if(nrow(d))shiny::p(paste(k,"observed MB dependency:",paste(d$CellLineName,collapse=", "),"| Chronos score:",format(d$dependency[1],digits=4)))
    }),if(x$counts$mb_verified_subgroup==0)shiny::p("Insufficient models for reliable subgroup comparison: no verified molecular subgroup assignments in this release."))
  })
  output$depmap_correlation<-shiny::renderUI({shiny::req(data());a<-data()$correlation;shiny::tagList(shiny::p(paste("N =",a$n,"| Spearman rho =",format(a$rho,digits=4),"| p =",format(a$p_value,digits=4))),shiny::p(a$status),shiny::p(a$method))})
  output$depmap_expression_plot<-shiny::renderPlot({shiny::req(data());p<-data()$plots$expression_dependency;shiny::validate(shiny::need(!is.null(p),"No matched expression/dependency values."));explorer_web_plot(p,"depmap_expression_plot",session)},res=120,execOnResize=TRUE)
  output$depmap_model_filter<-shiny::renderUI({shiny::req(data());shiny::tagList(shiny::p(data()$mb_filter),shiny::p(data()$retinoblastoma_filter))})
  output$depmap_metadata<-shiny::renderTable({shiny::req(data());d<-data()$tumor_models;d[,intersect(c("CellLineName","ModelID","cancer_type","OncotreeLineage","OncotreePrimaryDisease","OncotreeSubtype","OncotreeCode","mb_inclusion_reason","retinoblastoma_inclusion_reason","mb_subgroup","subgroup_source","subgroup_reference","dependency_missing","expression_missing"),names(d)),drop=FALSE]},striped=TRUE)
  output$depmap_overview<-shiny::renderUI({x<-data();if(is.null(x))return(shiny::p(if(state()$status=="failed")"DepMap is currently unavailable." else "DepMap is pending; see the DepMap tab."))
    s<-x$statistics[grepl("^Medulloblastoma",x$statistics$cohort),];strong<-x$strongest
    shiny::tagList(shiny::p(paste(x$canonical_symbol,"|",x$provenance$release,"| Medulloblastoma models:",x$counts$mb_models,"| Retinoblastoma models:",x$counts$retinoblastoma_models)),
      shiny::p(paste("Median MB dependency:",format(s$median,digits=4))),
      if(nrow(strong))shiny::p(paste("Strongest observed MB dependency:",paste(strong$CellLineName,collapse=", "),"| score:",format(strong$dependency[1],digits=4),"| best rank:",strong$dependency_rank[1],"of",strong$total_models[1],"| strictly-weaker percentile:",format(strong$dependency_percentile[1],digits=4))),
      shiny::p(paste("Expression-dependency correlation: N =",x$correlation$n,"| rho =",format(x$correlation$rho,digits=3),"| p =",format(x$correlation$p_value,digits=3))),shiny::p("Cell-line dependency and relative rank do not establish therapeutic efficacy or MB specificity."))
  })
  kinds<-c("dependency","expression","matched_expression_dependency","model_metadata","subgroup_evidence","all_models_dependency","MB_context_ranks","tumor_context_ranks","statistics","result.rds",paste0(rep(c("dependency","expression_dependency","all_models_ranked_dependency"),each=2),"_plot.",c("pdf","png")))
  for(kind in kinds)local({k<-kind;output[[paste0("depmap_download_",gsub("[.]","_",k))]]<-shiny::downloadHandler(
    filename=function(){shiny::req(data());paste0(data()$canonical_symbol,"_DepMap_",k,if(grepl("[.](pdf|png|rds)$",k))"" else ".csv")},
    content=function(file){shiny::req(data());core$write_depmap_download(data(),k,file)})})
  output$depmap_downloads<-shiny::renderUI({x<-data();shiny::req(x)
    available<-c("statistics","result.rds",if(nrow(x$tumor_models))c("model_metadata","subgroup_evidence"),if(nrow(core$depmap_tumor_display(x)))"dependency",if(any(is.finite(x$tumor_models$expression)))"expression",if(nrow(x$matched))"matched_expression_dependency",if(nrow(x$all_models))"all_models_dependency",if(nrow(x$mb_context))"MB_context_ranks",if(nrow(x$tumor_context))"tumor_context_ranks",unlist(lapply(names(x$plots),function(p)paste0(p,"_plot.",c("pdf","png")))))
    shiny::tagList(shiny::h4("DepMap downloads"),shiny::div(class="d-flex flex-wrap gap-3",lapply(available,function(k)shiny::downloadButton(paste0("depmap_download_",gsub("[.]","_",k)),gsub("_"," ",k)))))
  })
  data
}
