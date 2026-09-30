hpa_sections <- c(summary="Summary",tissue="Tissue",brain="Brain",single_cell="Single Cell",
  subcellular="Subcell",cancer="Cancer",blood="Blood",cell_line="Cell Line",structure="Structure",interaction="Interaction")
hpa_file_section <- function(k) switch(k,tissue="tissue_expression",brain="brain_expression",single_cell="single_cell_expression",interaction="interactions",k)
explorer_hpa_ui <- function() bslib::nav_panel("ProteinAtlas.org Profile",
  shiny::uiOutput("hpa_status"),
  do.call(bslib::navset_tab,c(list(id="hpa_section"),lapply(names(hpa_sections),function(key)
    bslib::nav_panel(hpa_sections[[key]],value=key,shiny::uiOutput(paste0("hpa_content_",key)))))))
start_hpa_worker <- function(cfg,gene) {
  callr::r_bg(function(cfg,gene){
    source(file.path(cfg$app_root,"R","load_core.R"));core<-load_r2_core(cfg$app_root)
    lock<-core$deployment_worker_slot(cfg);on.exit(filelock::unlock(lock),add=TRUE)
    core$get_hpa_profile(cfg,gene)
  },args=list(cfg=cfg,gene=gene),libpath=.libPaths(),supervise=TRUE)
}
explorer_hpa_server <- function(input,output,session,cfg,core,biology,worker_start=start_hpa_worker) {
  state<-shiny::reactiveVal(NULL);busy<-shiny::reactiveVal(FALSE);tick<-shiny::reactiveVal(0L);job<-NULL;started<-NULL
  shiny::observeEvent(input$run,{
    if(!is.null(job)&&job$is_alive())job$kill()
    job<<-NULL;state(NULL);busy(FALSE)
    gene<-tryCatch(core$biology_symbol(input$gene),error=function(e)NULL)
    if(is.null(gene))return()
    busy(TRUE);started<<-Sys.time()
    job<<-tryCatch(worker_start(cfg,input$gene),error=function(e)NULL)
    if(is.null(job)){busy(FALSE);state(list(status="unavailable"))}
    tick(shiny::isolate(tick())+1L)
  },ignoreInit=TRUE)
  shiny::observe({
    tick();if(is.null(job))return()
    if(explorer_worker_expired(started,cfg)){job$kill();job<<-NULL;busy(FALSE);state(list(status="unavailable"));return()}
    if(job$is_alive()){shiny::invalidateLater(500,session);return()}
    x<-tryCatch(job$get_result(),error=function(e)list(status="unavailable"));job<<-NULL;state(x);busy(FALSE)
  })
  session$onSessionEnded(function(){if(!is.null(job)&&job$is_alive())job$kill()})
  session$userData$hpa_loading_state<-shiny::reactive(if(busy())'loading'else if(is.null(state()))'idle'else state()$status)
  available<-shiny::reactive({x<-state();if(!is.null(x)&&identical(x$status,"complete"))x else NULL})
  link<-function(url,label)shiny::a(href=url,target="_blank",rel="noopener",label)
  output$hpa_status<-shiny::renderUI({
    x<-available()
    if(is.null(x))return(shiny::p(class="alert alert-info",if(busy())"Loading Human Protein Atlas profile..." else if(is.null(state()))"Run an analysis to view Human Protein Atlas." else "Human Protein Atlas data are currently unavailable."))
    shiny::tagList(shiny::h3(paste(x$identity$canonical_symbol,"- Human Protein Atlas")),
      link(x$links$profile,"View full profile on Human Protein Atlas"),
      shiny::p(class="text-muted small",paste("Human Protein Atlas | Entry release",x$metadata$release,"| Retrieved",x$metadata$retrieved_at,if(x$metadata$cache_hit)"(cached)" else "")),
      shiny::p(class="small","Independent research application. HPA expression and prognosis are separate from R2 MB cohorts and DepMap dependency."))
  })
  output$gene_glance<-shiny::renderUI({
    x<-available();fallback<-state()$identity
    if(is.null(x)&&is.null(fallback))return(bslib::card(class="dashboard-card gene-glance",bslib::card_header("Gene at a Glance"),
      shiny::p(if(busy())"Retrieving source-backed gene information..." else if(is.null(state()))"Enter a gene symbol and run an analysis to explore its biology." else "General gene information is currently unavailable.")))
    id<-if(!is.null(x))x$identity else fallback
    record<-if(!is.null(x))x$raw$json else list()
    value<-function(field)core$hpa_value(record[[field]])
    facts<-c("Biological function"=value("Molecular function"),"Biological processes"=value("Biological process"),
      "Protein class"=value("Protein class"),"Subcellular location"=value("Subcellular location"),
      "Tissue context"=value("RNA tissue specificity"),"Cancer / disease annotation"=value("Disease involvement"))
    b<-biology()
    if(!is.null(b)&&identical(b$gene,id$canonical_symbol))facts<-c(facts,"Transcription factor status"=if(isTRUE(b$tf$status)||identical(b$tf$status,"TRUE"))"Listed in the TF catalogue" else if(identical(b$tf$status,"FALSE"))"Not listed in the TF catalogue" else as.character(b$tf$status),"Family / domain"=b$tf$family)
    else if(grepl("Transcription factors",value("Protein class"),fixed=TRUE))facts<-c(facts,"Transcription factor status"="Listed as a transcription factor by HPA")
    facts<-facts[!is.na(facts)&nzchar(facts)]
    orientation<-if(nzchar(value("Molecular function")))paste0("HPA lists the molecular function(s): ",value("Molecular function"),".") else if(nzchar(value("Biological process")))paste0("HPA annotates this gene with the biological processes: ",value("Biological process"),".") else ""
    bslib::card(class="dashboard-card gene-glance",bslib::card_header("Gene at a Glance"),
      shiny::h2(id$canonical_symbol),shiny::h4(id$official_name),shiny::p(orientation),
      shiny::div(class="gene-facts",lapply(names(facts),function(n)shiny::div(shiny::strong(n),shiny::p(facts[[n]])))),
      shiny::tags$details(shiny::tags$summary("Aliases and identifiers"),
        shiny::p(paste("Original query:",id$original_query)),shiny::p(paste("Aliases:",paste(id$aliases,collapse=", "))),
        shiny::p(paste("Ensembl:",id$ensembl_id,"| Entrez:",id$gene_id,"| UniProt:",paste(id$uniprot_id,collapse=", ")))),
      shiny::div(class="source-badges",if(!is.null(x))link(x$links$profile,"HPA"),
        if(!is.null(id$gene_id))link(paste0("https://www.ncbi.nlm.nih.gov/gene/",id$gene_id),"NCBI"),
        if(length(id$uniprot_id))link(paste0("https://www.uniprot.org/uniprotkb/",id$uniprot_id[1]),"UniProt")))
  })
  output$hpa_overview<-shiny::renderUI({x<-available()
    if(is.null(x))return(shiny::p(if(busy())"HPA is loading." else "Human Protein Atlas data are currently unavailable."))
    shiny::p(paste("Human Protein Atlas release",x$metadata$release,"|",core$hpa_value(x$raw$json[["RNA tissue specificity"]]),"|",core$hpa_value(x$raw$json[["Subcellular location"]])),". See ProteinAtlas.org Profile for expression, annotations and source provenance.")
  })
  for(key in names(hpa_sections))local({k<-key
    output[[paste0("hpa_content_",k)]]<-shiny::renderUI({
      x<-available();if(is.null(x))return(NULL)
      tab<-core$hpa_section_table(x,k)
      if(!nrow(tab))return(shiny::p("No data available from Human Protein Atlas for this section."))
      s<-x[[k]];d<-if(k!="summary")s$expression else NULL
      shiny::tagList(
        if(k=="brain")shiny::p(class="alert alert-info","Human adult brain context; these data do not establish developmental cerebellar expression."),
        if(k=="cancer")shiny::p(class="alert alert-info","HPA cohort prognostic annotations are not R2 medulloblastoma survival results. HPA's reported prognostic labels are preserved."),
        if(k=="single_cell")shiny::p("Specificity-selected cell types from the per-gene record; this is not a complete single-cell expression matrix."),
        if(k=="cell_line")shiny::p("HPA cell-line expression/context, separate from DepMap dependency. Individual expression values may not be supplied in the per-gene record."),
        if(k=="interaction")shiny::p("The per-gene record supplies an interaction count, when available. Partner-level evidence is not supplied here; see the original HPA profile. Interaction evidence does not establish dependency."),
        if(!is.null(d)&&nrow(d))shiny::tagList(shiny::uiOutput(paste0("hpa_choices_",k)),shiny::uiOutput(paste0("hpa_chart_",k)),explorer_report_control(paste0('hpa_',k),input)),
        if(k!="summary")shiny::selectInput(paste0("hpa_table_kind_",k),"Table content",
          choices=setNames(names(Filter(function(z)is.data.frame(z)&&nrow(z)>0,s[c("annotations","expression","details")])),names(Filter(function(z)is.data.frame(z)&&nrow(z)>0,s[c("annotations","expression","details")])))),
        shiny::textInput(paste0("hpa_filter_",k),"Filter table",placeholder="Filter any field"),
        shiny::div(style="max-height:600px;overflow:auto",shiny::tableOutput(paste0("hpa_table_",k))),shiny::downloadButton(paste0("hpa_csv_",k),"Download complete section CSV"))
    })
    output[[paste0("hpa_choices_",k)]]<-shiny::renderUI({
      d<-available()[[k]]$expression;shiny::req(nrow(d)>0)
      choices<-unique(paste(d$assay,d$metric,sep=" | "))
      preferred<-if(k=="brain")"humanBrainRegional | nTPM" else if(k=="tissue")"consensusTissue | nTPM" else choices[1]
      shiny::selectInput(paste0("hpa_metric_",k),"Assay and expression metric",choices=choices,selected=if(preferred %in% choices)preferred else choices[1],selectize=FALSE)
    })
    output[[paste0("hpa_chart_",k)]]<-shiny::renderUI({
      shiny::req(available(),input[[paste0("hpa_metric_",k)]])
      p<-core$report_hpa_selected_plot(available(),k,input[[paste0("hpa_metric_",k)]])
      shiny::req(p)
      shiny::tagList(shiny::p(class="text-muted small","All available values for this assay are shown. Scroll horizontally to see every label on wide charts."),
        shiny::div(style="overflow-x:auto",shiny::plotOutput(paste0("hpa_plot_",k),
          width=paste0(round(attr(p,'hpa_width')*100),'px'),height="650px")))
    })
    output[[paste0("hpa_plot_",k)]]<-shiny::renderPlot({
      shiny::req(available(),input[[paste0("hpa_metric_",k)]])
      core$report_hpa_selected_plot(available(),k,input[[paste0("hpa_metric_",k)]])
    },res=120)
    output[[paste0("hpa_table_",k)]]<-shiny::renderTable({
      x<-available();shiny::req(x)
      if(k=="summary")tab<-x$general else {
        kind<-input[[paste0("hpa_table_kind_",k)]];shiny::req(kind);tab<-as.data.frame(x[[k]][[kind]])
      }
      filter<-input[[paste0("hpa_filter_",k)]]
      if(!is.null(filter)&&nzchar(filter))tab<-tab[apply(tab,1,function(row)any(grepl(tolower(filter),tolower(row),fixed=TRUE))),,drop=FALSE]
      tab
    },striped=TRUE,spacing="s",rownames=FALSE)
    output[[paste0("hpa_csv_",k)]]<-shiny::downloadHandler(
      filename=function()paste0(available()$identity$canonical_symbol,"_HPA_",hpa_file_section(k),".csv"),
      content=function(file){x<-available();shiny::req(x);d<-core$hpa_section_table(x,k);shiny::req(nrow(d)>0);data.table::fwrite(d,file)})
    output[[paste0("hpa_export_",k)]]<-shiny::downloadHandler(
      filename=function()paste0(available()$identity$canonical_symbol,"_HPA_",hpa_file_section(k),".csv"),
      content=function(file){x<-available();shiny::req(x);d<-core$hpa_section_table(x,k);shiny::req(nrow(d)>0);data.table::fwrite(d,file)})
  })
  output$hpa_downloads<-shiny::renderUI({x<-available();shiny::req(x)
    shiny::tagList(shiny::h4("Human Protein Atlas"),shiny::div(class="d-flex flex-wrap gap-2",
      lapply(names(hpa_sections),function(k)if(nrow(core$hpa_section_table(x,k)))shiny::downloadButton(paste0("hpa_export_",k),paste(hpa_sections[[k]],"CSV")))),
      shiny::downloadButton("hpa_metadata","HPA metadata JSON"))
  })
  output$hpa_metadata<-shiny::downloadHandler(filename=function()paste0(available()$identity$canonical_symbol,"_HPA_metadata.json"),
    content=function(file){x<-available();shiny::req(x);jsonlite::write_json(list(metadata=x$metadata,identity=x$identity,warnings=x$warnings),file,auto_unbox=TRUE,pretty=TRUE,na="null")})
  invisible(available)
}
