explorer_cerebellum_ui<-function() {
  panel<-function(species,kind) bslib::card(class='dashboard-card',
    bslib::card_header(if(kind=='expression')shiny::tagList('Gene expression: ',shiny::textOutput(paste0('cerebellum_title_',species),inline=TRUE)) else paste(tools::toTitleCase(species),'cell annotation')),
    shiny::uiOutput(paste0('cerebellum_',kind,'_state_',species)),
    plotly::plotlyOutput(paste0('cerebellum_',kind,'_',species),height='380px'),
    explorer_report_control(paste('cerebellum',kind,species,sep='_')))
  bslib::nav_panel('snRNA-seq Developing Cerebellum',
    shiny::h3('Developing Cerebellum \u2014 single-nucleus RNA-seq'),
    shiny::p('Explore expression of the queried gene across human and mouse cerebellar development.'),
    bslib::card(class='dashboard-card',shiny::textOutput('cerebellum_query'),
      shiny::selectInput('cerebellum_annotation','Colour UMAP by:',
        c('Broad lineage'='broad_lineage','Cell type'='cell_type','Cell subtype'='subtype','Cell state'='dev_state'),
        selected='broad_lineage',selectize=FALSE),
      shiny::p(class='small text-muted','Expression: exonic UMI counts | Species: Human + Mouse'),
      shiny::uiOutput('cerebellum_status')),
    shiny::p(class='small text-muted','Colour shows raw exonic UMI on a log(1 + count) scale, with a separate range for each species. Raw counts are not directly comparable across species. Grey: zero; pale grey: no exonic measurement.'),
    shiny::h4('Human developing cerebellum'),
    shiny::div(class='cerebellum-grid',panel('human','annotation'),panel('human','expression')),
    shiny::h4('Mouse developing cerebellum'),
    shiny::div(class='cerebellum-grid',panel('mouse','annotation'),panel('mouse','expression')),
    shiny::p(class='small text-muted','For responsiveness, overview and background points use a fixed sample stratified by cell type and subtype. Every expressing nucleus is retained, at its published coordinates. Plot downloads contain all nuclei.'),
    shiny::uiOutput('cerebellum_downloads'),
    shiny::div(class='small text-muted',
      shiny::p('Data source: Sepp M., Leiss K. et al., Cellular development and evolution of the mammalian cerebellum, Nature (2023). Kaessmann Lab / Pfister Lab. Integration adapted for Human and Mouse, exonic counts, and local indexed storage. Source app: CC BY 4.0.'),
      shiny::a('Original UMAP explorer',href='https://apps.kaessmannlab.org/sc-cerebellum-transcriptome/',target='_blank',rel='noopener'), ' | ',
      shiny::a('Publication',href='https://doi.org/10.1038/s41586-023-06884-x',target='_blank',rel='noopener'), ' | ',
      shiny::a('Source code',href='https://gitlab.com/kaessmannlab/shiny-mammalian-cerebellum',target='_blank',rel='noopener'), ' | ',
      shiny::a('CC BY 4.0',href='https://creativecommons.org/licenses/by/4.0/',target='_blank',rel='noopener')))
}

start_cerebellum_worker<-function(cfg,gene) callr::r_bg(function(cfg,gene){
  source(file.path(cfg$app_root,'R/load_core.R'));core<-load_r2_core(cfg$app_root)
  lock<-core$deployment_worker_slot(cfg);on.exit(filelock::unlock(lock),add=TRUE)
  core$get_cerebellum_profile(cfg,gene)
},args=list(cfg=cfg,gene=gene),libpath=.libPaths(),supervise=TRUE)

explorer_cerebellum_server<-function(input,output,session,cfg,core,worker_start=start_cerebellum_worker) {
  state<-shiny::reactiveVal(NULL);query<-shiny::reactiveVal(NULL);busy<-shiny::reactiveVal(FALSE)
  tick<-shiny::reactiveVal(0L);job<-NULL;started<-NULL
  pending<-shiny::reactiveVal(NULL)
  reference<-shiny::reactive({f<-session$userData$reference_status;x<-if(is.function(f))f()[['cerebellum']]else NULL;if(is.null(x))list(state='ready')else x})
  session$userData$cerebellum_loading_state<-shiny::reactive(if(!is.null(query())&&reference()$state!='ready'){if(reference()$state=='failed')'failed'else'loading'}else if(busy())'loading'else if(is.null(state()))'idle'else if(isTRUE(state()$error))'failed'else 'complete')
  assets<-shiny::reactive({if(reference()$state!='ready')return(NULL);tryCatch(core$cerebellum_assets(cfg),error=function(e){
    message('[Cerebellum assets] ',conditionMessage(e));NULL})})
  launch<-function(gene){
    busy(TRUE);started<<-Sys.time()
    job<<-tryCatch(worker_start(cfg,gene),error=function(e){message('[Cerebellum worker] ',conditionMessage(e));NULL})
    if(is.null(job)){busy(FALSE);state(list(error=TRUE))}
    tick(shiny::isolate(tick())+1L)
  }
  shiny::observeEvent(input$run,{
    if(!is.null(job)&&job$is_alive())job$kill()
    job<<-NULL;state(NULL);busy(FALSE);pending(NULL)
    gene<-tryCatch(core$biology_symbol(input$gene),error=function(e)NULL);query(gene)
    if(is.null(gene))return()
    if(reference()$state!='ready'){pending(gene);return()}
    launch(gene)
  },ignoreInit=TRUE)
  shiny::observe({gene<-pending();if(!is.null(gene)&&reference()$state=='ready'){pending(NULL);launch(gene)}})
  shiny::observe({
    tick();if(is.null(job))return()
    if(exists('explorer_worker_expired',mode='function')&&explorer_worker_expired(started,cfg)){job$kill();job<<-NULL;busy(FALSE);state(list(error=TRUE));return()}
    if(job$is_alive()){shiny::invalidateLater(250,session);return()}
    x<-tryCatch(job$get_result(),error=function(e){message('[Cerebellum query] ',conditionMessage(e));list(error=TRUE)})
    job<<-NULL;busy(FALSE)
    if(isTRUE(x$error)||identical(x$query,shiny::isolate(query())))state(x)
  })
  session$onSessionEnded(function(){if(!is.null(job)&&job$is_alive())job$kill()})
  output$cerebellum_query<-shiny::renderText(paste('Queried gene:',if(is.null(query()))'Run an analysis' else query()))
  output$cerebellum_expression_heading<-shiny::renderText(paste('Expression of',if(is.null(query()))'the queried gene' else query()))
  output$cerebellum_status<-shiny::renderUI({
    if(reference()$state!='ready')return(shiny::p(class='alert alert-info',if(reference()$state=='failed')'Reference download failed. Use Retry dataset downloads above; this is not a missing-gene result.'else'Developing cerebellum data are being prepared. Your latest gene query will run when the files are ready.'))
    if(busy())return(shiny::p(class='alert alert-info',shiny::icon('spinner',class='fa-spin'),' Loading developing cerebellum expression...'))
    if(is.null(query()))return(shiny::p('Run a gene analysis to explore this dataset.'))
    if(is.null(assets())||isTRUE(state()$error))return(shiny::p(class='alert alert-info','Developing cerebellum data are currently unavailable. Other analyses remain available.'))
    x<-state();shiny::req(x)
    shiny::tagList(shiny::p(paste('Human:',x$identity$human$symbol %||% 'unresolved','| Mouse:',x$identity$mouse$symbol %||% 'no 1:1 orthologue')),
      shiny::p(class='small text-muted','Gene mapping: saved Ensembl human-mouse one-to-one orthology. No paralogue substitution.'))
  })
  `%||%`<-core$`%||%`
  for(s in c('human','mouse'))local({species<-s;label<-tools::toTitleCase(s)
    output[[paste0('cerebellum_title_',s)]]<-shiny::renderText({x<-state();paste(label,'-',x$identity[[species]]$symbol %||% query() %||% 'gene expression')})
    output[[paste0('cerebellum_annotation_state_',s)]]<-shiny::renderUI({
      if(reference()$state!='ready')return(shiny::p('Reference data are not ready; see dataset preparation status above.'))
      a<-assets();if(is.null(a))return(shiny::p('Published embedding data are currently unavailable.'))
      x<-a$species[[species]];shiny::p(class='small text-muted',paste(format(nrow(x$cells),big.mark=','),'published nuclei;',format(length(x$display_cells),big.mark=','),'shown.'))
    })
    output[[paste0('cerebellum_expression_state_',s)]]<-shiny::renderUI({
      if(reference()$state!='ready')return(shiny::p('Waiting for reference data, not a gene-expression result.'))
      if(busy())return(shiny::p('Loading developing cerebellum expression...'))
      if(is.null(query()))return(shiny::p('Run a gene analysis to view expression.'))
      x<-state()$species[[species]]
      if(is.null(x))return(shiny::p('Expression data are currently unavailable.'))
      if(x$status!='complete')return(shiny::p(x$status))
      shiny::p(class='small text-muted',paste(format(x$n_expressing,big.mark=','),'expressing /',format(x$n_measured,big.mark=','),'measured nuclei.',
        if(!x$n_expressing)'No exonic expression detected.' else ''))
    })
    output[[paste0('cerebellum_annotation_',s)]]<-plotly::renderPlotly({
      a<-assets();shiny::req(a);field<-input$cerebellum_annotation %||% 'broad_lineage'
      shiny::req(field %in% unname(core$CEREBELLUM_ANNOTATIONS))
      core$cerebellum_plot(a$species[[species]],field=field,colors=core$cerebellum_colors(a,field))
    })
    output[[paste0('cerebellum_expression_',s)]]<-plotly::renderPlotly({
      a<-assets();x<-state()$species[[species]];shiny::req(a,!busy(),identical(x$status,'complete'))
      core$cerebellum_plot(a$species[[species]],expression=x)
    })
    for(k in c('expression','annotation'))for(f in c('png','pdf'))local({kind<-k;format<-f
      output[[paste0('cerebellum_download_',species,'_',kind,'_',format)]]<-shiny::downloadHandler(
        filename=function()paste0(query() %||% 'cerebellum','_',species,'_',kind,'_UMAP.',format),
        content=function(file){a<-assets();shiny::req(a);x<-if(kind=='expression')state()$species[[species]] else NULL
          if(kind=='expression')shiny::req(identical(x$status,'complete'))
          field<-input$cerebellum_annotation %||% 'broad_lineage'
          p<-core$cerebellum_export_plot(a$species[[species]],x,paste(label,if(kind=='expression')x$gene$symbol else 'cell annotation'),field,core$cerebellum_colors(a,field))
          ggplot2::ggsave(file,p,device=format,width=9,height=7,dpi=300)
        })
    })
  })
  output$cerebellum_downloads<-shiny::renderUI({shiny::req(assets())
    buttons<-list()
    for(s in c('human','mouse'))for(k in c('expression','annotation')) {
      if(k=='expression'&&!identical(state()$species[[s]]$status,'complete'))next
      for(f in c('png','pdf'))buttons[[length(buttons)+1L]]<-shiny::downloadButton(paste0('cerebellum_download_',s,'_',k,'_',f),paste(tools::toTitleCase(s),k,toupper(f)))
    }
    shiny::div(class='d-flex flex-wrap gap-2',buttons)
  })
  invisible(shiny::reactive({
    x<-state();if(is.null(x))return(NULL)
    a<-assets();if(!is.null(a)) {
      x$cells<-lapply(a$species,function(z)z$cells)
      x$display_cells<-lapply(a$species,function(z)z$display_cells)
    }
    x
  }))
}
