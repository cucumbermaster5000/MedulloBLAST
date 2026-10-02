# Presentation and orchestration only. Scientific computations live in R/.
explorer_ui <- function() {
  biology_ui <- explorer_biology_ui()
  future <- function(title) bslib::nav_panel(title,
    bslib::card(bslib::card_header(title), shiny::p("Not implemented yet")))
  bslib::page_fluid(
    title = "MedulloBLAST",
    shiny::tags$head(shiny::tags$link(rel="icon",type="image/svg+xml",sizes="any",href="medulloblast-assets/cerebellum.svg?v=20261001"),
      shiny::tags$script(src="medulloblast-assets/loading.js"),
      shiny::tags$script(src="medulloblast-assets/session-policy.js")),
    theme = bslib::bs_theme(version = 5, primary = "#285DA8"),
    explorer_style(),explorer_dashboard_style(),explorer_toolbar(),explorer_capacity_notice(),shiny::uiOutput("reference_status"),shiny::uiOutput("status"),
    shiny::div(class="dashboard-nav",bslib::navset_pill_list(id = "section",widths=c(2,10),well=FALSE,
      bslib::nav_panel("Overview", shiny::uiOutput("gene_glance"),shiny::uiOutput("overview"),shiny::uiOutput("report_overview")),
      bslib::nav_panel("Pediatric Pan-Cancer",shiny::uiOutput("pfister_summary"),
        shiny::plotOutput("pfister_plot",height=explorer_plot_height("large")),explorer_report_control('pfister_pan_cancer'),
        shiny::h4("Expression summaries"),shiny::tableOutput("pfister_summaries"),
        shiny::h4("Sample types"),shiny::tableOutput("pfister_counts"),
        shiny::downloadButton("pfister_csv","Pfister patient CSV")),
      bslib::nav_panel("MB Subgroups",
        shiny::uiOutput("plot_heading"),
        shiny::div(style = "overflow-x:auto", shiny::div(
          shiny::plotOutput("subgroup_plot", height = explorer_plot_height()))),
        explorer_report_control('r2_subgroups'),shiny::h4("Statistical comparison"),
        shiny::div(style = "overflow-x:auto", shiny::tableOutput("statistics")),
        shiny::p(class = "small text-muted", "The original omnibus result is unchanged. Additional pairwise evidence for clinical cohort selection is shown below."),
        shiny::uiOutput("selection_summary"),
        shiny::div(style="overflow-x:auto",shiny::tableOutput("selection_pairs")),
        shiny::h4("Subgroup summaries"), shiny::tableOutput("summaries")),
      bslib::nav_panel("MB Subtypes",
        shiny::uiOutput("subtype_summary"),
        shiny::div(style="overflow-x:auto", shiny::div(
          shiny::plotOutput("subtype_plot", height=explorer_plot_height("large")))),
        explorer_report_control('r2_subtypes'),shiny::h4("Descriptive statistics"),
        shiny::div(style="overflow-x:auto", shiny::tableOutput("subtype_descriptives")),
        shiny::h4("Overall comparison"),
        shiny::div(style="overflow-x:auto", shiny::tableOutput("subtype_overall")),
        shiny::h4("All pairwise subtype comparisons"),
        shiny::p("Positive rank-biserial effect means subtype A tends higher than B; negative means lower. Median difference is A minus B on the source expression scale. BH adjustment covers all subtype pairs for this gene. These results do not establish subtype specificity."),
        shiny::div(style="overflow:auto;max-height:600px", shiny::tableOutput("subtype_pairs"))),
      bslib::nav_panel("Metastasis",
        shiny::checkboxGroupInput("metastasis_cohorts","Additional metastasis subgroups (exploratory)",
          choices=c("WNT"="wnt","SHH"="shh","Group 3"="group3","Group 4"="group4")),
        shiny::p("Complete MB and automatically selected subgroups remain available. Manual selection is exploratory, not evidence of expression enrichment."),
        shiny::uiOutput("metastasis_content")),
      bslib::nav_panel("Survival",
        shiny::checkboxGroupInput("survival_cohorts","Additional survival subgroups (exploratory)",
          choices=c("WNT"="wnt","SHH"="shh","Group 3"="group3","Group 4"="group4")),
        shiny::p("Complete MB and automatically selected subgroups remain available. Manual selection is exploratory, not evidence of expression enrichment."),
        shiny::selectInput("survival_cutoff","Survival cutoff method",choices=c("Mean"="mean","Median"="median"),selected="mean"),
        shiny::p("High > cutoff; Low <= cutoff, including ties. Each cohort uses its own finite expression values before clinical exclusions."),
        shiny::uiOutput("survival_content")),
      biology_ui$functional, explorer_hpa_ui(), explorer_cerebellum_ui(), biology_ui$targets, explorer_depmap_ui(),
      explorer_publications_ui(),
      bslib::nav_panel("Downloads", shiny::uiOutput("report_downloads"),shiny::uiOutput("downloads"),shiny::uiOutput("biology_downloads"),shiny::uiOutput("hpa_downloads"),shiny::uiOutput("publication_downloads"),shiny::uiOutput("depmap_downloads"))
    )),shiny::tags$details(class="explorer-methods",shiny::tags$summary("About / Methods and data sources"),
      shiny::p("Research exploration of gene expression, clinical associations, cell-line dependency and published evidence. Interpret statistical results with sample sizes and effect sizes."),
      shiny::p("Sources: R2 Cavalli / GSE85217 and Pfister informp3; Enrichr GO and pathway libraries; documented human TF catalogues, TRRUST and binding evidence; current public DepMap through Breadbox; PubMed / NCBI. Source releases and retrieval dates are shown with each module. External databases may change."),
      shiny::p(paste("Runtime:",R.version.string,"| Package versions are pinned in renv.lock. Exported figures retain their full publication dimensions.")))
  )
}

explorer_error <- function(error, stage) {
  detail <- conditionMessage(error)
  if (grepl("gene symbol|Gene absent|found 0", detail))
    return("Gene not found in the selected R2 dataset.")
  if (grepl("Specify reporter|exact reporter", detail))
    return("This gene has multiple R2 reporters. Reporter selection is not available in Version 1.")
  if (grepl("subgroup|annotated subgroups", detail, ignore.case = TRUE))
    return("MB subgroup annotation is unavailable for this dataset.")
  if (stage == "retrieval") return("R2 data retrieval failed. Please try again later.")
  "The subgroup comparison could not be completed. Please try another gene."
}

# Uses the stored data/plot; these handlers never run statistical analysis.
write_explorer_download <- function(result, kind, file) {
  switch(kind,
    subtype_patients = data.table::fwrite(result$subtypes$data, file, na="NA"),
    subtype_summaries = data.table::fwrite(result$subtypes$summaries, file, na="NA"),
    subtype_pairwise = data.table::fwrite(result$subtypes$pairwise, file, na="NA"),
    subtype_overall = data.table::fwrite(result$subtypes$overall, file, na="NA"),
    subtype_pdf = ggplot2::ggsave(file, result$subtypes$plot, device=grDevices::pdf, width=13, height=7.5),
    subtype_png = ggplot2::ggsave(file, result$subtypes$plot, device="png", width=13, height=7.5, dpi=300),
    patients = data.table::fwrite(result$retrieval$data, file, na = "NA"),
    statistics = data.table::fwrite(result$analysis$statistics, file, na = "NA"),
    pdf = ggplot2::ggsave(file, result$analysis$plot, device = grDevices::pdf,
      width = 8, height = 5.8),
    png = ggplot2::ggsave(file, result$analysis$plot, device = "png",
      width = 8, height = 5.8, dpi = 300),
    stop("Unknown download type"))
  invisible(file)
}

explorer_server <- function(core, cfg, reference_service=NULL) {
  force(core); force(cfg)
  function(input, output, session) {
    if(exists("explorer_session_policy",mode="function"))explorer_session_policy(input,output,session)
    if(exists("explorer_session_config",mode="function"))cfg<-explorer_session_config(cfg,session)
    if(exists("explorer_provisioning_server",mode="function"))explorer_provisioning_server(input,output,session,cfg,core,reference_service)
    result <- shiny::reactiveVal(NULL)
    failure <- shiny::reactiveVal(NULL)
    busy <- shiny::reactiveVal(FALSE)
    main_job<-NULL;main_started<-NULL;main_gene<-NULL;main_tick<-shiny::reactiveVal(0L)
    session$onSessionEnded(function(){if(!is.null(main_job) && main_job$is_alive())main_job$kill()})
    biology <- if(exists("explorer_biology_server",mode="function")) explorer_biology_server(input,output,session,cfg,result,core) else shiny::reactive(NULL)
    hpa <- if(exists("explorer_hpa_server",mode="function"))explorer_hpa_server(input,output,session,cfg,core,biology) else shiny::reactive(NULL)
    cerebellum <- if(exists("explorer_cerebellum_server",mode="function"))explorer_cerebellum_server(input,output,session,cfg,core) else shiny::reactive(NULL)
    publications <- if(exists("explorer_publications_server",mode="function"))explorer_publications_server(input,output,session,cfg,core) else shiny::reactive(NULL)
    depmap <- if(exists("explorer_depmap_server",mode="function"))explorer_depmap_server(input,output,session,cfg,core) else shiny::reactive(NULL)
    shiny::observeEvent(input$run, {
      if (busy() && (!isTRUE(cfg$public_app) || identical(main_gene,toupper(trimws(input$gene))))) return()
      if(!is.null(main_job) && main_job$is_alive())main_job$kill()
      busy(TRUE)
      if(!isTRUE(cfg$public_app))on.exit(busy(FALSE), add = TRUE)
      result(NULL); failure(NULL); pfister(NULL)
      gene <- toupper(trimws(input$gene))
      if (length(gene) != 1L || is.na(gene) || nchar(gene)>64 || !grepl("^[A-Z][A-Z0-9._-]*$", gene)) {
        failure("Please enter one valid gene symbol, for example HLX or SLC2A1.")
        busy(FALSE)
        return()
      }
      if(exists("explorer_log",mode="function"))explorer_log("query",gene)
      if(isTRUE(cfg$public_app)) {
        main_started<<-Sys.time();main_gene<<-gene
        main_job<<-tryCatch(start_explorer_worker(cfg,gene),error=function(e){failure("Analysis could not start. Please retry later.");busy(FALSE);NULL})
        main_tick(shiny::isolate(main_tick())+1L)
        return()
      }
      stage <- "retrieval"
      tryCatch(shiny::withProgress(message = paste("Analyzing", gene), value = 0, {
        shiny::incProgress(.15, detail = "Retrieving R2 patient expression")
        retrieval <- core$get_r2_expression(cfg, gene)
        stage <- "analysis"
        shiny::incProgress(.55, detail = "Comparing broad MB subgroups")
        if (!"subgroup" %in% names(retrieval$data) || all(core$r2_missing(retrieval$data$subgroup)))
          stop("MB subgroup annotation unavailable")
        analysis <- core$analyze_r2_subgroups(retrieval)
        shiny::incProgress(.15, detail="Comparing molecular subtypes")
        subtype_failure <- NULL
        subtypes <- tryCatch(core$analyze_r2_subtypes(retrieval), error=function(e) {
          message(sprintf("[Gene Explorer] gene=%s stage=subtypes: %s", gene, conditionMessage(e)))
          subtype_failure <<- "Molecular subtype analysis could not be completed. Broad subgroup results remain available."
          NULL
        })
        clinical <- NULL
        if(is.function(core$analyze_r2_clinical_v3)) clinical <- tryCatch(
          core$analyze_r2_clinical_v3(cfg,retrieval,analysis),error=function(e) {
            message("[Gene Explorer clinical] ",conditionMessage(e)); NULL
          })
        if(!is.null(clinical) && identical(input$survival_cutoff,"median"))
          clinical <- core$recalculate_r2_survival(retrieval,clinical,"median")
        if(!is.null(clinical) && length(input$survival_cohorts))
          clinical <- core$recalculate_r2_survival(retrieval,clinical,input$survival_cutoff,input$survival_cohorts)
        if(!is.null(clinical) && length(input$metastasis_cohorts))
          clinical <- core$recalculate_r2_metastasis(retrieval,clinical,input$metastasis_cohorts)
        result(list(gene = retrieval$provenance$queried_gene, retrieval = retrieval, analysis = analysis,
          subtypes=subtypes, subtype_failure=subtype_failure,clinical=clinical))
        shiny::incProgress(.3, detail = "Preparing figure and tables")
      }), error = function(e) {
        # Technical details go to the R console/server log, never the web page.
        message(sprintf("[Gene Explorer] gene=%s stage=%s: %s", gene, stage, conditionMessage(e)))
        failure(explorer_error(e, stage))
      })
    }, ignoreInit = TRUE)

    shiny::observe({
      main_tick();if(is.null(main_job))return()
      if(explorer_worker_expired(main_started,cfg,main_job)){
        main_job$kill();main_job<<-NULL;busy(FALSE);failure("R2 analysis timed out. Please retry later.");return()
      }
      if(main_job$is_alive()){shiny::invalidateLater(500,session);return()}
      x<-tryCatch(main_job$get_result(),error=function(e){explorer_log_error("R2 analysis failed",e);NULL});main_job<<-NULL;busy(FALSE)
      if(is.null(x)){failure("R2 data or analysis are currently unavailable. Please retry later.");return()}
      if(!is.null(x$clinical)) {
        cutoff<-if(identical(input$survival_cutoff,"median"))"median" else "mean"
        x$clinical<-core$recalculate_r2_survival(x$retrieval,x$clinical,cutoff,input$survival_cohorts)
        x$clinical<-core$recalculate_r2_metastasis(x$retrieval,x$clinical,input$metastasis_cohorts)
      }
      result(x);explorer_log("analysis complete",x$gene)
    })

    shiny::observeEvent(input$survival_cutoff, {
      x <- result()
      if(is.null(x$clinical) || !input$survival_cutoff %in% c("mean","median")) return()
      x$clinical <- core$recalculate_r2_survival(x$retrieval,x$clinical,input$survival_cutoff)
      result(x)
    },ignoreInit=TRUE)

    for(endpoint in c("survival","metastasis")) local({
      ep <- endpoint
      shiny::observeEvent(input[[paste0(ep,"_cohorts")]], {
        x <- result(); if(is.null(x$clinical)) return()
        chosen <- input[[paste0(ep,"_cohorts")]]
        if(is.null(chosen)) chosen <- character()
        x$clinical <- if(ep=="survival") core$recalculate_r2_survival(x$retrieval,x$clinical,input$survival_cutoff,chosen) else
          core$recalculate_r2_metastasis(x$retrieval,x$clinical,chosen)
        result(x)
      },ignoreInit=TRUE,ignoreNULL=FALSE)
    })

    # Load once for each completed gene analysis so the integrated report is
    # complete without requiring the user to visit this tab first.
    pfister <- shiny::reactiveVal(NULL)
    pfister_job<-NULL;pfister_started<-NULL;pfister_gene<-NULL;pfister_tick<-shiny::reactiveVal(0L)
    session$onSessionEnded(function(){if(!is.null(pfister_job) && pfister_job$is_alive())pfister_job$kill()})
    shiny::observeEvent(result()$gene, {
      shiny::req(result())
      gene <- result()$gene
      if(identical(pfister()$gene,gene)) return()
      if(isTRUE(cfg$public_app)) {
        if(!is.null(pfister_job) && pfister_job$is_alive() && identical(pfister_gene,gene))return()
        if(!is.null(pfister_job) && pfister_job$is_alive())pfister_job$kill()
        pfister_gene<<-gene;pfister_started<<-Sys.time()
        pfister_job<<-tryCatch(start_explorer_worker(cfg,gene,"pfister"),error=function(e){pfister(list(gene=gene,analysis=NULL));NULL})
        pfister_tick(shiny::isolate(pfister_tick())+1L);return()
      }
      a <- tryCatch(core$analyze_pfister_pan_cancer(core$get_pfister_expression(cfg,gene)),error=function(e) {
        message("[Pfister] ",conditionMessage(e));NULL
      })
      pfister(list(gene=gene,analysis=a))
    },ignoreInit=TRUE)
    shiny::observe({
      pfister_tick();if(is.null(pfister_job))return()
      if(explorer_worker_expired(pfister_started,cfg,pfister_job)){
        pfister_job$kill();pfister_job<<-NULL;pfister(list(gene=pfister_gene,analysis=NULL));return()
      }
      if(pfister_job$is_alive()){shiny::invalidateLater(500,session);return()}
      a<-tryCatch(pfister_job$get_result(),error=function(e){explorer_log_error("Pfister retrieval failed",e);NULL});pfister_job<<-NULL
      if(identical(result()$gene,pfister_gene))pfister(list(gene=pfister_gene,analysis=a))
    })
    output$pfister_summary <- shiny::renderUI({
      shiny::req(result())
      if(!identical(pfister()$gene,result()$gene)) return(shiny::p("Loading Pfister expression for the current gene..."))
      a <- pfister()$analysis
      if(is.null(a)) return(shiny::p("Pfister expression could not be retrieved for this gene. MB results remain available. Run the gene again to retry."))
      shiny::tagList(shiny::h3(paste(result()$gene,"Pfister pan-cancer expression")),shiny::p(a$provenance$dataset_label),
        shiny::p(sprintf("%d samples; %d cancer types. All sample types retained; descriptive comparison only.",nrow(a$data),nrow(a$summaries))),
        shiny::p("Expression is log2(1 + FPKM). Do not directly compare its numerical scale to the Cavalli microarray dataset."))
    })
    pfister_current <- function() {shiny::req(result(),identical(pfister()$gene,result()$gene),pfister()$analysis);pfister()$analysis}
    output$pfister_plot <- shiny::renderPlot({explorer_web_plot(pfister_current()$plot,"pfister_plot",session)},res=120,execOnResize=TRUE)
    output$pfister_summaries <- shiny::renderTable({pfister_current()$summaries},digits=3)
    output$pfister_counts <- shiny::renderTable({pfister_current()$counts})
    output$pfister_csv <- shiny::downloadHandler(filename=function()paste0(result()$gene,"_Pfister_patient_data.csv"),
      content=function(file)data.table::fwrite(pfister_current()$data,file,na="NA"))

    if(exists("explorer_report_server",mode="function"))
      explorer_report_server(input,output,session,cfg,core,result,pfister,biology,hpa,depmap,publications,cerebellum=cerebellum)

    loading_summary<-shiny::reactive({
      states<-c(R2=if(busy())'loading'else if(!is.null(failure()))'failed'else if(is.null(result()))'idle'else 'complete',
        `Pediatric pan-cancer`=if(is.null(result())){if(!is.null(failure()))'unavailable'else 'idle'}else if(!identical(pfister()$gene,result()$gene))'loading'else if(is.null(pfister()$analysis))'unavailable'else 'complete')
      modules<-c(biology='Functional biology',hpa='Protein Atlas',cerebellum='Developing cerebellum',publications='Publications',depmap='DepMap')
      for(k in names(modules)) {
        state_fn<-session$userData[[paste0(k,'_loading_state')]]
        if(is.function(state_fn))states[modules[[k]]]<-state_fn()
      }
      explorer_loading_summary(states)
    })
    output$status <- shiny::renderUI({
      if(is.null(input$run)||input$run==0)return(shiny::div(class='alert alert-info','Enter a gene and click Run Analysis to begin.'))
      if(!grepl('^[A-Za-z][A-Za-z0-9._-]{0,63}$',trimws(input$gene)))return(shiny::div(class='alert alert-danger',role='alert','Please enter one valid gene symbol.'))
      s<-loading_summary()
      shiny::tagList(
        if(!is.null(failure()))shiny::div(class='alert alert-warning',role='alert',failure()),
        if(length(s$pending))shiny::div(class='alert alert-info',role='status',
          shiny::icon('spinner',class='fa-spin'),shiny::strong(' Loading analysis and data...'),
          shiny::div(paste('Still loading:',paste(s$pending,collapse=', '))))
        else shiny::tagList(
          shiny::div(class='alert alert-info analysis-rendering',role='status',shiny::icon('spinner',class='fa-spin'),' Preparing visualizations...'),
          shiny::div(class=paste('analysis-ready alert',if(length(s$unavailable))'alert-warning'else 'alert-success'),role='status',
            if(length(s$unavailable))paste('Loading finished. Available results are ready to visualize. Unavailable:',paste(s$unavailable,collapse=', '))
            else paste('Data ready to visualize:',result()$gene))))
    })
    output$overview <- shiny::renderUI({
      x <- result()
      if (is.null(x)) return(shiny::tagList(shiny::p("No expression analysis yet."),shiny::h4("DepMap"),shiny::uiOutput("depmap_overview"),shiny::h4("Publications"),shiny::uiOutput("publications_overview")))
      s <- x$analysis$statistics
      metadata <- names(x$retrieval$metadata)
      shiny::tagList(
        if(exists("explorer_summary_cards",mode="function"))explorer_summary_cards(x) else
          bslib::layout_columns(bslib::value_box("Queried gene",x$gene),bslib::value_box("Patient samples retrieved",nrow(x$retrieval$data)),bslib::value_box("Samples analyzed",s$n)),
        shiny::h4("Dataset"), shiny::p(x$retrieval$provenance$dataset_label),
        shiny::p("GSE85217 - R2 log2 expression - WNT, SHH, Group 3 and Group 4"),
        shiny::h4("Broad subgroup result"),
        shiny::p(sprintf("Kruskal-Wallis comparison across %d subgroups: H(%d) = %.3f, raw p = %s. %d samples excluded for missing expression or subgroup annotation.",
          s$groups, s$df, s$statistic, format(s$p_value, digits = 4, scientific = TRUE), s$excluded_n)),
        shiny::div(class = "alert alert-warning", "This is an unadjusted comparison of expression distributions. It does not establish subgroup specificity or therapeutic benefit; age, batch and other confounders are not adjusted."),
        shiny::h4("Molecular subtype result"),
        shiny::p(if(is.null(x$subtypes)) "Unavailable" else sprintf("%d subtypes; overall p = %s",nrow(x$subtypes$summaries),format(x$subtypes$overall$p_value,digits=4))),
        shiny::h4("Pediatric pan-cancer"),
        shiny::p(if(is.null(pfister()$analysis))"Open Pediatric Pan-Cancer to load the independent Pfister cohort." else sprintf("%d samples across %d cancer types; descriptive comparison on log2(1 + FPKM).",nrow(pfister()$analysis$data),nrow(pfister()$analysis$summaries))),
        shiny::h4("Clinical results"),
        if(is.null(x$clinical)) shiny::p("Clinical analyses unavailable; expression results remain available.") else shiny::tagList(
          shiny::p(paste("Selected subgroups:",if(length(x$clinical$selection$selected)) paste(x$clinical$selection$selected,collapse=", ") else "None; complete MB cohort only")),
          shiny::p(x$clinical$selection$criteria),
          lapply(c("metastasis","survival"),function(endpoint) shiny::tagList(
            shiny::h5(endpoint),lapply(x$clinical[[endpoint]],function(a) {
              s <- a$statistics
              if(a$status!="complete") return(shiny::p(paste(a$cohort,paste(a$warnings,collapse=" "),sep=": ")))
              if(endpoint=="metastasis") shiny::p(sprintf("%s: metastatic N=%d; M0 N=%d; raw p=%.4g; rank-biserial effect=%.3f.",a$cohort,s$n_metastatic,s$n_m0,s$p_value,s$rank_biserial_metastatic_vs_m0)) else
                shiny::tagList(shiny::p(sprintf("%s: survival N=%d; events=%d; %s cutoff=%.10g; log-rank p=%.4g; HR High/Low=%.3g (95%% CI %.3g-%.3g).",a$cohort,s$n,s$events,s$cutoff_method,s$cutoff,s$logrank_p,s$hazard_ratio_high_vs_low,s$ci_lower,s$ci_upper)),
                  lapply(a$warnings[-1],shiny::p))
            }))),shiny::p(class="alert alert-warning",x$clinical$warnings)),
        shiny::h4("Functional Biology and TF evidence"),
        if(is.null(biology())) shiny::p("Functional evidence is loading or unavailable; see Functional Biology for status.") else shiny::tagList(
          shiny::p(biology()$annotation$summary),shiny::p(paste("Transcription factor:",biology()$tf$status,"| Family/domain:",biology()$tf$family)),
          shiny::p(paste("Unique target identifiers from sources (before enrichment symbol validation):",nrow(biology()$targets$summary))),
          if(!is.null(biology()$enrichment[[input$biology_target_set]])) {
            e<-biology()$enrichment[[input$biology_target_set]]
            if(nrow(e$rows))shiny::p(paste("Top target-associated terms (selected target set; adjusted p <= 0.05):",
              paste(head(e$rows$term[e$rows$adjusted_p_value<=.05],3),collapse="; "))) else shiny::p(paste(e$warnings,collapse=" "))
          }),
        shiny::h4("Human Protein Atlas"),shiny::uiOutput("hpa_overview"),
        shiny::h4("DepMap"),shiny::uiOutput("depmap_overview"),
        shiny::h4("Publications"),shiny::uiOutput("publications_overview"),
        shiny::h4("Available metadata fields"), shiny::p(paste(metadata, collapse = ", ")),
        shiny::p(class = "small text-muted", paste("Retrieved:", x$retrieval$provenance$retrieved_at,
          if (isTRUE(x$retrieval$provenance$cache_hit)) "- Existing R2 cache reused" else "- Retrieved from R2"))
      )
    })
    output$plot_heading <- shiny::renderUI({
      shiny::req(result()); shiny::h3(paste(result()$gene, "expression across MB subgroups"))
    })
    output$subgroup_plot <- shiny::renderPlot({
      shiny::req(result()); explorer_web_plot(result()$analysis$plot,"subgroup_plot",session)
    }, res=120,execOnResize=TRUE)
    output$statistics <- shiny::renderTable({
      shiny::req(result()); s <- result()$analysis$statistics
      data.frame(Comparison = s$test, N = s$n, `Effect size` = signif(s$epsilon_squared, 4),
        Statistic = signif(s$statistic, 5), `Raw p-value` = format(s$p_value, scientific = TRUE, digits = 4),
        `Adjusted p-value` = "Not applicable (one test)", check.names = FALSE)
    }, striped = TRUE, bordered = FALSE, spacing = "m", digits = 4)
    output$summaries <- shiny::renderTable({
      shiny::req(result()); result()$analysis$summaries
    }, striped = TRUE, digits = 3)
    output$subtype_summary <- shiny::renderUI({
      x <- result()
      if (is.null(x)) return(shiny::p("Run an analysis to view molecular subtypes."))
      if (is.null(x$subtypes)) return(shiny::div(class="alert alert-warning", x$subtype_failure))
      a <- x$subtypes
      shiny::tagList(shiny::h3(paste(x$gene,"molecular subtype expression")),
        shiny::p(a$metadata$classification),
        shiny::p(sprintf("Metadata field: %s | Annotated: %d | Missing subtype: %d | Subtypes: %d | Analyzed: %d | Excluded: %d",
          a$metadata$field,a$metadata$annotated_n,a$metadata$missing_n,nrow(a$summaries),a$overall$n,a$missing$excluded_n)),
        shiny::p(a$diagnostics$rationale),
        lapply(a$warnings, function(w) shiny::div(class="alert alert-warning",w)))
    })
    output$subtype_plot <- shiny::renderPlot({
      shiny::req(result()$subtypes); explorer_web_plot(result()$subtypes$plot,"subtype_plot",session)
    }, res=120,execOnResize=TRUE)
    output$subtype_descriptives <- shiny::renderTable({
      shiny::req(result()$subtypes); result()$subtypes$summaries
    }, striped=TRUE, digits=3)
    # Format only at presentation time; CSVs retain full numeric precision.
    subtype_table <- function(x) {
      for (field in intersect(c("p_value","adjusted_p_value"), names(x)))
        x[[field]] <- ifelse(is.na(x[[field]]), "Not applicable", format(x[[field]],scientific=TRUE,digits=4))
      x
    }
    output$subtype_overall <- shiny::renderTable({
      shiny::req(result()$subtypes); subtype_table(result()$subtypes$overall)
    }, striped=TRUE, digits=4)
    output$subtype_pairs <- shiny::renderTable({
      shiny::req(result()$subtypes); subtype_table(result()$subtypes$pairwise)
    }, striped=TRUE, digits=4)
    output$downloads <- shiny::renderUI({
      if (is.null(result())) return(shiny::p("Run an analysis to enable downloads."))
      shiny::tagList(shiny::h3(paste("Download", result()$gene, "results")),
        shiny::p("Exports contain the same patient data, statistics and figure shown in this analysis."),
        shiny::div(class = "d-flex flex-wrap gap-3",
          shiny::downloadButton("patients", "Patient expression CSV"),
          shiny::downloadButton("stats_csv", "Subgroup statistics CSV"),
          shiny::downloadButton("figure_pdf", "Figure PDF"),
          shiny::downloadButton("figure_png", "Figure PNG")),
        if (!is.null(result()$subtypes)) shiny::tagList(shiny::h4("Molecular subtypes"),
          shiny::div(class="d-flex flex-wrap gap-3",
            shiny::downloadButton("subtype_patients", "Subtype patient CSV"),
            shiny::downloadButton("subtype_summaries", "Subtype descriptive CSV"),
            shiny::downloadButton("subtype_pairwise", "Subtype pairwise CSV"),
            shiny::downloadButton("subtype_overall_csv", "Subtype overall test CSV"),
            shiny::downloadButton("subtype_pdf", "Subtype PDF"),
            shiny::downloadButton("subtype_png", "Subtype PNG (300 dpi)"))),
        shiny::uiOutput("clinical_downloads"))
    })
    downloads <- list(patients = c("patients", "csv", "text/csv"),
      stats_csv = c("statistics", "csv", "text/csv"),
      figure_pdf = c("pdf", "pdf", "application/pdf"),
      figure_png = c("png", "png", "image/png"),
      subtype_patients=c("subtype_patients","csv","text/csv"),
      subtype_summaries=c("subtype_summaries","csv","text/csv"),
      subtype_pairwise=c("subtype_pairwise","csv","text/csv"),
      subtype_overall_csv=c("subtype_overall","csv","text/csv"),
      subtype_pdf=c("subtype_pdf","pdf","application/pdf"),
      subtype_png=c("subtype_png","png","image/png"))
    for (id in names(downloads)) local({
      spec <- downloads[[id]]
      output[[id]] <- shiny::downloadHandler(
        filename = function() { shiny::req(result()); paste0(result()$gene, "_", spec[1], ".", spec[2]) },
        contentType = spec[3],
        content = function(file) { shiny::req(result()); write_explorer_download(result(), spec[1], file) })
    })
    output$selection_summary <- shiny::renderUI({
      shiny::req(result()$clinical); s <- result()$clinical$selection
      shiny::tagList(shiny::h4("Clinical cohort selection"),shiny::p(s$criteria),
        shiny::p(paste("Selected:",if(length(s$selected)) paste(s$selected,collapse=", ") else "None")))
    })
    output$selection_pairs <- shiny::renderTable({
      shiny::req(result()$clinical); result()$clinical$selection$pairwise
    },digits=6,striped=TRUE)
    for(endpoint in c("metastasis","survival")) local({
      ep <- endpoint
      output[[paste0(ep,"_content")]] <- shiny::renderUI({
        x <- result()$clinical
        if(is.null(x)) return(shiny::p("Clinical results are not available. Run an analysis or retry retrieval."))
        v <- x$verification[[ep]]
        shiny::tagList(shiny::h3(paste(result()$gene,ep)),shiny::h4("Metadata verification"),
          shiny::p(v$evidence),
          if(ep=="metastasis") shiny::p(sprintf("Field: %s | Type: %s | Missing: %d | Annotated: %d",v$field,v$source_type,v$missing_n,v$usable_n)) else
            shiny::p(sprintf("Time: %s | Event: %s | Units: %s | Missing event: %d | Missing time: %d | Verified usable: %d",v$time_field,v$event_field,v$units,v$missing_event_n,v$missing_time_n,v$usable_n)),
          shiny::tableOutput(paste0(ep,"_codes")),
          shiny::p(class="alert alert-warning",x$warnings),
          lapply(names(x[[ep]]),function(cohort) {
            a <- x[[ep]][[cohort]]; id <- paste(ep,cohort,sep="_")
            shiny::tagList(shiny::h4(if(cohort=="all") "Complete MB cohort" else paste("Subgroup:",cohort)),
              if(cohort %in% x[[paste0("manual_",ep)]]) shiny::p("Manually requested exploratory cohort; selection is not evidence of clinical significance."),
              lapply(a$warnings,function(w)shiny::p(w)),
              shiny::div(style="overflow-x:auto",shiny::tableOutput(paste0(id,"_counts"))),
              shiny::div(style="overflow-x:auto",shiny::tableOutput(paste0(id,"_stats"))),
              if(!is.null(a$plot)) shiny::tagList(shiny::plotOutput(paste0(id,"_plot"),height=explorer_plot_height()),explorer_report_control(id,input)),
              shiny::div(style="overflow-x:auto",shiny::tableOutput(paste0(id,"_details"))))
          }))
      })
      output[[paste0(ep,"_codes")]] <- shiny::renderTable({
        shiny::req(result()$clinical)
        v <- result()$clinical$verification[[ep]]
        if(ep=="metastasis") v$values else v$event_values
      })
      for(cohort in c("all","wnt","shh","group3","group4")) local({
        co <- cohort; id <- paste(ep,co,sep="_")
        get_analysis <- function() {shiny::req(result()$clinical[[ep]][[co]]);result()$clinical[[ep]][[co]]}
        output[[paste0(id,"_plot")]] <- shiny::renderPlot({explorer_web_plot(get_analysis()$plot,paste0(id,"_plot"),session)},res=120,execOnResize=TRUE)
        output[[paste0(id,"_counts")]] <- shiny::renderTable({get_analysis()$missing})
        output[[paste0(id,"_stats")]] <- shiny::renderTable({
          d <- get_analysis()$statistics
          d <- d[,setdiff(names(d),c("warnings","coding_evidence")),drop=FALSE]
          for(f in names(d)[vapply(d,is.numeric,logical(1))]) d[[f]] <- format(d[[f]],digits=16,trim=TRUE)
          d
        },striped=TRUE)
        output[[paste0(id,"_details")]] <- shiny::renderTable({
          a <- get_analysis(); if(ep=="metastasis") a$summaries else a$risk_table
        },digits=4,striped=TRUE)
        for(kind in c("patient_data","statistics","plot_pdf","plot_png")) local({
          k <- kind; ext <- if(k=="plot_pdf") "pdf" else if(k=="plot_png") "png" else "csv"
          output[[paste0(id,"_download_",k)]] <- shiny::downloadHandler(
            filename=function()paste0(result()$gene,if(co=="all") "" else paste0("_",co),"_",ep,"_",sub("_pdf$|_png$","",k),".",ext),
            content=function(file)core$write_r2_clinical_download(get_analysis(),k,file))
        })
      })
    })
    output$clinical_downloads <- shiny::renderUI({
      shiny::req(result()$clinical)
      shiny::tagList(shiny::h4("Clinical analyses"),lapply(c("metastasis","survival"),function(ep)
        lapply(names(result()$clinical[[ep]]),function(co) {
          a <- result()$clinical[[ep]][[co]]
          if(is.null(a$data)) return(NULL)
          kinds <- c("patient_data","statistics",if(!is.null(a$plot)) c("plot_pdf","plot_png"))
          shiny::tagList(shiny::h5(paste(ep,co)),shiny::div(class="d-flex flex-wrap gap-3",
            lapply(kinds,function(k)shiny::downloadButton(paste(ep,co,"download",k,sep="_"),gsub("_"," ",k)))))
        })))
    })
  }
}
