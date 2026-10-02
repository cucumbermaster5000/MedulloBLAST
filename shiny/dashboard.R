explorer_loading_summary <- function(states) {
  list(pending=names(states)[is.na(states)|states%in%c('idle','loading','running')],
    unavailable=names(states)[!is.na(states)&!states%in%c('idle','loading','running','complete')])
}
explorer_toolbar <- function() {
  release<-read.dcf(file.path(root,"DESCRIPTION"),fields=c("Version","Date","Build","Build-Time"))
  has_build<-!is.na(release[1,"Build"])
  shiny::div(class="dashboard-toolbar",
  shiny::div(class="toolbar-query",shiny::textInput("gene","Gene symbol",value="",placeholder="e.g. HLX"),
    bslib::input_task_button("run","Run Analysis",label_busy="Analyzing...")),
  shiny::div(class="toolbar-context",
    shiny::span("MedulloBLAST. A medulloblastoma gene explorer."),
    shiny::tags$small("Designed by Dave Ng with the Pouponnot Lab 2026"),
    shiny::tags$small(id="app-release",paste0("Version ",release[1,"Version"]," | Latest update: ",format(as.Date(release[1,"Date"]),"%d %B %Y"),if(has_build)" (UTC)" else "")),
    if(has_build)shiny::tags$small(id="app-build",title=release[1,"Build-Time"],paste("Build",release[1,"Build"]))))
}
explorer_dashboard_style <- function() shiny::tags$style(shiny::HTML(paste(readLines(file.path(root,"shiny","dashboard.css"),warn=FALSE),collapse="\n")))
explorer_summary_cards <- function(x) {
  card<-function(title,value,note)shiny::div(class="summary-card",shiny::strong(title),shiny::div(class="stat-value",value),shiny::tags$small(note))
  fmt<-function(v)if(length(v)==1L&&is.finite(v))format(signif(v,3),trim=TRUE) else "Unavailable"
  s<-x$clinical$survival$all$statistics;m<-x$clinical$metastasis$all$statistics
  sub<-x$subtypes$summaries
  top<-if(!is.null(sub)&&nrow(sub))sub$subtype[which.max(sub$median)] else "Unavailable"
  shiny::div(class="summary-grid",
    card("Patient samples",nrow(x$retrieval$data),paste(x$analysis$statistics$n,"included in subgroup analysis")),
    card("Broad MB association",paste("p =",fmt(x$analysis$statistics$p_value)),"Omnibus raw p; see MB Subgroups"),
    card("Highest median subtype",top,"Descriptive ranking; see MB Subtypes"),
    card("Survival: complete MB",paste("HR",fmt(s$hazard_ratio_high_vs_low)),paste(s$cutoff_method,"cutoff; log-rank p =",fmt(s$logrank_p))),
    card("Metastasis: complete MB",paste("p =",fmt(m$p_value)),"Unadjusted association; see Metastasis"),
    card("Gene context",x$gene,"HPA, Functional Biology, DepMap and Publications below"))
}
