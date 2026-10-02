# Cavalli plot exports use completed session results only; no retrieval or fitting.
cavalli_export_analysis <- function(result,kind,cohort='all') {
  if(!kind %in% c('subgroups','subtypes','metastasis','survival'))stop('Unknown Cavalli plot')
  if(!cohort %in% c('all','wnt','shh','group3','group4'))stop('Unknown Cavalli cohort')
  switch(kind,subgroups=result$analysis,subtypes=result$subtypes,result$clinical[[kind]][[cohort]])
}
cavalli_export_payload <- function(result,kind,cohort='all',app_root='.') {
  provenance<-result$retrieval$provenance
  if(!identical(provenance$dataset_table,'ps_avgpres_gse85217geo763_hugene11t'))stop('This export supports the Cavalli dataset only')
  a<-cavalli_export_analysis(result,kind,cohort);p<-a$plot
  if(is.null(p)||!inherits(p,'ggplot')||!nrow(p$data))stop('No completed plot is available for export')
  audit<-switch(kind,subgroups=result$retrieval$data,a$data)
  plotted<-if(kind=='survival')a$data[which(a$data$included),,drop=FALSE]else p$data
  if(!nrow(plotted)||anyDuplicated(plotted$sample_id))stop('Invalid plotted patient identifiers')
  audit$included<-audit$sample_id %in% plotted$sample_id
  if(!'exclusion_reason' %in% names(audit))audit$exclusion_reason<-ifelse(audit$included,'',ifelse(!is.finite(audit$expression),'Non-finite expression','Missing subgroup'))
  columns<-c('sample_id','expression','expression_source_value','subgroup','subtype','plot_subtype','plot_subgroup','metastasis_raw','metastasis_label','time','event','expression_group','cutoff_method','cutoff_value')
  clean<-function(d,extra=character()) {
    d<-as.data.frame(d[,intersect(c(columns,extra),names(d)),drop=FALSE])
    for(k in names(d))if(is.factor(d[[k]]))d[[k]]<-as.character(d[[k]])
    d$gene<-result$gene;d$cohort<-cohort;rownames(d)<-NULL;d
  }
  tables<-list(plot_data=clean(plotted),patient_audit=clean(audit,c('included','exclusion_reason')))
  for(k in c('statistics','summaries','overall','pairwise','risk_table'))if(is.data.frame(a[[k]]))tables[[k]]<-a[[k]]
  if(kind %in% c('metastasis','survival')) {
    for(k in c('evidence','pairwise'))if(is.data.frame(result$clinical$selection[[k]]))tables[[paste0('cohort_selection_',k)]]<-result$clinical$selection[[k]]
  }
  if(kind=='survival')tables$curves<-as.data.frame(p$data[,c('time','survival','lower','upper','censored','group'),drop=FALSE])
  # Read the plot's trained scales, not a second implementation of its selection rules.
  built<-ggplot2::ggplot_build(p)
  if(kind=='survival') {
    scale<-built$plot$scales$get_scales('colour');groups<-as.character(scale$get_breaks())
    group_labels<-as.character(scale$get_labels(groups));colors<-as.character(scale$map(groups))
  } else {
    scale<-built$layout$panel_params[[1]]$x;groups<-as.character(scale$get_breaks())
    group_labels<-as.character(scale$get_labels());colors<-as.character(built$plot$scales$get_scales('fill')$map(groups))
  }
  labels<-lapply(p$labels[intersect(names(p$labels),c('title','subtitle','caption','x','y','colour','fill'))],as.character)
  release<-read.dcf(file.path(app_root,'DESCRIPTION'),fields=c('Version','Date','Build'))
  meta<-list(schema_version=1L,gene=result$gene,plot=kind,cohort=cohort,
    app=list(version=release[1,'Version'],release_date=release[1,'Date'],build=release[1,'Build']),
    source=c(list(attribution='Cavalli et al.; GSE85217. Expression and annotations retrieved through R2.',
      accession_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85217'),
      provenance[intersect(names(provenance),c('dataset_table','accession','dataset_label','queried_gene','reporter','transformation','retrieved_at','sample_join'))]),
    clinical=if(kind %in% c('metastasis','survival'))list(verified=a$verification$verified,
      units=a$verification$units,field=a$verification$field,time_field=a$verification$time_field,event_field=a$verification$event_field,
      coding_evidence=a$verification$evidence,cutoff=a$cutoff,
      selection=if(cohort=='all')'Complete MB cohort'else if(cohort %in% result$clinical[[paste0('manual_',kind)]])'Manually requested exploratory cohort'else 'Automatically selected exploratory cohort',
      cohort_selection_criteria=result$clinical$selection$criteria)else NULL,
    status=a$status %||% 'complete',warnings=unique(c(a$warnings,if(kind %in% c('metastasis','survival'))result$clinical$warnings)),
    plotting=list(labels=labels,groups=groups,group_labels=group_labels,colors=colors,jitter_seed=85217L,
      width=if(kind=='subtypes')13 else if(kind=='subgroups')8 else 9,height=if(kind=='subtypes')7.5 else if(kind=='subgroups')5.8 else 6),
    export=list(created_at=utc_now(),plotted_patients=nrow(plotted),audit_patients=nrow(audit),
      missing_value='NA',expression='Already on the source transformation recorded above; do not log-transform again.',
      survival='time is verified follow-up in years; event 1=death, 0=censored. High > cutoff; Low <= cutoff, including ties. The cutoff uses finite cohort expression before clinical exclusions.',
      reproduction='Uses exported plot inputs and stored curves; does not retrieve data or refit statistical models.'))
  list(tables=tables,metadata=meta)
}
cavalli_export_name <- function(x,format) {
  stopifnot(format %in% c('csv','zip'))
  m<-x$metadata;parts<-c(m$gene,'Cavalli',m$plot,m$cohort,if(m$plot=='survival')m$clinical$cutoff$method)
  paste0(paste(gsub('[^A-Za-z0-9._-]','_',parts),collapse='_'),'.',format)
}
cavalli_export_write_csv <- function(data,file) {
  data.table::fwrite(data,file,na='NA',quote='auto',row.names=FALSE)
  invisible(file)
}
write_cavalli_export <- function(x,file,format='csv',app_root='.') {
  if(format=='csv')return(cavalli_export_write_csv(x$tables$plot_data,file))
  if(format!='zip')stop('Unknown Cavalli export format')
  work<-tempfile('cavalli-export-');dir.create(work);on.exit(unlink(work,recursive=TRUE),add=TRUE)
  folder<-file.path(work,sub('[.]zip$','',cavalli_export_name(x,'zip')));dir.create(folder)
  for(k in names(x$tables))cavalli_export_write_csv(x$tables[[k]],file.path(folder,paste0(k,'.csv')))
  jsonlite::write_json(x$metadata,file.path(folder,'metadata.json'),auto_unbox=TRUE,pretty=TRUE,digits=NA,na='null',null='null')
  stopifnot(file.copy(file.path(app_root,'resources','cavalli-replot.R'),file.path(folder,'replot.R')))
  descriptions<-c(sample_id='Original R2 sample identifier (read as text).',gene='Queried gene.',cohort='Complete MB (all) or selected subgroup.',
    expression='Expression on the recorded R2 source scale; already transformed where specified.',expression_source_value='Original expression value string returned by R2.',
    subgroup='Original broad R2 subgroup label.',subtype='Original R2 molecular subtype label.',plot_subtype='Subtype label ordered by metadata.json plotting.groups.',plot_subgroup='Observed subtype-to-subgroup membership.',
    metastasis_raw='Verified source code: 0=M0; 1=metastatic.',metastasis_label='Verified metastasis label used on the plot.',
    time='Verified overall-survival follow-up in years.',event='Verified event: 1=death; 0=censored.',expression_group='High > cutoff; Low <= cutoff.',
    cutoff_method='Mean or median of finite cohort expression before clinical exclusions.',cutoff_value='Numeric threshold on the expression source scale.')
  columns<-names(x$tables$plot_data)
  notes<-c('MedulloBLAST - Cavalli plot reproduction package','',
    paste('Gene:',x$metadata$gene,'| Plot:',x$metadata$plot,'| Cohort:',x$metadata$cohort),
    x$metadata$source$attribution,'Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85217',
    'See metadata.json for the R2 dataset/reporter, retrieval date, transformation, release/build, options and coding evidence.','',
    'RECREATE THE PLOT',
    'Extract the ZIP. Install R and the ggplot2 and jsonlite packages if needed (the script does not install anything).',
    'Run: Rscript path/to/replot.R [output-directory]',
    'Or, in R: source("path/to/replot.R"); p <- cavalli_replot("path/to/extracted-package"); print(p)',
    'The standalone script writes reproduced_plot.png and reproduced_plot.pdf. It reads only this package, with no network or webapp access.',
    'Colours, group ordering, jitter settings and curve values match the app. Fonts and rendering may differ by platform/package version.',
    'Customize cavalli_replot() or add ggplot2 layers to p to create your own figures.','',
    'FILES',
    'plot_data.csv: only the patient records used for this plot; identical to the direct CSV download.',
    'patient_audit.csv: all candidate patients in this cohort, included flag and exclusion_reason. These records are not all plotted.',
    'Other CSVs: available analysis statistics, descriptive summaries and pairwise/cohort-selection evidence. Values are not rounded for presentation.',
    'Survival packages also contain curves.csv (time, survival, lower/upper pointwise 95% limits, censored count and group) and risk_table.csv (group, time_years, n_at_risk).',
    'Survival is drawn from the stored curves; no model is refitted. Origin rows at time zero and censoring marks are preserved.',
    'UTF-8, comma-separated; missing values are NA. Read sample IDs as text and preserve numeric precision.','',
    'PLOT DATA COLUMNS',paste0(columns,': ',unname(descriptions[columns])),'',
    'METHODS AND LIMITATIONS',x$metadata$export$expression,
    if(x$metadata$plot=='survival')x$metadata$export$survival,
    'Unadjusted exploratory comparisons do not establish causality or subgroup specificity. See statistics tables and metadata warnings for applicable test/adjustment details.',
    x$metadata$warnings)
  writeLines(enc2utf8(notes),file.path(folder,'README.txt'),useBytes=TRUE)
  archive<-file.path(work,'package.zip');report_zip_dir(folder,archive)
  if(!file.copy(archive,file,overwrite=TRUE))stop('Could not write Cavalli ZIP')
  invisible(file)
}
