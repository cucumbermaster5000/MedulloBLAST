# Deterministic, versioned aggregation of existing session results. No I/O or tests.
REPORT_EVIDENCE_VERSION <- '1.0.0'
report_source_details <- function(citation) {
  m<-citation$metadata
  if(!is.list(m))return('Full provenance is retained in the report metadata package.')
  keys<-intersect(c('dataset_label','dataset_id','release','version','retrieved_at','quantification','source_commit'),names(m))
  values<-vapply(keys,function(k) {
    v<-m[[k]]
    if(!is.atomic(v)||!length(v)||all(is.na(v)))return('')
    paste(gsub('_',' ',k),paste(v,collapse=', '),sep=': ')
  },character(1))
  paste(c(values[nzchar(values)],'Full provenance is retained in the report metadata package.'),collapse='; ')
}
report_number <- function(x) length(x)==1L && is.numeric(x) && is.finite(x)
report_citations <- function(s) {
  entry <- function(id,label,url,metadata) list(id=id,label=label,url=url,metadata=metadata)
  out <- list(R2=entry('R2',report_value(s$main$retrieval$provenance$dataset_label,'R2 cohort'),
    'https://r2.amc.nl/',s$main$retrieval$provenance))
  if(!is.null(s$pfister))out$Pfister<-entry('Pfister','R2 Pfister / informp3','https://r2.amc.nl/',s$pfister$provenance)
  if(!is.null(s$depmap))out$DepMap<-entry('DepMap','DepMap','https://depmap.org/portal/',s$depmap$provenance)
  if(!is.null(s$hpa))out$HPA<-entry('HPA','Human Protein Atlas',s$hpa$links$profile %||% 'https://www.proteinatlas.org/',s$hpa$metadata)
  if(!is.null(s$cerebellum))out$Cerebellum<-entry('Cerebellum','Sepp, Leiss et al., Nature (2023)',
    'https://doi.org/10.1038/s41586-023-06884-x',s$cerebellum$manifest)
  if(!is.null(s$biology))out$Regulation<-entry('Regulation','TF catalogue, TRRUST and binding evidence',
    s$biology$tf$source_url %||% 'https://www.grnpedia.org/trrust/',s$biology$targets$sources)
  tf<-s$biology$targets$evidence
  if(is.data.frame(tf)&&all(c('source','source_url')%in%names(tf))) {
    rows<-unique(tf[,c('source','source_url'),drop=FALSE])
    rows<-rows[!is.na(rows$source_url)&grepl('^https?://',rows$source_url),,drop=FALSE]
    for(i in seq_len(nrow(rows))) {
      id<-paste0('TFSource',i)
      out[[id]]<-entry(id,rows$source[i],rows$source_url[i],list())
    }
  }
  if(!is.null(s$publications))out$PubMed<-entry('PubMed','NCBI PubMed','https://pubmed.ncbi.nlm.nih.gov/',s$sources$PubMed)
  out
}

report_gene_evidence <- function(s) {
  obs<-list();missing<-character();g<-s$gene
  add<-function(id,source,tab,text,type='observation',strength='Moderate evidence',statistics=list(),flags=list()) {
    obs[[id]]<<-list(id=id,gene=g,source=source,tab=tab,text=text,type=type,strength=strength,statistics=statistics,flags=flags)
  }
  a<-s$main$analysis;st<-a$statistics;d<-a$summaries
  if(!is.null(st)&&report_number(st$p_value)&&is.data.frame(d)&&all(c('subgroup','median','n')%in%names(d))) {
    ok<-is.finite(d$median)&is.finite(d$n)&d$n>0;v<-d[ok,,drop=FALSE]
    if(nrow(v)>1) {
      top<-paste(v$subgroup[v$median==max(v$median)],collapse=', ')
      sig<-st$p_value<.05
      add('r2_subgroups','R2','MB Subgroups',sprintf('%s has its highest observed subgroup median in %s in %s (Kruskal-Wallis raw p = %.3g; n = %s%s). The omnibus test does not establish which subgroup pairs differ.',
        g,top,report_value(s$main$retrieval$provenance$dataset_label,'the analysed R2 cohort'),st$p_value,report_value(st$n),
        if(report_number(st$epsilon_squared))sprintf('; epsilon-squared = %.3g',st$epsilon_squared)else ''),
        'statistical association',if(sig)'Strong evidence (within this cohort)'else'Moderate evidence',as.list(st),list(difference=sig,highest=top))
    }
  } else missing<-c(missing,'R2 subgroup statistics were unavailable.')
  if(!is.null(s$pfister)&&is.data.frame(s$pfister$summaries)&&nrow(s$pfister$summaries))
    add('pfister','Pfister','Pediatric Pan-Cancer',paste(g,'has a separate descriptive comparison in the Pfister pan-cancer cohort; its expression scale is not directly comparable with the main R2 cohort.'),statistics=list(summary=s$pfister$summaries))
  for(ep in c('metastasis','survival')) {
    z<-s$main$clinical[[ep]]$all
    if(identical(z$status,'complete')) {
      v<-z$statistics;p<-if(ep=='survival')v$logrank_p else v$p_value
      if(report_number(p))add(ep,'R2',if(ep=='survival')'Survival'else'Metastasis',
        if(ep=='survival')sprintf('%s overall-survival association uses the %s cutoff: raw log-rank p = %.3g; n = %s. This is an exploratory, unadjusted association.',g,report_value(v$cutoff_method),p,report_value(v$n)) else
        sprintf('%s expression by metastatic status has raw p = %.3g; rank-biserial effect = %s. This is exploratory and unadjusted for clinical covariates.',g,p,report_value(v$rank_biserial_metastatic_vs_m0)),
        'statistical association','Moderate evidence',as.list(v))
    }
  }
  dm<-s$depmap$mb_models
  if(is.data.frame(dm)&&'dependency'%in%names(dm)&&any(is.finite(dm$dependency))) {
    v<-dm$dependency[is.finite(dm$dependency)];n<-length(v);k<-sum(v<=-.5)
    add('dependency','DepMap','DepMap',sprintf('%s has gene effect <= -0.5 in %d of %d measured medulloblastoma models (range %.3g to %.3g). This descriptive screening threshold is not a significance test or evidence of clinical benefit.',g,k,n,min(v),max(v)),
      statistics=list(n=n,n_below_threshold=k,threshold=-.5,min=min(v),max=max(v),release=s$depmap$provenance),flags=list(dependent=k>0,none=k==0))
  } else missing<-c(missing,'Medulloblastoma DepMap dependency evidence was unavailable; this is not a negative result.')
  cor<-s$depmap$correlation
  if(report_number(cor$rho)&&report_number(cor$p_value))add('dependency_correlation','DepMap','DepMap',
    sprintf('Matched medulloblastoma models show an exploratory expression/gene-effect Spearman correlation of %.3g (raw p = %.3g; n = %s). Lower gene effect indicates greater dependency; correlation does not establish causality.',cor$rho,cor$p_value,report_value(cor$n)),
    'correlation','Moderate evidence',cor)
  for(species in c('human','mouse')) {
    z<-s$cerebellum$species[[species]]
    if(identical(z$status,'complete')&&report_number(z$n_measured)&&z$n_measured>0&&report_number(z$n_expressing)) {
      frac<-z$n_expressing/z$n_measured
      add(paste0('development_',species),'Cerebellum','snRNA-seq Developing Cerebellum',sprintf('%s exonic transcripts were detected in %s of %s measured %s cerebellar nuclei (%.2f%%). Detection depends on sequencing depth and does not establish biological absence.',z$gene$symbol,z$n_expressing,z$n_measured,species,100*frac),
        statistics=list(n=z$n_measured,positive=z$n_expressing,fraction=frac,gene_id=z$gene$id),flags=list(low_detection=frac<=.01))
    } else missing<-c(missing,paste(tools::toTitleCase(species),'developing-cerebellum expression evidence was unavailable.'))
  }
  h<-s$hpa
  if(is.null(h))missing<-c(missing,'Protein-level evidence was unavailable for this gene in the current report.') else {
    # HPA RNA measurements are not protein evidence. Preserve source modality.
    loc<-hpa_value(h$raw$json[['Subcellular location']])
    if(length(loc)&&nzchar(loc))add('hpa_localisation','HPA','ProteinAtlas.org Profile',paste(g,'has the following HPA subcellular annotation:',loc,'. This annotation is not a matched tumour protein abundance measurement.'),statistics=list(location=loc,release=h$metadata$release))
    missing<-c(missing,'Matched medulloblastoma RNA/protein abundance evidence was unavailable; no RNA/protein concordance is inferred from HPA RNA or localisation annotations.')
  }
  tf<-s$biology$targets$evidence
  if(is.data.frame(tf)&&nrow(tf)&&'evidence_type'%in%names(tf)) {
    ncur<-sum(tf$evidence_type=='CURATED TARGET',na.rm=TRUE);nbind<-sum(tf$evidence_type=='TF BINDING',na.rm=TRUE)
    add('regulation','Regulation','TF Targets',sprintf('%s has %d curated regulatory evidence records and %d binding-associated records. Curated records retain their references and contexts; binding alone does not demonstrate regulation in medulloblastoma.',g,ncur,nbind),
      'regulatory evidence','Moderate evidence',list(curated_records=ncur,binding_records=nbind,references=unique(tf$reference)))
  }
  list(gene=g,version=REPORT_EVIDENCE_VERSION,observations=obs,missing=unique(missing),citations=report_citations(s))
}

report_detect_patterns <- function(e) {
  o<-e$observations;out<-list();r<-o$r2_subgroups;d<-o$dependency;h<-o$development_human
  add<-function(id,title,text,refs)out[[id]]<<-list(id=id,title=title,text=text,evidence=refs,strength='Hypothesis-generating observation',type='cross-dataset hypothesis')
  if(isTRUE(r$flags$difference)&&isTRUE(h$flags$low_detection))add('tumour_development','Tumour / development',
    'Subgroup-associated tumour expression co-occurs with detection in at most 1% of measured human developmental nuclei. This supports exploring tumour/developmental context differences, but the unmatched platforms do not establish tumour-specific expression or a cell of origin.',c(r$id,h$id))
  if(isTRUE(r$flags$difference)&&isTRUE(d$flags$dependent))add('expression_dependency','Expression / dependency',
    'Tumour subgroup expression differences coexist with dependency-screen signals in medulloblastoma models. This suggests context-dependent biology worth investigating; the patient groups and models are not matched and expression is not shown to cause dependency.',c(r$id,d$id))
  if(isTRUE(r$flags$difference)&&isTRUE(d$flags$none))add('expression_dependency_discordance','Expression / dependency discordance',
    'Despite subgroup expression differences in tumours, none of the measured medulloblastoma models meets the report gene-effect threshold. This discordance may reflect context dependence; it does not rule out a biological role.',c(r$id,d$id))
  m<-o$development_mouse
  if(!is.null(h)&&!is.null(m)&&isTRUE(h$flags$low_detection)!=isTRUE(m$flags$low_detection))add('species_difference','Developmental detection differs by species',
    'The human gene and its one-to-one mouse orthologue fall on different sides of the descriptive 1% detection threshold. Species, sampled stages, cell composition and sequencing depth may contribute; this is not evidence of evolutionary divergence.',c(h$id,m$id))
  out
}

report_interpretation <- function(s) {
  e<-report_gene_evidence(s);p<-report_detect_patterns(e)
  # Each sentence retains its observation ID; do not pad sparse evidence to a word target.
  priority<-c('r2_subgroups','dependency','development_human','hpa_localisation','regulation')
  chosen<-head(intersect(priority,names(e$observations)),3)
  summary<-lapply(chosen,function(id)list(text=e$observations[[id]]$text,evidence=id))
  if(length(p))summary<-c(summary,list(list(text=p[[1]]$text,evidence=p[[1]]$evidence)))
  list(evidence=e,summary=summary,patterns=p,
    fallback=if(!length(p))'No strong cross-dataset pattern was identified from the currently available analyses.'else NULL)
}

report_hpa_selected_plot <- function(profile,section,metric=NULL) {
  d<-profile[[section]]$expression
  if(!is.data.frame(d)||!nrow(d))return(NULL)
  choices<-unique(paste(d$assay,d$metric,sep=' | '))
  preferred<-switch(section,brain='humanBrainRegional | nTPM',tissue='consensusTissue | nTPM',choices[1])
  if(is.null(metric))metric<-if(preferred%in%choices)preferred else choices[1]
  if(!metric%in%choices)return(NULL)
  d<-d[paste(d$assay,d$metric,sep=' | ')==metric&is.finite(d$expression),,drop=FALSE]
  if(!nrow(d))return(NULL)
  d$category<-as.character(d$label)
  for(field in c('lineage','context'))if(field%in%names(d)) {
    present<-!is.na(d[[field]])&nzchar(d[[field]])
    d$category[present]<-d[[field]][present]
  }
  organ_order<-c('Brain','Eye','Endocrine tissues','Respiratory system',
    'Proximal digestive tract','Gastrointestinal tract','Liver & Gallbladder',
    'Pancreas','Kidney & Urinary bladder','Male tissues','Female tissues',
    'Muscle & Vascular tissue','Connective & Soft tissue','Skin','Bone marrow & Lymphoid tissues')
  label_order<-tolower(c('Cerebral cortex','Cerebellum','Basal ganglia','Hypothalamus',
    'Midbrain','Amygdala','Choroid plexus','Hippocampal formation','Spinal cord','Retina',
    'Thyroid gland','Parathyroid gland','Adrenal gland','Pituitary gland','Lung',
    'Salivary gland','Esophagus','Tongue','Stomach','Duodenum','Small intestine','Colon',
    'Rectum','Liver','Gallbladder','Pancreas','Kidney','Urinary bladder','Testis',
    'Epididymis','Seminal vesicle','Prostate','Vagina','Ovary','Fallopian tube',
    'Endometrium','Cervix','Cervix, uterine','Placenta','Breast','Heart muscle',
    'Smooth muscle','Skeletal muscle','Adipose tissue','Skin','Appendix','Spleen',
    'Lymph node','Tonsil','Bone marrow','Thymus'))
  groups<-c(organ_order,setdiff(unique(d$category),organ_order))
  # Keep every measurement, grouping source organ/lineage annotations together.
  # Unlisted structures follow listed structures within their source group.
  d<-d[order(match(d$category,groups),match(tolower(d$label),label_order),seq_len(nrow(d)),na.last=TRUE),,drop=FALSE]
  d$category<-factor(d$category,levels=unique(d$category))
  # Each source row has its own position: repeated labels must never be summed.
  d$position<-seq_len(nrow(d))
  colors<-hpa_category_colors(as.character(unique(d$category)))
  p<-ggplot2::ggplot(d,ggplot2::aes(position,expression,fill=category))+
    ggplot2::geom_col(width=.8)+
    ggplot2::scale_x_continuous(breaks=d$position,labels=d$label,expand=ggplot2::expansion(add=.6))+
    ggplot2::scale_fill_manual(values=colors,name='Organ / cell group')+
    ggplot2::labs(x=NULL,y=unique(d$metric)[1],title=paste(profile$identity$canonical_symbol,section),
      subtitle=paste(nrow(d),'available values |',metric))+
    ggplot2::theme_classic(base_size=11)+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=60,hjust=1),legend.position='bottom')+
    ggplot2::guides(fill=ggplot2::guide_legend(ncol=4))
  attr(p,'hpa_width')<-max(10,nrow(d)*.23)
  p
}

hpa_category_colors <- function(labels) {
  fixed<-c('Brain'='#F5BE28','Endocrine tissues'='#CEC7E0',
    'Bone marrow & Lymphoid tissues'='#D52C2C','Connective & Soft tissue'='#EC862B',
    'Muscle & Vascular tissue'='#805441','Female tissues'='#906CAB',
    'Gastrointestinal tract'='#347EB8','Male tissues'='#2D9261',
    'Proximal digestive tract'='#EF9799','Liver & Gallbladder'='#499DA0',
    'Kidney & Urinary bladder'='#F3BCD8','Respiratory system'='#8BBC86',
    'Pancreas'='#584599','Eye'='#EBDF51','Skin'='#41B6C7')
  # A label's fallback colour is stable across genes, subsets and sub-tabs.
  palette<-grDevices::hcl.colors(36,'Dynamic')
  out<-vapply(labels,function(label) {
    if(label%in%names(fixed))return(unname(fixed[label]))
    palette[1+sum(utf8ToInt(label)*seq_along(utf8ToInt(label)))%%length(palette)]
  },character(1))
  stats::setNames(out,labels)
}

report_figure_catalog <- function(s,materialize=FALSE) {
  out<-list();g<-s$gene
  add<-function(id,title,source,tab,plot=NULL,detail='',kind=NULL) {
    if(is.null(plot)&&is.null(kind))return()
    out[[id]]<<-list(id=id,title=title,source=source,tab=tab,plot=plot,kind=kind,
      selected=!identical(s$figure_selection[[id]],FALSE),caption=paste0(g,': ',title,'. ',detail))
  }
  stats<-function(z) {
    if(is.null(z))return('')
    keys<-intersect(c('test','n','low_n','high_n','p_value','epsilon_squared','logrank_p','hazard_ratio_high_vs_low','cutoff_method','cutoff'),names(z))
    labels<-c(test='Test',n='n',low_n='Low n',high_n='High n',p_value='raw p',epsilon_squared='epsilon-squared',
      logrank_p='raw log-rank p',hazard_ratio_high_vs_low='HR High/Low',cutoff_method='cutoff method',cutoff='cutoff')
    paste(vapply(keys,function(k) {
      value<-z[[k]]
      if(is.numeric(value))value<-format(value,digits=if(k=='cutoff')10 else 4,trim=TRUE)
      paste(labels[[k]],paste(value,collapse=', '))
    },character(1)),collapse='; ')
  }
  add('pfister_pan_cancer','Pfister / informp3 pan-cancer expression','Pfister','Pediatric Pan-Cancer',s$pfister$plot,'Descriptive log2(1 + FPKM); no cross-platform abundance comparison.')
  add('r2_subgroups','R2 medulloblastoma subgroup expression','R2','MB Subgroups',s$main$analysis$plot,paste(report_value(s$main$retrieval$provenance$dataset_label),stats(s$main$analysis$statistics)))
  add('r2_subtypes','R2 molecular-subtype expression','R2','MB Subtypes',s$main$subtypes$plot,stats(s$main$subtypes$overall))
  for(ep in c('metastasis','survival'))for(co in names(s$main$clinical[[ep]])) {
    a<-s$main$clinical[[ep]][[co]]
    add(paste(ep,co,sep='_'),paste(ep,co,'R2 cohort'),'R2',if(ep=='survival')'Survival'else'Metastasis',a$plot,stats(a$statistics))
  }
  dv<-if(identical(s$settings$depmap_view,'all'))'all_models_ranked_dependency'else'dependency'
  add('depmap_dependency',paste('DepMap dependency:',if(dv=='dependency')'medulloblastoma and retinoblastoma'else'all models'),'DepMap','DepMap',s$depmap$plots[[dv]],paste('Release',report_value(s$depmap$provenance$release)))
  add('depmap_expression_dependency','DepMap matched expression and dependency','DepMap','DepMap',s$depmap$plots$expression_dependency,'Only models with both measurements are shown; correlation is not causation.')
  for(k in c('tissue','brain','single_cell','subcellular','cancer','blood','cell_line','structure','interaction')) {
    d<-s$hpa[[k]]$expression
    if(is.data.frame(d)&&nrow(d))add(paste0('hpa_',k),paste('Human Protein Atlas',k),'HPA','ProteinAtlas.org Profile',
      report_hpa_selected_plot(s$hpa,k,s$settings$hpa_metrics[[k]]),paste('Selected assay:',s$settings$hpa_metrics[[k]] %||% 'section default','; all available values; organ / cell group colours; assay units preserved.'))
  }
  for(species in c('human','mouse'))for(kind in c('annotation','expression')) {
    z<-s$cerebellum$species[[species]];cells<-s$cerebellum$cells[[species]]
    if(is.null(cells)||(kind=='expression'&&!identical(z$status,'complete')))next
    id<-paste('cerebellum',kind,species,sep='_');field<-s$settings$cerebellum_annotation %||% 'broad_lineage'
    p<-NULL
    if(materialize&&!identical(s$figure_selection[[id]],FALSE)) {
      assets<-list(species=lapply(s$cerebellum$cells,function(x)list(cells=x)))
      p<-cerebellum_export_plot(list(cells=cells),if(kind=='expression')z else NULL,paste(tools::toTitleCase(species),if(kind=='expression')z$gene$symbol else field),field,cerebellum_colors(assets,field))
    }
    add(id,paste(tools::toTitleCase(species),'developing cerebellum',kind,'UMAP'),'Cerebellum','snRNA-seq Developing Cerebellum',p,
      paste('Published coordinates; full-data export;',nrow(cells),'nuclei;',if(kind=='expression')paste(z$n_expressing,'expressing /',z$n_measured,'measured; exonic UMI, log1p colour scale.')else paste('annotation:',field)),kind)
  }
  set<-s$settings$biology_target_set %||% 'curated'
  add('tf_target_enrichment',paste('TF target enrichment:',set),'Regulation','Functional Biology',s$biology$enrichment[[set]]$plot,'Target-set enrichment is not independent tumour pathway activity evidence; see source evidence and adjusted statistics.')
  add('functional_annotation','Functional annotation','Regulation','Functional Biology',s$biology$annotation$plot,'Source library memberships; not all annotations have experimental evidence codes.')
  out
}
