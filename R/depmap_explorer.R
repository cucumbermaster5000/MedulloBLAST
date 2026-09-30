# Cell-line dependency analysis; independent of Shiny and the legacy report.
depmap_explorer_files <- function(cfg) {
  manifest_path<-cfg$depmap_explorer_manifest %||% Sys.getenv("DEPMAP_EXPLORER_MANIFEST","")
  if(nzchar(manifest_path)) {
    m<-jsonlite::fromJSON(manifest_path,simplifyVector=FALSE)
    if(is.null(m$release) || is.null(m$files$model) || is.null(m$files$dependency))stop("DepMap manifest needs release, model and dependency files")
    if(!identical(m$dependency_metric,"Chronos gene effect"))stop("This module requires documented Chronos gene effect")
    if(!is.null(m$files$expression) && !identical(m$expression_metric,"log2(TPM + 1)"))stop("Unsupported expression metric; verify release documentation")
    paths<-list();audit<-list()
    for(role in names(m$files)) {
      f<-m$files[[role]]
      if(!identical(f$release,m$release))stop("DepMap release mismatch in manifest")
      if(is.null(f$md5) || !file.exists(f$path) || unname(tools::md5sum(f$path))!=f$md5)stop("DepMap file missing or checksum/release mismatch: ",role)
      paths[[role]]<-normalizePath(f$path,winslash="/")
      audit[[role]]<-list(name=basename(f$path),md5=f$md5,path=paths[[role]],source_url=f$source_url %||% m$source_url)
    }
    return(list(paths=paths,provenance=list(release=m$release,source=m$source_url,files=audit,license=m$license,
      manifest=normalizePath(manifest_path,winslash="/"),dependency_metric=m$dependency_metric,expression_metric=m$expression_metric,
      verified_at=utc_now(),release_note="Explicit local release manifest; file checksums and matching declared releases verified."),warnings=character()))
  }
  base<-depmap_public_files(cfg)
  extra<-tryCatch(depmap_public_files(cfg,c("OmicsExpressionProteinCodingGenesTPMLogp1.csv","README.txt")),error=function(e)NULL)
  p<-list(model=unname(base$paths[["Model.csv"]]),dependency=unname(base$paths[["CRISPRGeneEffect.csv"]]))
  if(!is.null(extra))p$expression<-unname(extra$paths[["OmicsExpressionProteinCodingGenesTPMLogp1.csv"]])
  m<-jsonlite::fromJSON(base$provenance$manifest)
  audit<-lapply(p,function(path){f<-m$files[m$files$name==basename(path),];list(name=f$name,md5=f$computed_md5,path=path,source_url=f$download_url,
    file_modified_at=format(file.info(path)$mtime,tz="UTC",usetz=TRUE))})
  list(paths=p,provenance=list(release="DepMap 24Q4 Public v1",source=m$url,doi=m$doi,manifest=base$provenance$manifest,
    files=audit,license=m$license,dependency_metric="Chronos gene effect",expression_metric="log2(TPM + 1)",verified_at=utc_now(),
    release_note="Pinned archived 24Q4 v1, not the latest release. Latest public catalogue checked separately; current download URLs require browser verification. All three datasets use this same archive."),
    warnings=if(is.null(extra))"Same-release expression is unavailable; dependency analyses remain available." else character())
}

depmap_breadbox_base <- function(cfg) sub("/+$","",cfg$depmap_breadbox_url %||% "https://depmap.org/portal/breadbox")

depmap_breadbox_request <- function(cfg,path,body=NULL) {
  url<-paste0(depmap_breadbox_base(cfg),if(startsWith(path,"/"))path else paste0("/",path))
  last<-NULL
  for(attempt in seq_len(3L)) {
    response<-tryCatch({
      h<-curl::new_handle();curl::handle_setheaders(h,"User-Agent"="r2-depmap/0.4 (+local research app)",Accept="application/json")
      if(!is.null(body)) {
        curl::handle_setheaders(h,"User-Agent"="r2-depmap/0.4 (+local research app)",Accept="application/json","Content-Type"="application/json")
        curl::handle_setopt(h,customrequest="POST",postfields=jsonlite::toJSON(body,auto_unbox=TRUE,null="null"))
      }
      curl::curl_fetch_memory(url,h)
    },error=function(e){last<<-e;NULL})
    if(!is.null(response) && response$status_code>=200L && response$status_code<300L)
      return(jsonlite::fromJSON(rawToChar(response$content),simplifyVector=FALSE))
    if(!is.null(response))last<-simpleError(paste("DepMap Breadbox returned HTTP",response$status_code))
    if(attempt<3L)Sys.sleep(attempt)
  }
  stop("Official DepMap Breadbox request failed: ",conditionMessage(last),call.=FALSE)
}

depmap_breadbox_source <- function(cfg) {
  catalog<-depmap_breadbox_request(cfg,"/datasets/?feature_type=gene&sample_type=depmap_model")
  pick<-function(id) {
    hit<-Filter(function(x)identical(x$given_id,id),catalog)
    if(length(hit)!=1L)stop("Official DepMap dataset is missing or ambiguous: ",id)
    hit[[1]]
  }
  dependency<-pick("Chronos_Combined");expression<-pick("expression")
  release_of<-function(x)x$dataset_metadata$download_file_info$release_name %||% NA_character_
  releases<-c(release_of(dependency),release_of(expression))
  if(anyNA(releases) || length(unique(releases))!=1L)stop("DepMap dependency and expression releases do not match")
  file_info<-function(x)list(name=x$name,given_id=x$given_id,taiga_id=x$taiga_id,units=x$units,
    file=x$dataset_metadata$download_file_info$file_name,release=release_of(x),dataset_id=x$id)
  list(dependency=dependency,expression=expression,provenance=list(release=releases[[1]],source=depmap_breadbox_base(cfg),
    files=list(model=list(name="DepMap model metadata",given_id="depmap_model_metadata",release=releases[[1]]),
      dependency=file_info(dependency),expression=file_info(expression)),dependency_metric=dependency$units,
    expression_metric=expression$units,verified_at=utc_now(),
    release_note="Current official public DepMap dependency, expression, and model metadata queried through Breadbox."),warnings=character())
}

depmap_breadbox_models <- function(cfg) {
  fields<-c("CellLineName","OncotreeLineage","OncotreePrimaryDisease","OncotreeSubtype","OncotreeCode")
  raw<-depmap_breadbox_request(cfg,"/datasets/tabular/depmap_model_metadata",list(columns=fields))
  ids<-sort(unique(unlist(lapply(raw[intersect(fields,names(raw))],names),use.names=FALSE)))
  if(!length(ids))stop("Official DepMap model metadata returned no ModelIDs")
  value<-function(field,id) {z<-raw[[field]][[id]];if(is.null(z)||!length(z))NA_character_ else as.character(z[[1]])}
  out<-data.frame(ModelID=ids,stringsAsFactors=FALSE)
  for(field in fields)out[[field]]<-vapply(ids,function(id)value(field,id),character(1))
  out
}

depmap_breadbox_gene <- function(cfg,dataset_id,identity) {
  symbol<-identity$canonical_symbol
  raw<-depmap_breadbox_request(cfg,paste0("/datasets/matrix/",dataset_id),list(features=list(symbol),feature_identifier="label"))
  values<-raw[[symbol]]
  if(is.null(values))return(list(data=data.frame(ModelID=character(),value=numeric()),column=NA_character_,status="Gene absent from matrix",profiles_excluded=0L))
  ids<-names(values);number<-vapply(values,function(x)if(is.null(x)||!length(x))NA_real_ else suppressWarnings(as.numeric(x[[1]])),numeric(1))
  list(data=data.frame(ModelID=ids,value=number,stringsAsFactors=FALSE),column=symbol,entrez_id=identity$entrez_id,
    status="available",profiles_excluded=0L)
}

depmap_tumor_display <- function(x) {
  models<-x$tumor_models %||% data.frame()
  out<-data.frame(CellLineName=character(),ModelID=character(),cancer_type=character(),dependency=numeric(),expression=numeric(),
    dependency_rank=numeric(),total_models=numeric(),dependency_percentile=numeric(),dependency_metric=character(),source=character(),
    release=character(),mb_subgroup=character(),subgroup_source=character(),stringsAsFactors=FALSE)
  if(nrow(models))out<-data.frame(CellLineName=models$CellLineName,ModelID=models$ModelID,cancer_type=models$cancer_type,
    dependency=models$dependency,expression=models$expression,dependency_rank=x$tumor_context$dependency_rank[match(models$ModelID,x$tumor_context$ModelID)],
    total_models=x$tumor_context$total_models[match(models$ModelID,x$tumor_context$ModelID)],dependency_percentile=x$tumor_context$dependency_percentile[match(models$ModelID,x$tumor_context$ModelID)],
    dependency_metric=models$dependency_metric,source=x$provenance$source %||% NA_character_,release=models$release,
    mb_subgroup=models$mb_subgroup,subgroup_source=models$subgroup_source,stringsAsFactors=FALSE)
  out[order(out$cancer_type,out$dependency,out$ModelID,na.last=TRUE),,drop=FALSE]
}

depmap_mb_metadata <- function(models,provenance) {
  require_columns(models,c("ModelID","CellLineName","OncotreePrimaryDisease","OncotreeSubtype"))
  if(anyNA(models$ModelID) || anyDuplicated(models$ModelID))stop("Model metadata identifiers must be present and unique")
  fields<-c("OncotreePrimaryDisease","OncotreeSubtype")
  flags<-lapply(fields,function(f)!is.na(models[[f]]) & grepl("^medulloblastoma(,|$)",trimws(models[[f]]),ignore.case=TRUE))
  models$is_medulloblastoma<-Reduce(`|`,flags)
  rb_flags<-lapply(fields,function(f)!is.na(models[[f]]) & grepl("^retinoblastoma(,|$)",trimws(models[[f]]),ignore.case=TRUE))
  models$is_retinoblastoma<-Reduce(`|`,rb_flags)
  models$cancer_type<-ifelse(models$is_medulloblastoma,"Medulloblastoma",ifelse(models$is_retinoblastoma,"Retinoblastoma","Other"))
  models$mb_inclusion_reason<-vapply(seq_len(nrow(models)),function(i)paste(vapply(which(vapply(flags,`[`,logical(1),i)),function(j)paste0(fields[j],"=",models[[fields[j]]][i]),character(1)),collapse="; "),character(1))
  models$retinoblastoma_inclusion_reason<-vapply(seq_len(nrow(models)),function(i)paste(vapply(which(vapply(rb_flags,`[`,logical(1),i)),function(j)paste0(fields[j],"=",models[[fields[j]]][i]),character(1)),collapse="; "),character(1))
  models$mb_subgroup<-"Unknown";models$subgroup_source<-"Not supplied/verified in native metadata"
  models$subgroup_reference<-NA_character_;models$subgroup_evidence_type<-"Unknown"
  for(i in which(models$is_medulloblastoma)) {
    label<-paste(models$OncotreePrimaryDisease[i],models$OncotreeSubtype[i])
    patterns<-c(WNT="\\bWNT\\b",SHH="\\bSHH\\b",`Group 3`="\\bGroup[ -]?3\\b",`Group 4`="\\bGroup[ -]?4\\b")
    found<-names(patterns)[vapply(patterns,function(p)grepl(p,label,ignore.case=TRUE,perl=TRUE),logical(1))]
    if(length(found)==1L) {
      models$mb_subgroup[i]<-found;models$subgroup_source[i]<-"Native DepMap OncoTree disease/subtype metadata"
      models$subgroup_reference[i]<-provenance$doi %||% provenance$source
      models$subgroup_evidence_type[i]<-"Native metadata annotation"
    }
  }
  # Explicit project annotations requested by the user; retain their provenance.
  model_name<-toupper(gsub('[^A-Za-z0-9]','',models$CellLineName))
  group3<-models$is_medulloblastoma & (models$ModelID=='ACH-001349' | model_name%in%c('HDMB03','D283','D283MED'))
  models$mb_subgroup[group3]<-'Group 3'
  models$subgroup_source[group3]<-'User-specified project annotation: HD-MB03 and D283 as Group 3'
  models$subgroup_reference[group3]<-NA_character_
  models$subgroup_evidence_type[group3]<-'Project annotation'
  models$release<-provenance$release
  models
}

depmap_gene_identity <- function(cfg,gene,dependency_path=NULL) {
  original<-gene;symbol<-biology_symbol(gene)
  x<-tryCatch(resolve_publication_gene(cfg,gene),error=function(e)NULL)
  if(!is.null(x))return(list(original_query=original,canonical_symbol=x$canonical_symbol,entrez_id=x$gene_id,
    source="NCBI human Gene (same identity resolver as Publications)",warnings=x$warnings))
  if(is.null(dependency_path))return(list(original_query=original,canonical_symbol=symbol,entrez_id=NA_character_,
    source="Exact submitted official symbol; NCBI identity service unavailable",warnings="Current nomenclature could not be verified; the exact submitted symbol was queried in DepMap."))
  mapping<-depmap_gene_columns(dependency_path);hit<-mapping[toupper(mapping$gene)==symbol,,drop=FALSE]
  if(nrow(hit)!=1L)stop("Gene identity unavailable or ambiguous; use an official symbol and retry")
  id<-if(grepl("\\([0-9]+\\)$",hit$column))sub("^.*\\(([0-9]+)\\)$","\\1",hit$column) else NA_character_
  list(original_query=original,canonical_symbol=hit$gene,entrez_id=id,source="Unique exact symbol in pinned DepMap release header; NCBI unavailable",
    warnings="Current nomenclature could not be verified; exact release symbol used. No alias guessed.")
}

depmap_gene_vector <- function(path,identity,expression=FALSE) {
  columns<-depmap_header(path);idcol<-depmap_id_column(columns)
  ids<-ifelse(grepl("\\([0-9]+\\)$",columns),sub("^.*\\(([0-9]+)\\)$","\\1",columns),NA_character_)
  byid<-which(!is.na(ids) & !is.na(identity$entrez_id) & ids==identity$entrez_id)
  bysymbol<-which(toupper(clean_symbol(columns))==toupper(identity$canonical_symbol))
  hit<-if(length(byid))byid else bysymbol
  if(length(hit)>1L)stop("Ambiguous DepMap matrix gene mapping")
  if(!length(hit))return(list(data=data.frame(ModelID=character(),value=numeric()),column=NA_character_,status="Gene absent from matrix",profiles_excluded=0L))
  if(!length(byid) && !is.na(identity$entrez_id) && !is.na(ids[hit]) && ids[hit]!=identity$entrez_id)stop("Gene symbol/Entrez identifier conflict in release")
  extra<-if(expression)intersect(c("ProfileID","is_default_entry"),columns) else character()
  d<-data.table::fread(path,select=unique(c(idcol,extra,columns[hit])),data.table=FALSE,check.names=FALSE,showProgress=FALSE)
  excluded<-0L
  if(expression && "is_default_entry" %in% names(d)) {
    keep<-!is.na(d$is_default_entry) & tolower(as.character(d$is_default_entry)) %in% c("true","1")
    excluded<-sum(!keep);d<-d[keep,,drop=FALSE]
  }
  if(anyNA(d[[idcol]]) || anyDuplicated(d[[idcol]]))stop("Model IDs must be unique after default-profile selection; no arbitrary replicate averaging")
  values<-suppressWarnings(as.numeric(d[[columns[hit]]]))
  if(any(!is.na(d[[columns[hit]]]) & is.na(values)))stop("Non-numeric values in DepMap gene column")
  out<-data.frame(ModelID=as.character(d[[idcol]]),value=values)
  if("ProfileID" %in% names(d))out$ProfileID<-d$ProfileID
  list(data=out,column=columns[hit],entrez_id=ids[hit],status="available",profiles_excluded=excluded)
}

depmap_rank_context <- function(data) {
  d<-data[is.finite(data$dependency),,drop=FALSE];n<-nrow(d)
  d$dependency_rank<-rank(d$dependency,ties.method="min")
  d$total_models<-rep(n,n)
  # Strictly weaker means a larger score. Ties receive the same rank/percentile.
  d$dependency_percentile<-if(n)100*(n-rank(d$dependency,ties.method="max"))/n else numeric()
  d[order(d$dependency,d$ModelID),,drop=FALSE]
}

depmap_descriptive <- function(values,cohort,total) {
  x<-values[is.finite(values)];n<-length(x)
  data.frame(cohort=cohort,total_models=total,n=n,mean=if(n)mean(x) else NA_real_,median=if(n)stats::median(x) else NA_real_,
    sd=if(n>1)stats::sd(x) else NA_real_,iqr=if(n)stats::IQR(x) else NA_real_,minimum=if(n)min(x) else NA_real_,maximum=if(n)max(x) else NA_real_)
}

analyze_depmap_profile <- function(models,dependency,expression,identity,provenance) {
  models<-depmap_mb_metadata(models,provenance)
  allids<-union(models$ModelID,dependency$data$ModelID)
  d<-models[match(allids,models$ModelID),,drop=FALSE];d$missing_model_metadata<-is.na(d$ModelID);d$ModelID<-allids
  d$CellLineName[is.na(d$CellLineName) | !nzchar(d$CellLineName)]<-d$ModelID[is.na(d$CellLineName) | !nzchar(d$CellLineName)]
  d$dependency<-dependency$data$value[match(d$ModelID,dependency$data$ModelID)]
  d$expression<-expression$data$value[match(d$ModelID,expression$data$ModelID)]
  d$dependency_missing<-!is.finite(d$dependency);d$expression_missing<-!is.finite(d$expression)
  d$gene<-identity$canonical_symbol;d$release<-provenance$release
  d$dependency_metric<-provenance$dependency_metric;d$expression_metric<-provenance$expression_metric
  d$mb_subgroup[is.na(d$mb_subgroup)]<-"Unknown"
  mb<-d[which(d$is_medulloblastoma %in% TRUE),,drop=FALSE]
  rb<-d[which(d$is_retinoblastoma %in% TRUE),,drop=FALSE]
  tumor<-d[which(d$is_medulloblastoma %in% TRUE | d$is_retinoblastoma %in% TRUE),,drop=FALSE]
  pan<-depmap_rank_context(d);mb_context<-pan[which(pan$is_medulloblastoma %in% TRUE),,drop=FALSE]
  rb_context<-pan[which(pan$is_retinoblastoma %in% TRUE),,drop=FALSE]
  tumor_context<-pan[which(pan$is_medulloblastoma %in% TRUE | pan$is_retinoblastoma %in% TRUE),,drop=FALSE]
  stats<-rbind(depmap_descriptive(d$dependency,"All DepMap models",nrow(d)),depmap_descriptive(mb$dependency,"Medulloblastoma",nrow(mb)),
    depmap_descriptive(rb$dependency,"Retinoblastoma",nrow(rb)))
  matched<-mb[is.finite(mb$dependency) & is.finite(mb$expression),,drop=FALSE]
  correlation<-list(n=nrow(matched),method="Spearman, two-sided; asymptotic p-value (ties permitted)",rho=NA_real_,p_value=NA_real_,status="Insufficient matched MB models for correlation (minimum 5)",
    dependency_metric=provenance$dependency_metric,expression_metric=provenance$expression_metric)
  if(nrow(matched)>=5 && length(unique(matched$dependency))>1 && length(unique(matched$expression))>1) {
    test<-stats::cor.test(matched$expression,matched$dependency,method="spearman",exact=FALSE)
    correlation$rho<-unname(test$estimate);correlation$p_value<-test$p.value;correlation$status<-"Exploratory association; correlation does not establish causality"
  } else if(nrow(matched)>=5)correlation$status<-"Constant expression or dependency; correlation is undefined"
  warnings<-c(identity$warnings,if(!nrow(mb))"No models have explicit medulloblastoma disease metadata.",
    if(!nrow(mb_context))"No medulloblastoma models have valid queried-gene dependency values.",
    if(!nrow(rb))"No models have explicit retinoblastoma disease metadata.",
    if(nrow(rb) && !nrow(rb_context))"No retinoblastoma models have valid queried-gene dependency values.",
    if(dependency$status!="available")dependency$status,if(expression$status!="available")paste("Expression:",expression$status),
    if(any(d$missing_model_metadata))paste(sum(d$missing_model_metadata),"dependency models lack matching metadata; retained in pan-DepMap context."),
    "Subgroup comparisons are not performed. Subgroups use native metadata plus explicit project Group 3 annotations for HD-MB03 and D283; no subgroup is inferred from marker genes.")
  list(queried_gene=identity$original_query,canonical_symbol=identity$canonical_symbol,identity=identity,provenance=provenance,
    mapping=list(dependency_column=dependency$column,expression_column=expression$column,expression_profiles_excluded=expression$profiles_excluded),
    model_metadata=models,mb_models=mb,retinoblastoma_models=rb,tumor_models=tumor,all_models=pan,mb_context=mb_context,retinoblastoma_context=rb_context,tumor_context=tumor_context,matched=matched,statistics=stats,correlation=correlation,
    counts=list(mb_models=nrow(mb),mb_dependency=sum(is.finite(mb$dependency)),mb_expression=sum(is.finite(mb$expression)),mb_matched=nrow(matched),mb_verified_subgroup=sum(mb$mb_subgroup!="Unknown"),
      retinoblastoma_models=nrow(rb),retinoblastoma_dependency=sum(is.finite(rb$dependency)),retinoblastoma_expression=sum(is.finite(rb$expression)),all_dependency=nrow(pan)),
    strongest=mb_context[mb_context$dependency==min(c(mb_context$dependency,Inf)),,drop=FALSE],
    weakest=mb_context[mb_context$dependency==max(c(mb_context$dependency,-Inf)),,drop=FALSE],
    rank_definition="Rank 1 is the lowest (most negative) Chronos score. Ties use minimum rank.",
    percentile_definition="100 x number of valid DepMap models with strictly larger (weaker) scores / total valid models. Equal scores are not counted as weaker. Higher means relatively stronger dependency; the maximum is below 100%.",
    mb_filter="Trimmed OncotreePrimaryDisease OR OncotreeSubtype matches ^medulloblastoma(,|$), case-insensitive; all models, including HD-MB03, come from the same official release.",
    retinoblastoma_filter="Trimmed OncotreePrimaryDisease OR OncotreeSubtype matches ^retinoblastoma(,|$), case-insensitive; no cell-line-name matching.",
    warnings=warnings,retrieved_at=utc_now())
}

plot_depmap_profile <- function(x) {
  tumor<-depmap_tumor_display(x);tumor<-tumor[is.finite(tumor$dependency),,drop=FALSE];pan<-x$all_models;matched<-x$matched
  base<-ggplot2::theme_classic(base_size=12);p<-list()
  colors<-c(Medulloblastoma="#b43b47",Retinoblastoma="#2369a8")
  if(nrow(tumor))p$dependency<-ggplot2::ggplot(tumor,ggplot2::aes(dependency,reorder(CellLineName,dependency),color=cancer_type))+
    ggplot2::geom_vline(xintercept=0,linetype=3,color="grey50")+ggplot2::geom_point(size=3)+base+
    ggplot2::scale_color_manual(values=colors)+
    ggplot2::labs(title=paste(x$canonical_symbol,"dependency in medulloblastoma and retinoblastoma models"),subtitle=x$provenance$release,x="DepMap gene effect",y="Model",color="Cancer type",caption="Lower scores indicate stronger in-vitro dependency. Every highlighted model uses the same official DepMap release.")
  if(nrow(pan))p$all_models_ranked_dependency<-ggplot2::ggplot(pan,ggplot2::aes(dependency_rank,dependency))+
    ggplot2::geom_point(color="#b8c1c8",size=1.2,alpha=.65)+ggplot2::geom_hline(yintercept=0,linetype=3,color="grey50")+
    ggplot2::geom_point(data=x$tumor_context,ggplot2::aes(color=cancer_type),size=3)+
    ggrepel::geom_text_repel(data=x$tumor_context,ggplot2::aes(label=CellLineName,color=cancer_type),seed=17,max.overlaps=Inf,max.iter=20000,max.time=3,min.segment.length=0,size=3.5,box.padding=.6)+base+
    ggplot2::scale_color_manual(values=colors)+
    ggplot2::labs(title=paste(x$canonical_symbol,"dependency across all DepMap models"),subtitle=paste(x$provenance$release,"| N =",nrow(pan)),x="Rank (1 = strongest dependency; ties share minimum rank)",y="Chronos gene effect",color="Cancer type",caption="Every current-release score is included. Red marks medulloblastoma and blue marks retinoblastoma.")
  if(nrow(matched))p$expression_dependency<-ggplot2::ggplot(matched,ggplot2::aes(expression,dependency,color=mb_subgroup))+
    ggplot2::scale_color_manual(values=category_colors(unique(matched$mb_subgroup)))+
    ggplot2::geom_point(size=3)+ggrepel::geom_text_repel(ggplot2::aes(label=CellLineName),seed=17,max.overlaps=Inf,size=3.5)+base+
    ggplot2::labs(title=paste(x$canonical_symbol,"expression and dependency in MB models"),
      subtitle=paste("N =",nrow(matched),"| Spearman rho =",format(x$correlation$rho,digits=3),"| p =",format(x$correlation$p_value,digits=3)),
      x=x$provenance$expression_metric,y="Chronos gene effect",color="MB subgroup",caption="Same-release expression and dependency joined by ModelID. HD-MB03 and D283: project Group 3 annotations. Exploratory association, not causation.")
  p
}

get_depmap_profile <- function(cfg,gene,refresh=FALSE) {
  biology_symbol(gene)
  source<-depmap_breadbox_source(cfg)
  identity<-depmap_gene_identity(cfg,gene)
  key<-paste0("depmap_explorer_v5_",object_digest(list(source$provenance$release,source$provenance$files,identity$canonical_symbol,identity$entrez_id)))
  lock<-filelock::lock(file.path(cfg$cache_dir,paste0(key,".fetch.lock")),timeout=120000)
  if(is.null(lock))stop("DepMap gene analysis is busy; retry later")
  on.exit(filelock::unlock(lock),add=TRUE)
  cached<-if(!refresh)cache_get(cfg,key) else NULL
  if(!is.null(cached)){cached$queried_gene<-gene;cached$identity<-identity;return(cached)}
  models<-depmap_breadbox_models(cfg)
  tumor_list<-depmap_mb_metadata(models,source$provenance)
  tumor_list<-tumor_list[tumor_list$is_medulloblastoma | tumor_list$is_retinoblastoma,,drop=FALSE]
  list_path<-file.path(cfg$cache_dir,paste0("depmap_mb_retinoblastoma_models_",object_digest(source$provenance$files$model),".csv"))
  metadata_lock<-filelock::lock(paste0(list_path,".lock"),timeout=10000)
  if(is.null(metadata_lock))stop("DepMap metadata cache is busy")
  tryCatch({if(!file.exists(list_path))data.table::fwrite(tumor_list,list_path,na="NA")},finally=filelock::unlock(metadata_lock))
  dependency<-depmap_breadbox_gene(cfg,source$dependency$given_id,identity)
  expression<-tryCatch(depmap_breadbox_gene(cfg,source$expression$given_id,identity),error=function(e)
    list(data=data.frame(ModelID=character(),value=numeric()),column=NA_character_,profiles_excluded=0L,status=conditionMessage(e)))
  x<-analyze_depmap_profile(models,dependency,expression,identity,source$provenance)
  x$warnings<-unique(c(x$warnings,source$warnings));x$provenance$tumor_model_list_file<-normalizePath(list_path,winslash="/")
  x$plots<-plot_depmap_profile(x);x$artifact<-cache_path(cfg,key)
  if(expression$status=="available")cache_put(cfg,key,x) else {
    x$artifact<-tempfile("depmap_partial_",tmpdir=cfg$cache_dir,fileext=".rds")
    saveRDS(x,x$artifact)
  }
  x
}

write_depmap_download <- function(x,kind,file) {
  display<-depmap_tumor_display(x)
  tables<-list(dependency=display[is.finite(display$dependency),,drop=FALSE],expression=x$tumor_models[is.finite(x$tumor_models$expression),,drop=FALSE],
    matched_expression_dependency=x$matched,model_metadata=x$tumor_models,subgroup_evidence=x$mb_models[,c("ModelID","CellLineName","mb_subgroup","subgroup_source","subgroup_reference","subgroup_evidence_type","release")],
    all_models_dependency=x$all_models,MB_context_ranks=x$mb_context,tumor_context_ranks=x$tumor_context,statistics=x$statistics)
  if(kind %in% names(tables)) {if(!nrow(tables[[kind]]))stop("No applicable data for this download");data.table::fwrite(tables[[kind]],file,na="NA")} else if(kind=="result.rds")saveRDS(x,file) else {
    plot<-sub("_plot[.](pdf|png)$","",kind);format<-sub("^.*[.]","",kind)
    if(is.null(x$plots[[plot]]) || !format %in% c("pdf","png"))stop("No applicable plot for this download")
    ggplot2::ggsave(file,x$plots[[plot]],device=format,width=13,height=8,dpi=300,limitsize=FALSE)
  }
  invisible(file)
}
