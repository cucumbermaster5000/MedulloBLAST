classify_human_tf <- function(cfg,gene,refresh=FALSE) {
  gene <- biology_symbol(gene)
  url <- "https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv"
  out <- list(gene=gene,status="UNKNOWN",family=NA_character_,identifier=NA_character_,source="Lambert human TF catalogue v1.01; PMID 29425488",source_url=url,warnings=character())
  tryCatch({
    source <- biology_resource(cfg,"Lambert-full-v1.01",url,refresh)
    d <- data.table::fread(text=source$raw,data.table=FALSE,check.names=FALSE)
    require_columns(d,c("HGNC symbol","Is TF?","DBD","Ensembl ID"))
    if(nrow(d)<1000)stop("Incomplete TF catalogue")
    row <- d[toupper(d[["HGNC symbol"]])==gene,,drop=FALSE]
    values <- unique(row[["Is TF?"]])
    out$status <- if(!nrow(row)) "FALSE" else if(identical(values,"Yes")) "TRUE" else if(identical(values,"No")) "FALSE" else "UNKNOWN"
    out$family <- if(nrow(row))paste(unique(row$DBD),collapse="; ") else NA_character_
    out$identifier <- if(nrow(row))paste(unique(row[["Ensembl ID"]]),collapse="; ") else NA_character_
    out$evidence <- if(nrow(row))row else data.frame()
    out$provenance <- source$provenance
    out$note <- if(!nrow(row)) "Not listed as a TF in this selected catalogue; not a universal assertion of absence of regulatory function." else "Catalogue assessment retained; DBD is the reported family/domain."
    out
  },error=function(e) {message("[TF classification] ",conditionMessage(e));out$warnings<-"Transcription-factor status could not be determined reliably.";out})
}

empty_tf_evidence <- function() data.frame(tf=character(),target_gene=character(),source=character(),evidence_type=character(),
  experiment=character(),context=character(),direction=character(),score=numeric(),reference=character(),source_url=character(),retrieved_at=character())

get_supported_tf_targets <- function(cfg,gene,classification=NULL,refresh=FALSE) {
  gene <- biology_symbol(gene)
  if(is.null(classification))classification<-classify_human_tf(cfg,gene,refresh)
  evidence <- empty_tf_evidence();warnings <- character();sources <- list()
  if(classification$status!="TRUE") return(list(gene=gene,classification=classification,evidence=evidence,summary=data.frame(),
    target_sets=list(curated=character(),binding=character()),sources=sources,warnings=if(classification$status=="FALSE")
      "This gene is not classified as a transcription factor in the selected reference resource." else "Transcription-factor status could not be determined reliably."))
  for(lib in c("ChEA_2022","ENCODE_TF_ChIP-seq_2015")) {
    rows <- tryCatch({
      x <- biology_library(cfg,lib,refresh);sources[[lib]]<-x$provenance
      labels <- names(x$sets)
      tf_token <- sub("[_ ].*$","",labels)
      human <- if(lib=="ChEA_2022")grepl("(^|[_ ])HUMAN($|[_ ])",toupper(labels)) else grepl("(^| )hg19($| )",labels)
      terms <- labels[toupper(tf_token)==gene & human]
      if(!length(terms))return_rows <- NULL else return_rows <- data.table::rbindlist(lapply(terms,function(term) {
        genes <- x$sets[[term]]
        tokens <- strsplit(term,"[_ ]+")[[1]]
        pmid <- tokens[grepl("^[0-9]{7,9}$",tokens)]
        data.frame(tf=gene,target_gene=genes,source=lib,evidence_type="TF BINDING",experiment=term,
          context=term,direction="Not supplied; binding does not establish regulation",score=NA_real_,
          reference=if(length(pmid))paste(pmid,collapse=";") else NA_character_,source_url=x$provenance$url,retrieved_at=x$provenance$retrieved_at)
      }),fill=TRUE)
      return_rows
    },error=function(e) {message("[TF binding source] ",conditionMessage(e));warnings<<-c(warnings,paste("Target source unavailable:",lib));NULL})
    if(!is.null(rows))evidence<-rbind(evidence,as.data.frame(rows))
  }
  curated <- tryCatch({
    url <- "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
    x <- biology_resource(cfg,"TRRUST-v2-human",url,refresh);sources$TRRUST_v2<-x$provenance
    d <- data.table::fread(text=x$raw,header=FALSE,data.table=FALSE)
    if(ncol(d)!=4 || nrow(d)<1000)stop("TRRUST schema changed")
    names(d)<-c("tf","target","direction","pmids")
    d<-d[d$tf==gene,,drop=FALSE]
    if(!nrow(d))NULL else do.call(rbind,lapply(seq_len(nrow(d)),function(i) {
      refs<-strsplit(d$pmids[i],";",fixed=TRUE)[[1]]
      data.frame(tf=gene,target_gene=d$target[i],source="TRRUST_v2_human",evidence_type="CURATED TARGET",
        experiment=paste0("PMID:",refs),context="Not supplied in TRRUST edge download",direction=d$direction[i],score=NA_real_,
        reference=refs,source_url=url,retrieved_at=x$provenance$retrieved_at)
    }))
  },error=function(e) {message("[TRRUST] ",conditionMessage(e));warnings<<-c(warnings,"Curated target source TRRUST is currently unavailable.");NULL})
  if(!is.null(curated))evidence<-rbind(evidence,curated)
  evidence <- unique(evidence)
  summary <- if(nrow(evidence))do.call(rbind,lapply(split(evidence,evidence$target_gene),function(d)
    data.frame(target_gene=d$target_gene[1],sources=paste(sort(unique(d$source)),collapse=";"),source_count=length(unique(d$source)),
      evidence_records=nrow(d),source_experiments=length(unique(paste(d$source,d$experiment))),
      binding=any(d$evidence_type=="TF BINDING"),curated=any(d$evidence_type=="CURATED TARGET")))) else
    data.frame(target_gene=character(),sources=character(),source_count=integer(),evidence_records=integer(),source_experiments=integer(),binding=logical(),curated=logical())
  rownames(summary)<-NULL
  list(gene=gene,classification=classification,evidence=evidence,summary=summary,sources=sources,
    target_sets=list(curated=sort(unique(evidence$target_gene[evidence$evidence_type=="CURATED TARGET"])),
      binding=sort(unique(evidence$target_gene[evidence$evidence_type=="TF BINDING"]))),
    warnings=c(warnings,"Binding is not proof of regulation. Sources/experiments can overlap and are not necessarily independent. No perturbation or predicted-target library is used. Evidence is not assumed to be medulloblastoma-specific."))
}

get_gene_biology <- function(cfg,gene,refresh=FALSE) {
  gene<-biology_symbol(gene)
  latest_key<-paste0("biology_current_v3_",gene)
  cached<-if(!refresh)cache_get(cfg,latest_key) else NULL
  if(!is.null(cached) && is.finite(as.numeric(difftime(Sys.time(),as.POSIXct(cached$retrieved_at,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC"),units="hours"))) &&
      as.numeric(difftime(Sys.time(),as.POSIXct(cached$retrieved_at,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC"),units="hours"))<24)return(cached)
  annotation<-get_functional_annotation(cfg,gene,refresh)
  classification<-classify_human_tf(cfg,gene,refresh)
  targets<-get_supported_tf_targets(cfg,gene,classification,refresh)
  enrichment<-list()
  definitions<-c(curated="Human TRRUST curated regulatory targets; union across references; directions retained in evidence, not mixed to infer activation",
    binding="Human ChEA 2022 and ENCODE hg19 binding-associated targets; union across experiments; binding only")
  if(classification$status=="TRUE") for(set in names(definitions)) {
    enrichment[[set]]<-tryCatch(enrich_tf_gene_set(cfg,targets$target_sets[[set]],gene,definitions[[set]],refresh),error=function(e) {
      message("[Target enrichment] ",conditionMessage(e))
      list(status="unavailable",definition=definitions[[set]],rows=empty_biology_enrichment(),plot=NULL,
        submitted_genes=targets$target_sets[[set]],warnings="Enrichr is currently unavailable; target evidence remains available.")
    })
  }
  out<-list(gene=gene,annotation=annotation,tf=classification,targets=targets,enrichment=enrichment,
    retrieved_at=utc_now(),default_target_set="curated",libraries=biology_libraries())
  key<-paste0("biology_result_v1_",object_digest(list(gene,out$retrieved_at)))
  cache_put(cfg,key,out)
  out$artifact<-cache_path(cfg,key)
  if(annotation$status=="complete" && classification$status!="UNKNOWN" && !any(grepl("unavailable",targets$warnings,fixed=TRUE)) &&
      !any(vapply(enrichment,function(a)a$status %in% c("unavailable","partial"),logical(1))))cache_put(cfg,latest_key,out)
  out
}

write_biology_download <- function(result,kind,file,target_set="curated") {
  if(kind=="functional_annotation") data.table::fwrite(result$annotation$rows,file,na="NA") else
  if(kind %in% paste0(names(biology_libraries()),"_annotation")) {
    category<-sub("_annotation$","",kind);data.table::fwrite(result$annotation$rows[result$annotation$rows$category==category,],file,na="NA")
  } else if(kind=="TF_target_evidence") data.table::fwrite(result$targets$evidence,file,na="NA") else
  if(kind %in% c("TF_targets","TF_target_summary")) data.table::fwrite(result$targets$summary,file,na="NA") else {
    a<-result$enrichment[[target_set]]
    if(is.null(a))stop("No target enrichment for this gene")
    if(kind=="TF_target_enrichment_plot.pdf")ggplot2::ggsave(file,a$plot,device=grDevices::pdf,width=12,height=9) else
    if(kind=="TF_target_enrichment_plot.png")ggplot2::ggsave(file,a$plot,device="png",width=12,height=9,dpi=300) else
    if(kind=="TF_target_input_audit")data.table::fwrite(a$input_audit %||% data.frame(gene=a$submitted_genes),file,na="NA") else
    if(kind=="TF_target_library_recognition")data.table::fwrite(a$library_recognition %||% data.frame(gene=character(),library=character(),recognized_in_library=logical()),file,na="NA") else
    if(startsWith(kind,"TF_target_enrichment_")) {
      category<-sub("TF_target_enrichment_","",kind)
      rows<-a$rows;if(nrow(rows))rows<-rows[rows$category==category,,drop=FALSE]
      data.table::fwrite(rows,file,na="NA")
    } else stop("Unknown biological download")
  }
  invisible(file)
}
