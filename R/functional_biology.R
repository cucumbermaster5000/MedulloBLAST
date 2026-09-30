# Source-backed functional biology; all functions are independent of Shiny.
biology_libraries <- function() c(GO_BP="GO_Biological_Process_2026",GO_MF="GO_Molecular_Function_2026",
  GO_CC="GO_Cellular_Component_2026",Reactome="Reactome_Pathways_2024",KEGG="KEGG_2021_Human",WikiPathways="WikiPathways_2024_Human")

biology_symbol <- function(gene) {
  if(length(gene)!=1L || is.na(gene) || nchar(trimws(gene))>64L || !grepl("^[A-Za-z][A-Za-z0-9._-]*$",trimws(gene))) stop("Provide one gene symbol")
  toupper(trimws(gene))
}

# A content-addressed raw snapshot accompanies every external response.
biology_resource <- function(cfg,key,url,refresh=FALSE,post_genes=NULL) {
  key <- paste0("biology_raw_v1_",object_digest(list(key,url,post_genes)))
  cached <- if(!refresh) cache_get(cfg,key) else NULL
  if(!is.null(cached)) return(cached)
  if(is.null(post_genes)) raw <- enrichr_get(url) else {
    h <- curl::new_handle(timeout=120,connecttimeout=20,followlocation=TRUE)
    curl::handle_setform(h,list=paste(post_genes,collapse="\n"),description="Human TF target evidence; Gene Explorer")
    response <- curl::curl_fetch_memory(url,handle=h)
    if(response$status_code %in% c(429,500,502,503,504)) {
      Sys.sleep(1)
      response <- curl::curl_fetch_memory(url,handle=h)
    }
    if(response$status_code!=200) stop("Enrichr list submission failed: ",response$status_code)
    raw <- rawToChar(response$content)
  }
  hash <- object_digest(raw)
  folder <- file.path(cfg$cache_dir,"biology_raw");dir.create(folder,showWarnings=FALSE)
  path <- file.path(folder,paste0(hash,".txt"));writeBin(charToRaw(raw),path)
  x <- list(raw=raw,provenance=list(url=url,retrieved_at=utc_now(),content_hash=hash,raw_file=path))
  cache_put(cfg,key,x);x
}

biology_library <- function(cfg,library,refresh=FALSE) {
  catalog <- list_enrichr_libraries(cfg,refresh=refresh)
  if(!library %in% catalog$libraries$libraryName) stop("Pinned Enrichr library is unavailable: ",library)
  # Reuse the validated GMT parser and cache; preserve raw GMT for this module.
  url <- paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=",utils::URLencode(library,reserved=TRUE))
  raw <- biology_resource(cfg,paste0("library:",library),url,refresh)
  list(sets=parse_gmt(raw$raw),provenance=c(raw$provenance,list(library=library,catalog_retrieved_at=catalog$provenance$retrieved_at)))
}

functional_memberships <- function(gene,inputs) {
  empty <- data.frame(queried_gene=character(),category=character(),library=character(),term=character(),
    term_identifier=character(),evidence_type=character(),source=character(),source_url=character(),retrieved_at=character())
  rows <- lapply(names(inputs),function(category) {
    x <- inputs[[category]];if(is.null(x$sets))return(NULL)
    terms <- names(x$sets)[vapply(x$sets,function(s)gene %in% s,logical(1))]
    if(!length(terms))return(NULL)
    ids <- regmatches(terms,regexpr("GO:[0-9]+|R-HSA-[0-9]+|WP[0-9]+",terms))
    id <- rep(NA_character_,length(terms));has <- grepl("GO:[0-9]+|R-HSA-[0-9]+|WP[0-9]+",terms);id[has]<-ids
    data.frame(queried_gene=gene,category=category,library=x$provenance$library,term=terms,term_identifier=id,
      evidence_type=if(startsWith(category,"GO_")) "GENE ANNOTATION" else "PATHWAY MEMBERSHIP",
      source="Enrichr library membership (source annotation evidence codes not supplied)",source_url=x$provenance$url,retrieved_at=x$provenance$retrieved_at)
  })
  rows <- Filter(Negate(is.null),rows)
  if(length(rows))do.call(rbind,rows) else empty
}

get_functional_annotation <- function(cfg,gene,refresh=FALSE) {
  gene <- biology_symbol(gene);libraries <- biology_libraries();inputs <- list();warnings <- character()
  for(category in names(libraries)) inputs[[category]] <- tryCatch(biology_library(cfg,libraries[[category]],refresh),error=function(e) {
    message("[Functional annotation] ",conditionMessage(e));warnings<<-c(warnings,paste("Enrichr is currently unavailable for",libraries[[category]]));list()
  })
  signatures <- lapply(inputs,function(x)x$provenance$content_hash)
  key <- paste0("biology_annotation_v1_",object_digest(list(gene,libraries,signatures)))
  cached <- if(!refresh && !length(warnings)) cache_get(cfg,key) else NULL
  if(!is.null(cached))return(cached)
  rows <- functional_memberships(gene,inputs)
  summary <- if(nrow(rows)) paste0(gene," has ",nrow(rows)," documented memberships in the selected libraries. Examples: ",
    paste(head(unique(rows$term),3),collapse="; "),". These are annotations, not evidence of pathway activation.") else
    if(length(warnings)) "Functional annotation could not be fully retrieved." else "No memberships were found in the selected libraries; this does not imply no biological function."
  plot <- NULL
  if(nrow(rows)) {
    counts <- as.data.frame(table(rows$category));names(counts)<-c("category","annotations")
    plot <- ggplot2::ggplot(counts,ggplot2::aes(reorder(category,annotations),annotations))+
      ggplot2::geom_col(fill="#216b75")+ggplot2::coord_flip()+ggplot2::theme_classic(base_size=12)+
      ggplot2::labs(title=paste("Documented biological annotations for",gene),x="Annotation category",y="Number of memberships",
        caption="Descriptive annotation counts, not enrichment statistics or independent biological effects.")
  }
  out <- list(gene=gene,status=if(length(warnings)) "partial" else "complete",rows=rows,summary=summary,plot=plot,
    libraries=libraries,provenance=lapply(inputs,function(x)x$provenance),warnings=warnings)
  if(!length(warnings))cache_put(cfg,key,out)
  out
}

biology_gene_audit <- function(genes) {
  raw <- as.character(genes);normalized <- toupper(trimws(raw))
  valid <- !is.na(normalized) & grepl("^[A-Z][A-Z0-9._-]*$",normalized)
  duplicate <- valid & duplicated(normalized)
  data.frame(original=raw,normalized=normalized,submitted=valid & !duplicate,
    reason=ifelse(!valid,"Invalid/empty symbol",ifelse(duplicate,"Duplicate normalized symbol","")))
}

empty_biology_enrichment <- function() {
  data.frame(library=character(),rank=integer(),term=character(),p_value=numeric(),adjusted_p_value=numeric(),
    api_z_score=numeric(),odds_ratio=numeric(),combined_score=numeric(),overlap_n=integer(),contributing_genes=character())
}

parse_biology_enrichment <- function(raw,library) {
  payload <- jsonlite::fromJSON(raw,simplifyVector=FALSE)
  if(!library %in% names(payload))stop("Enrichr response omitted the requested library")
  if(!length(payload[[library]]))return(empty_biology_enrichment())
  do.call(rbind,lapply(payload[[library]],function(row) {
    if(length(row)<7)stop("Enrichr result schema changed")
    genes <- unlist(row[[6]],use.names=FALSE)
    data.frame(library=library,rank=as.integer(row[[1]]),term=as.character(row[[2]]),p_value=as.numeric(row[[3]]),
      adjusted_p_value=as.numeric(row[[7]]),api_z_score=as.numeric(row[[4]]),odds_ratio=NA_real_,
      combined_score=as.numeric(row[[5]]),overlap_n=length(genes),contributing_genes=paste(genes,collapse=";"))
  }))
}

enrich_tf_gene_set <- function(cfg,genes,tf,target_set,refresh=FALSE) {
  audit <- biology_gene_audit(genes);submitted <- audit$normalized[audit$submitted]
  empty <- list(status="insufficient",definition=target_set,submitted_genes=submitted,input_audit=audit,
    accepted_genes=character(),rows=empty_biology_enrichment(),plot=NULL,warnings="At least five unique supported targets are required; no enrichment performed.")
  if(length(submitted)<5)return(empty)
  libraries <- biology_libraries();libraries <- libraries[names(libraries)!="GO_CC"]
  allrows <- list();inputs <- list();warnings <- character();provenance <- list()
  add <- biology_resource(cfg,list("target-list",tf,target_set,submitted),"https://maayanlab.cloud/Enrichr/addList",refresh,submitted)
  id <- jsonlite::fromJSON(add$raw)$userListId
  if(length(id)!=1L || is.na(id))stop("Enrichr did not return a gene-list identifier")
  view <- biology_resource(cfg,list("accepted-list",tf,target_set,id),paste0("https://maayanlab.cloud/Enrichr/view?userListId=",id),refresh)
  accepted <- unique(toupper(unlist(jsonlite::fromJSON(view$raw)$genes,use.names=FALSE)))
  audit$service_retained <- audit$normalized %in% accepted
  audit$service_note <- ifelse(audit$submitted & !audit$service_retained,"Not returned by Enrichr view",ifelse(audit$service_retained,"Retained in service list; annotation recognition is library-specific",audit$reason))
  if(length(accepted)<5) {empty$accepted_genes<-accepted;empty$input_audit<-audit;return(empty)}
  for(category in names(libraries)) {
    lib <- libraries[[category]]
    allrows[[category]] <- tryCatch({
      source <- biology_library(cfg,lib,refresh)
      recognized <- accepted %in% unique(unlist(source$sets,use.names=FALSE))
      # Recorded recognition is library membership, not evidence that absent symbols are invalid genes.
      inputs[[category]] <- data.frame(gene=accepted,library=lib,recognized_in_library=recognized)
      url <- paste0("https://maayanlab.cloud/Enrichr/enrich?userListId=",id,"&backgroundType=",utils::URLencode(lib,reserved=TRUE))
      response <- biology_resource(cfg,list("target-enrichment",tf,target_set,submitted,lib,source$provenance$content_hash),url,refresh)
      tab <- parse_biology_enrichment(response$raw,lib)
      tab$category <- rep(category,nrow(tab));tab$tf <- rep(tf,nrow(tab));tab$target_set <- rep(target_set,nrow(tab))
      tab$target_set_size <- rep(length(submitted),nrow(tab));tab$evidence_type <- rep("TARGET-GENE ENRICHMENT",nrow(tab))
      tab$retrieved_at <- rep(response$provenance$retrieved_at,nrow(tab));tab$source_url <- rep(url,nrow(tab))
      provenance[[category]] <- response$provenance
      tab
    },error=function(e) {message("[Enrichr target enrichment] ",conditionMessage(e));warnings<<-c(warnings,paste("Enrichr is currently unavailable for",lib));NULL})
  }
  rows <- if(length(Filter(Negate(is.null),allrows)))do.call(rbind,allrows) else empty_biology_enrichment()
  service_status <- if(length(warnings) && !nrow(rows)) "unavailable" else if(length(warnings)) "partial" else "complete"
  plot <- NULL
  if(nrow(rows)) {
    rows <- rows[order(rows$adjusted_p_value,rows$p_value,-rows$combined_score),]
    show <- head(rows[is.finite(rows$adjusted_p_value) & rows$adjusted_p_value<=.05,],15)
    if(nrow(show)) {
      show$label <- paste(show$category,show$term,sep=": ");show$label<-factor(show$label,levels=rev(unique(show$label)))
      plot <- ggplot2::ggplot(show,ggplot2::aes(combined_score,label,size=overlap_n,color=adjusted_p_value))+
        ggplot2::geom_point()+ggplot2::scale_color_gradient(low="#216b75",high="#d18a40")+ggplot2::theme_classic(base_size=11)+
        ggplot2::theme(plot.title.position="plot")+
        ggplot2::scale_y_discrete(labels=function(x)vapply(x,function(s)paste(strwrap(s,70),collapse="\n"),character(1)))+
        ggplot2::labs(title=paste(tf,"target-gene enrichment"),subtitle=paste(strwrap(target_set,100),collapse="\n"),x="Enrichr combined score",y=NULL,
          size="Overlap genes",color="Adjusted p",caption="Top 15 terms by Enrichr adjusted p <= 0.05. Target enrichment does not prove pathway activation.")
    } else warnings <- c(warnings,"No TF-target pathways met the predefined enrichment significance threshold.")
  }
  if(!nrow(rows) && service_status=="complete")warnings<-c(warnings,"No TF-target pathways were returned by the selected enrichment libraries.")
  list(status=service_status,
    definition=target_set,submitted_genes=submitted,accepted_genes=accepted,input_audit=audit,
    library_recognition=if(length(inputs))do.call(rbind,inputs) else data.frame(),rows=rows,plot=plot,
    provenance=list(add_list=add$provenance,view=view$provenance,libraries=provenance),warnings=warnings,
    method="Enrichr API default background; service adjusted p-values (BH per library) preserved. No across-library/target-set correction. API fourth score retained as api_z_score per official client; odds ratio not supplied separately.")
}
