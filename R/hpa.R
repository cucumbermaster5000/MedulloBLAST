# Official HPA structured records; no HTML scraping or invented biological facts.
hpa_fetch <- function(url) {
  h<-curl::new_handle(timeout=60,connecttimeout=15,useragent="MedulloblastomaGeneExplorer/0.5")
  r<-curl::curl_fetch_memory(url,h)
  if(r$status_code!=200L)stop("HPA endpoint unavailable")
  rawToChar(r$content)
}
hpa_value <- function(x) {
  if(is.null(x)||!length(x))return("")
  if(is.list(x) && !is.null(names(x)))return(paste(vapply(names(x),function(n)paste0(n,": ",hpa_value(x[[n]])),character(1)),collapse="; "))
  paste(as.character(unlist(x)),collapse="; ")
}
hpa_fields <- function(record,pattern) {
  fields<-grep(pattern,names(record),value=TRUE,ignore.case=TRUE)
  values<-vapply(record[fields],hpa_value,character(1))
  keep<-!is.na(values)&nzchar(values)&values!="NA"
  data.frame(field=fields[keep],value=unname(values[keep]),source=rep("Human Protein Atlas",sum(keep)),stringsAsFactors=FALSE)
}
hpa_parse <- function(raw_xml,record,identity,endpoint,retrieved_at=utc_now()) {
  doc<-xml2::read_xml(if(is.raw(raw_xml))raw_xml else charToRaw(raw_xml),options="NONET")
  entry<-xml2::xml_find_all(doc,"/proteinAtlas/entry")
  if(length(entry)!=1L)stop("HPA entry is ambiguous")
  txt<-function(node,path)xml2::xml_text(xml2::xml_find_first(node,path))
  ensembl<-xml2::xml_attr(xml2::xml_find_first(entry,"identifier"),"id")
  symbol<-txt(entry,"name")
  if(!identical(toupper(symbol),toupper(identity$canonical_symbol)) || !identical(ensembl,record$Ensembl))stop("HPA identity mismatch")
  release<-xml2::xml_attr(entry,"version")
  identity$ensembl_id<-ensembl
  identity$entrez_id<-xml2::xml_attr(xml2::xml_find_first(entry,"identifier/xref[@db='NCBI GeneID']"),"id")
  if(!is.na(identity$entrez_id) && !is.null(identity$gene_id) && identity$gene_id!=identity$entrez_id)stop("HPA/NCBI identifier mismatch")
  identity$uniprot_id<-xml2::xml_attr(xml2::xml_find_all(entry,"identifier/xref[contains(@db,'Uniprot')]"),"id")
  identity$aliases<-unique(c(identity$aliases_used,identity$aliases_excluded,xml2::xml_text(xml2::xml_find_all(entry,"synonym"))))
  identity$hgnc_id<-record[["HGNC ID"]]
  if(is.null(identity$official_name)||!nzchar(identity$official_name))identity$official_name<-record[["Gene description"]]
  rows<-list()
  for(assay in xml2::xml_find_all(entry,"rnaExpression")) {
    kind<-xml2::xml_attr(assay,"assayType")
    # Exclude animal data rather than silently treating it as human expression.
    if(grepl("^(mouse|pig)",kind))next
    for(d in xml2::xml_find_all(assay,"data")) {
      label<-xml2::xml_find_first(d,"tissue|immuneCell|cellLine|cellType")
      for(level in xml2::xml_find_all(d,"level[@expRNA]"))rows[[length(rows)+1L]]<-data.frame(
        label=xml2::xml_text(label),context=xml2::xml_attr(label,"organ"),
        lineage=xml2::xml_attr(label,"lineage"),assay=kind,
        expression=suppressWarnings(as.numeric(xml2::xml_attr(level,"expRNA"))),
        metric=xml2::xml_attr(level,"unitRNA"),level_type=xml2::xml_attr(level,"type"),
        source=xml2::xml_attr(assay,"source"),stringsAsFactors=FALSE)
    }
  }
  rna<-if(length(rows))do.call(rbind,rows) else data.frame()
  expression<-function(kinds)if(nrow(rna))rna[rna$assay %in% kinds,,drop=FALSE] else data.frame()
  section<-function(pattern,kinds=character())list(annotations=hpa_fields(record,pattern),expression=expression(kinds),details=data.frame())
  out<-list(status="complete",identity=identity,metadata=list(source="Human Protein Atlas",release=release,
    retrieved_at=retrieved_at,ensembl_id=ensembl,canonical_symbol=symbol,original_query=identity$original_query,
    endpoints=endpoint,expression_metrics=if(nrow(rna))unique(rna$metric) else character(),cache_hit=FALSE),
    general=hpa_fields(record,"^(Gene|Ensembl|Uniprot|Chromosome|Position|Protein class|Biological process|Molecular function|Disease involvement|Evidence|HPA evidence|UniProt evidence|Subcellular|Secretome|Interactions)"),
    tissue=section("^(RNA tissue|Protein tissue|Reliability \\(IH)",c("consensusTissue","tissue")),
    brain=section("^RNA (brain|single nuclei brain)",c("humanBrainRegional","humanBrain")),
    single_cell=section("^RNA single cell"),subcellular=section("^(Subcellular|Reliability \\(IF|CCD|Antibody)"),
    cancer=section("^(RNA cancer|Cancer prognostics)"),blood=section("^(RNA blood|Blood|Secretome)",c("immuneCell","immuneCellLineage")),
    cell_line=section("^(RNA cell line|Cell line)","cellLine"),structure=section("^(Uniprot|Protein length|Molecular mass|Signal peptide|Transmembrane)"),
    interaction=section("^Interactions"),warnings=character())
  # Preserve all available IHC levels and their tissue/cell-type context.
  ihc<-lapply(xml2::xml_find_all(entry,"tissueExpression/data"),function(d) {
    levels<-xml2::xml_find_all(d,"level")
    if(!length(levels))return(NULL)
    data.frame(tissue=txt(d,"tissue"),cell_type=txt(d,"tissueCell"),level=xml2::xml_text(levels),
      type=xml2::xml_attr(levels,"type"),reliability=txt(entry,"tissueExpression/verification"),source="HPA IHC")
  })
  out$tissue$details<-if(length(Filter(Negate(is.null),ihc)))do.call(rbind,ihc) else data.frame()
  chains<-xml2::xml_find_all(entry,"proteinstructure/structure/chain")
  if(length(chains))out$structure$details<-do.call(rbind,lapply(chains,function(n){
    parent<-xml2::xml_parent(n)
    data.frame(protein_id=xml2::xml_attr(n,"ensembl_peptide_id"),transcript_id=xml2::xml_attr(n,"ensembl_transcript_id"),
      length_aa=xml2::xml_attr(n,"length"),method=xml2::xml_attr(parent,"method"),structure_url=xml2::xml_attr(parent,"url"),source="HPA")
  }))
  # JSON supplies selected single-cell expression, not a complete cell atlas.
  for(field in grep("^RNA single cell type specific (nCPM|nTPM)$",names(record),value=TRUE)) {
    values<-record[[field]]
    if(length(values)&&!is.null(names(values)))out$single_cell$expression<-data.frame(label=names(values),
      expression=as.numeric(unlist(values)),metric=sub(".* ","",field),assay="HPA specificity-selected cell types",source="HPA")
  }
  out$links<-list(profile=paste0("https://www.proteinatlas.org/",ensembl),
    ncbi=paste0("https://www.ncbi.nlm.nih.gov/gene/",identity$gene_id))
  out$raw<-list(xml=raw_xml,json=record)
  out
}
get_hpa_profile <- function(cfg,gene,refresh=FALSE) {
  original<-gene;symbol<-biology_symbol(gene)
  hcfg<-cfg;hcfg$cache_dir<-file.path(cfg$cache_dir,"hpa");dir.create(hcfg$cache_dir,recursive=TRUE,showWarnings=FALSE)
  index<-cache_key("hpa-v1",symbol)
  cached<-if(!refresh)cache_get(hcfg,index) else NULL
  if(!is.null(cached) && pubmed_fresh(list(retrieved_at=cached$metadata$retrieved_at),168)) {
    cached$identity$original_query<-original;cached$metadata$original_query<-original;cached$metadata$cache_hit<-TRUE;return(cached)
  }
  identity<-tryCatch(resolve_publication_gene(cfg,gene,refresh),error=function(e)NULL)
  tryCatch({
    canonical<-if(is.null(identity))symbol else identity$canonical_symbol
    search_url<-paste0("https://www.proteinatlas.org/api/search_download.php?search=",utils::URLencode(canonical,reserved=TRUE),"&format=json&columns=g,gs,eg&compress=no")
    records<-jsonlite::fromJSON(hpa_fetch(search_url),simplifyVector=FALSE)
    if(!is.null(records$Gene))records<-list(records)
    hits<-Filter(function(r)identical(toupper(r$Gene),canonical),records)
    if(length(hits)!=1L||!grepl("^ENSG[0-9]+$",hits[[1]]$Ensembl))stop("HPA identity not uniquely resolved")
    ens<-hits[[1]]$Ensembl;base<-paste0("https://www.proteinatlas.org/",ens)
    raw_xml<-hpa_fetch(paste0(base,".xml"));record<-jsonlite::fromJSON(hpa_fetch(paste0(base,".json")),simplifyVector=FALSE)
    if(is.null(record$Gene)&&length(record)==1L)record<-record[[1]]
    if(is.null(identity))identity<-list(original_query=original,canonical_symbol=canonical,official_name=record[["Gene description"]],source="Human Protein Atlas")
    x<-hpa_parse(raw_xml,record,identity,c(search_url,paste0(base,c(".xml",".json"))))
    cache_put(hcfg,cache_key("hpa-v1",c(canonical,ens,x$metadata$release)),x)
    cache_put(hcfg,index,x);x
  },error=function(e)list(status="unavailable",identity=identity,metadata=list(original_query=original,retrieved_at=utc_now()),
    warnings="Human Protein Atlas data are currently unavailable."))
}
hpa_section_table <- function(x,section) {
  if(section=="summary")return(x$general)
  s<-x[[section]];if(is.null(s))return(data.frame())
  chunks<-Filter(function(d)is.data.frame(d)&&nrow(d)>0,s[c("annotations","expression","details")])
  if(!length(chunks))return(data.frame())
  data.table::rbindlist(chunks,fill=TRUE)
}
