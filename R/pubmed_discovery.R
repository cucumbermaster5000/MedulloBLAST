# Relevance-first discovery, separate from the legacy newest-paper report.
pubmed_fresh <- function(x,hours=24) {
  if(is.null(x))return(FALSE)
  age<-as.numeric(difftime(Sys.time(),as.POSIXct(x$retrieved_at,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC"),units="hours"))
  length(age)==1L && is.finite(age) && age>=0 && age<hours
}

pubmed_request <- function(cfg,endpoint,params) {
  if(!endpoint %in% c("esearch.fcgi","esummary.fcgi","efetch.fcgi"))stop("Unsupported NCBI endpoint")
  dir.create(cfg$cache_dir,recursive=TRUE,showWarnings=FALSE)
  params$tool<-"medulloblastoma_gene_explorer"
  for(pair in list(c("email","NCBI_EMAIL"),c("api_key","NCBI_API_KEY"))) {
    value<-Sys.getenv(pair[2]);if(nzchar(value))params[[pair[1]]]<-value
  }
  # All discovery workers sharing this cache serialize requests at <=2/sec.
  lock<-filelock::lock(file.path(cfg$cache_dir,"pubmed-discovery.lock"),timeout=180000)
  if(is.null(lock))stop("PubMed request queue is busy; please retry")
  on.exit(filelock::unlock(lock),add=TRUE)
  stamp<-file.path(cfg$cache_dir,"pubmed-discovery-last-request.rds")
  for(attempt in 1:3) {
    last<-if(file.exists(stamp))tryCatch(readRDS(stamp),error=function(e)0) else 0
    delay<-.55-(as.numeric(Sys.time())-last);if(is.finite(delay) && delay>0)Sys.sleep(delay)
    saveRDS(as.numeric(Sys.time()),stamp)
    h<-curl::new_handle(timeout=60,connecttimeout=15,postfields=encode_form(params),
      useragent="medulloblastoma-gene-explorer/1.0")
    curl::handle_setheaders(h,"Content-Type"="application/x-www-form-urlencoded")
    response<-tryCatch(curl::curl_fetch_memory(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",endpoint),h),
      error=function(e)NULL) # Never surface request bodies or credentials in errors.
    if(!is.null(response) && response$status_code==200L)return(rawToChar(response$content))
    retry<-is.null(response) || response$status_code %in% c(429L,500L,502L,503L,504L)
    if(!retry || attempt==3L)stop("NCBI request unavailable",if(!is.null(response))paste0(" (HTTP ",response$status_code,")"),call.=FALSE)
    Sys.sleep(attempt*2)
  }
}

resolve_publication_gene <- function(cfg,gene,refresh=FALSE) {
  original<-gene;symbol<-biology_symbol(gene)
  key<-paste0("pubmed_identity_v1_",symbol)
  x<-if(!refresh)cache_get(cfg,key) else NULL
  if(pubmed_fresh(x,168)) {x$original_query<-original;return(x)}
  raw_search<-pubmed_request(cfg,"esearch.fcgi",list(db="gene",term=paste0(symbol,"[Gene Name] AND 9606[Taxonomy ID]"),retmax=100,retmode="json"))
  s<-r2_json(raw_search,"NCBI Gene search")$esearchresult
  if(is.null(s$count) || is.null(s$idlist) || !is.null(s$errorlist))stop("NCBI Gene identity search could not be validated")
  if(!length(s$idlist))stop("Gene identity not found in NCBI human Gene; use an official symbol")
  raw_summary<-pubmed_request(cfg,"esummary.fcgi",list(db="gene",id=paste(s$idlist,collapse=","),retmode="json"))
  records<-r2_json(raw_summary,"NCBI Gene summary")$result
  records<-lapply(as.character(s$idlist),function(id)records[[id]])
  records<-Filter(function(r)!is.null(r$name) && identical(as.character(r$organism$taxid),"9606"),records)
  aliases<-function(r)trimws(strsplit(r$otheraliases %||% "",",",fixed=TRUE)[[1]])
  official<-Filter(function(r)!is.null(r$nomenclaturesymbol) && nzchar(r$nomenclaturesymbol) && toupper(r$nomenclaturesymbol)==symbol,records)
  exact<-if(length(official))official else Filter(function(r)toupper(r$name)==symbol,records)
  hits<-if(length(exact))exact else Filter(function(r)symbol %in% toupper(aliases(r)),records)
  if(length(hits)!=1L || (!length(exact) && as.integer(s$count)>100L))stop("Ambiguous human gene symbol; use an official canonical symbol")
  hit<-hits[[1]]
  canonical<-hit$nomenclaturesymbol %||% hit$name
  if(!nzchar(canonical))canonical<-hit$name
  name<-hit$nomenclaturename %||% hit$description
  if(is.null(name) || !nzchar(name))name<-hit$description
  a<-unique(aliases(hit));used<-a[nchar(a)>=3 & grepl("^[A-Za-z][A-Za-z0-9._-]*$",a)]
  x<-list(original_query=original,canonical_symbol=toupper(canonical),official_name=name,gene_id=as.character(hit$uid),
    aliases_used=setdiff(used,canonical),aliases_excluded=setdiff(a,used),source="NCBI Gene, human (taxid 9606)",
    retrieved_at=utc_now(),warnings=if(symbol!=toupper(canonical))paste("Input",original,"resolved to",canonical,"by NCBI Gene; the mapping is explicit.") else character(),
    raw=list(search=raw_search,summary=raw_summary))
  cache_put(cfg,key,x);x
}

publication_query <- function(identity,scope="general") {
  scope<-match.arg(scope,c("general","medulloblastoma"))
  quote_term<-function(x)paste0('"',gsub('["\\\r\n]',' ',x),'"[Title/Abstract]')
  symbols<-unique(c(identity$canonical_symbol,identity$aliases_used))
  context_terms<-c("gene","protein","expression","transcription","mutation","signaling","signalling","enzyme","receptor","transporter","oncogene","amplification","amplified","overexpression","kinase","factor","homeobox")
  biological<-paste0("(",paste(paste0(context_terms,"[Title/Abstract]"),collapse=" OR "),")")
  symbol_clause<-function(symbol) {
    if(grepl("^[A-Za-z]{1,3}$",symbol))return(paste0("(",paste(paste0('"',symbol," ",context_terms,'"[Title/Abstract:~3]'),collapse=" OR "),")"))
    quote_term(symbol)
  }
  # Symbol/alias mentions need biological context; the verified full name is an alternative.
  q<-paste0("((",paste(vapply(symbols,symbol_clause,character(1)),collapse=" OR "),") AND ",biological,")")
  if(!is.null(identity$official_name) && nzchar(identity$official_name))q<-paste0("(",q," OR ",quote_term(identity$official_name),")")
  if(scope=="medulloblastoma")q<-paste0(q,' AND ("medulloblastoma"[Title/Abstract] OR "Medulloblastoma"[MeSH Terms])')
  q
}

parse_publication_records <- function(raw) {
  doc<-xml2::read_xml(charToRaw(raw),options="NONET")
  if(xml2::xml_name(doc)!="PubmedArticleSet" || length(xml2::xml_find_all(doc,"//ERROR")))stop("Malformed PubMed EFetch response")
  nodes<-xml2::xml_find_all(doc,"/PubmedArticleSet/PubmedArticle | /PubmedArticleSet/PubmedBookArticle")
  lapply(nodes,function(node) {
    book<-xml2::xml_name(node)=="PubmedBookArticle"
    citation<-if(book)"./BookDocument" else "./MedlineCitation"
    article<-if(book)citation else paste0(citation,"/Article")
    value<-function(path)xml_value(node,path)
    pmid<-value(paste0(citation,"/PMID"));if(is.na(pmid) || !grepl("^[0-9]+$",pmid))stop("PubMed record has no valid PMID")
    title<-value(paste0(article,"/ArticleTitle"));if(is.na(title) && book)title<-value(paste0(article,"/Book/BookTitle"))
    authors<-xml2::xml_find_all(node,paste0(article,"/AuthorList/Author"))
    names<-vapply(authors,function(a) {
      collective<-xml_value(a,"./CollectiveName");if(!is.na(collective))return(collective)
      parts<-c(xml_value(a,"./ForeName"),xml_value(a,"./LastName"));paste(parts[!is.na(parts)],collapse=" ")
    },character(1));names<-names[nzchar(names)]
    doi<-value(if(book)'./PubmedBookData/ArticleIdList/ArticleId[@IdType="doi"]' else './PubmedData/ArticleIdList/ArticleId[@IdType="doi"]')
    if(is.na(doi))doi<-value(paste0(article,'/ELocationID[@EIdType="doi"]'))
    if(!is.na(doi) && !grepl("^10\\.[0-9]{4,9}/[^[:space:]]+$",doi))doi<-NA_character_
    date<-pubmed_date(xml2::xml_find_first(node,paste0(article,if(book)"/Book/PubDate" else "/Journal/JournalIssue/PubDate")))
    electronic<-pubmed_date(xml2::xml_find_first(node,paste0(article,'/ArticleDate[@DateType="Electronic"]')))
    year<-if(!is.na(date) && grepl("[12][0-9]{3}",date))regmatches(date,regexpr("[12][0-9]{3}",date)) else NA_character_
    blocks<-xml2::xml_find_all(node,paste0(article,"/Abstract/AbstractText"))
    abstract<-paste(vapply(blocks,function(b){label<-xml2::xml_attr(b,"Label");paste0(if(!is.na(label))paste0(label,": ") else "",xml2::xml_text(b))},character(1)),collapse="\n\n")
    types<-xml2::xml_text(xml2::xml_find_all(node,paste0(article,"/PublicationTypeList/PublicationType")))
    notices<-xml2::xml_find_all(node,paste0(citation,'/CommentsCorrectionsList/CommentsCorrections[@RefType="RetractionIn" or @RefType="ExpressionOfConcernIn" or @RefType="ErratumIn"]'))
    list(pmid=pmid,title=title,authors=names,first_author=if(length(names))names[1] else NA_character_,journal=value(paste0(article,if(book)"/Book/BookTitle" else "/Journal/Title")),
      publication_year=year,publication_date=date,electronic_date=electronic,doi=doi,publication_types=types,abstract=abstract,
      pubmed_url=paste0("https://pubmed.ncbi.nlm.nih.gov/",pmid,"/"),doi_url=if(!is.na(doi))paste0("https://doi.org/",utils::URLencode(doi,reserved=TRUE)) else NA_character_,
      notices=xml2::xml_attr(notices,"RefType"),retracted="Retracted Publication" %in% types || any(xml2::xml_attr(notices,"RefType")=="RetractionIn"),
      warnings=if(is.na(title))"Title not supplied in PubMed." else character())
  })
}

get_pubmed_publications <- function(cfg,gene,n=10L,sort="relevance",scope="general",refresh=FALSE) {
  sort<-match.arg(sort,c("relevance","recent"));scope<-match.arg(scope,c("general","medulloblastoma"))
  if(length(n)!=1L || !is.numeric(n) || !is.finite(n) || n!=floor(n) || n<1 || n>30)stop("n must be 1-30")
  identity<-resolve_publication_gene(cfg,gene,refresh)
  query<-publication_query(identity,scope)
  key<-paste0("pubmed_discovery_v1_",object_digest(list(identity$canonical_symbol,query,sort,scope,n)))
  cached<-if(!refresh)cache_get(cfg,key) else NULL
  if(pubmed_fresh(cached)) {
    cached$cache_hit<-TRUE;cached$queried_gene<-gene
    cached$warnings<-unique(c(setdiff(cached$warnings,cached$identity$warnings),identity$warnings))
    identity$raw<-NULL;cached$identity<-identity;return(cached)
  }
  api_sort<-if(sort=="recent")"pub_date" else "relevance"
  raw_search<-pubmed_request(cfg,"esearch.fcgi",list(db="pubmed",term=query,sort=api_sort,retmax=as.integer(n),retmode="json"))
  s<-r2_json(raw_search,"PubMed ESearch");if(!is.null(s$error))stop("PubMed rejected the query")
  s<-s$esearchresult
  if(is.null(s$count) || is.null(s$idlist))stop("Malformed PubMed ESearch result")
  if(!is.null(s$errorlist))stop("PubMed rejected one or more query terms")
  count<-suppressWarnings(as.numeric(s$count));ids<-as.character(s$idlist)
  if(length(count)!=1L || !is.finite(count) || count<0 || any(!grepl("^[0-9]+$",ids)) || anyDuplicated(ids))stop("Malformed PubMed counts or identifiers")
  if(length(ids)!=min(count,n))stop("PubMed returned an unexpected number of identifiers")
  articles<-list();raw_xml<-NULL
  if(length(ids)) {
    raw_xml<-pubmed_request(cfg,"efetch.fcgi",list(db="pubmed",id=paste(ids,collapse=","),retmode="xml"))
    articles<-parse_publication_records(raw_xml);returned<-vapply(articles,`[[`,character(1),"pmid")
    if(anyDuplicated(returned) || !setequal(ids,returned))stop("PubMed search/fetch mismatch; retry retrieval")
    articles<-articles[match(ids,returned)] # Preserve ESearch rank exactly, without local screening or reranking.
    articles<-lapply(seq_along(articles),function(i)c(list(rank=i),articles[[i]]))
  }
  folder<-file.path(cfg$cache_dir,"pubmed_discovery_sources",basename(tempfile("retrieval-")));dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  for(pair in list(c("search.json",raw_search),c("articles.xml",raw_xml),c("gene-search.json",identity$raw$search),c("gene-summary.json",identity$raw$summary)))
    if(length(pair)==2L)writeBin(charToRaw(pair[2]),file.path(folder,pair[1]))
  identity$raw<-NULL
  x<-list(queried_gene=gene,canonical_symbol=identity$canonical_symbol,identity=identity,query=query,query_translation=s$querytranslation,
    scope=scope,sort=sort,api_sort=api_sort,retmax=as.integer(n),retrieved_at=utc_now(),total_hits=count,returned_count=length(articles),
    pmids=ids,articles=articles,cache_hit=FALSE,raw_source_dir=normalizePath(folder,winslash="/"),
    warnings=c(identity$warnings,unlist(s$warninglist,use.names=FALSE),unlist(lapply(articles,`[[`,"warnings"),use.names=FALSE)),
    limitations="Literature discovery, not evidence grading or consensus. Human gene identity does not restrict study species. Symbol/alias ambiguity can remain; biological-context guards may miss papers. No local exclusion of reviews or reranking.")
  cache_put(cfg,key,x);x
}

publication_download_name <- function(x,extension="csv")paste0(x$canonical_symbol,"_PubMed_",if(x$scope=="medulloblastoma")"MB_" else "","top",x$retmax,"_",x$sort,".",extension)

publication_download_rows <- function(x) {
  fields<-c("rank","pmid","title","authors","first_author","journal","publication_year","publication_date","electronic_date","doi","publication_types","abstract","pubmed_url","doi_url","notices","retracted")
  if(length(x$articles))d<-do.call(rbind,lapply(x$articles,function(a)as.data.frame(setNames(lapply(fields,function(f)paste(a[[f]],collapse="; ")),fields),stringsAsFactors=FALSE))) else d<-as.data.frame(setNames(rep(list(character()),length(fields)),fields))
  for(f in c("queried_gene","canonical_symbol","query","scope","sort","api_sort","retrieved_at","total_hits","returned_count"))d[[f]]<-rep(x[[f]],nrow(d))
  d$official_gene_name<-rep(x$identity$official_name,nrow(d));d$aliases_used<-rep(paste(x$identity$aliases_used,collapse="; "),nrow(d));d
}

write_publication_download <- function(x,file,format="csv") {
  if(format=="csv")data.table::fwrite(publication_download_rows(x),file,na="NA") else if(format=="json")
    writeLines(jsonlite::toJSON(x,auto_unbox=TRUE,pretty=TRUE,null="null",na="null"),file,useBytes=TRUE) else stop("Unknown publication format")
  invisible(file)
}
