# Independent implementation informed by Kaessmann Lab's CC BY 4.0 explorer.
# Source: https://gitlab.com/kaessmannlab/shiny-mammalian-cerebellum
# Publication: https://doi.org/10.1038/s41586-023-06884-x (Sepp, Leiss et al.)
# See scripts/prepare-cerebellum.R and each asset manifest for provenance.
CEREBELLUM_ANNOTATIONS<-c('Broad lineage'='broad_lineage','Cell type'='cell_type',
  'Cell subtype'='subtype','Cell state'='dev_state')
.cerebellum_cache<-new.env(parent=emptyenv())

cerebellum_dir<-function(cfg) cfg$cerebellum_data_dir %||%
  Sys.getenv('CEREBELLUM_DATA_DIR',file.path(cfg$app_root %||% '.', 'data','cerebellum','processed'))

cerebellum_ids<-function(ids,label) {
  if(is.null(ids)||anyNA(ids)||any(!nzchar(ids))||anyDuplicated(ids))
    stop(label,' identifiers must be present and unique')
  invisible(ids)
}

cerebellum_align<-function(cells,expression_ids) {
  cerebellum_ids(cells$cell_id,'Embedding cell');cerebellum_ids(expression_ids,'Exonic cell')
  map<-match(expression_ids,cells$cell_id)
  if(anyNA(map))stop('Exonic cell IDs missing from published embedding')
  if(!identical(cells$cell_id[map],expression_ids))stop('Cell ID alignment failed')
  map
}

cerebellum_orthology<-function(table) {
  require_columns(table,c('human_id','human_symbol','mouse_id','mouse_symbol','homology_type'))
  table<-unique(table)
  valid<-!is.na(table$mouse_id)&nzchar(table$mouse_id)&table$homology_type=='ortholog_one2one'
  pairs<-unique(table[which(valid),c('human_id','mouse_id')])
  bad_h<-unique(pairs$human_id[duplicated(pairs$human_id)|duplicated(pairs$human_id,fromLast=TRUE)])
  bad_m<-unique(pairs$mouse_id[duplicated(pairs$mouse_id)|duplicated(pairs$mouse_id,fromLast=TRUE)])
  table[which(valid & !table$human_id %in% bad_h & !table$mouse_id %in% bad_m),,drop=FALSE]
}

cerebellum_assets<-function(cfg) {
  folder<-normalizePath(cerebellum_dir(cfg),winslash='/',mustWork=TRUE)
  path<-file.path(folder,'manifest.rds')
  if(!file.exists(path))stop('Cerebellum assets have not been prepared')
  key<-paste(folder,file.info(path)$mtime,file.info(path)$size)
  if(exists(key,.cerebellum_cache,inherits=FALSE))return(get(key,.cerebellum_cache))
  manifest<-readRDS(path)
  if(!identical(manifest$schema,1L)||!identical(manifest$quantification,'exonic umi'))stop('Unsupported cerebellum assets')
  species<-lapply(c('human','mouse'),function(s){
    x<-readRDS(file.path(folder,paste0(s,'-index.rds')))
    cerebellum_ids(x$cells$cell_id,'Cached cell')
    cerebellum_ids(x$genes$gene_id,'Cached gene')
    x$binary<-file.path(folder,paste0(s,'-counts.bin'))
    if(!file.exists(x$binary)||file.info(x$binary)$size!=x$bytes)stop('Incomplete cerebellum expression file')
    x
  });names(species)<-c('human','mouse')
  x<-list(species=species,genes=readRDS(file.path(folder,'gene-map.rds')),manifest=manifest)
  # Only one prepared atlas is retained in this process.
  rm(list=ls(.cerebellum_cache),envir=.cerebellum_cache)
  assign(key,x,.cerebellum_cache);x
}

cerebellum_resolve<-function(query,gene_map,human_dataset_ids=NULL) {
  symbol<-toupper(trimws(query));query_id<-sub('\\.[0-9]+$','',symbol)
  hit<-gene_map$human[gene_map$human$human_id==query_id |
    (!is.na(gene_map$human$human_symbol)&toupper(gene_map$human$human_symbol)==symbol),,drop=FALSE]
  ids<-unique(hit$human_id)
  # Resolve duplicate modern symbols only when one exact ID occurs in the
  # published human matrix. Otherwise keep the mapping explicitly ambiguous.
  if(length(ids)>1L&&!is.null(human_dataset_ids)) {
    in_data<-intersect(ids,human_dataset_ids)
    if(length(in_data)==1L){ids<-in_data;hit<-hit[hit$human_id==ids,,drop=FALSE]}
  }
  if(length(ids)!=1L)return(list(human=NULL,mouse=NULL,status='No developing cerebellum snRNA-seq data were found for this gene.'))
  human<-list(id=ids[[1]],symbol=hit$human_symbol[1])
  if(is.na(human$symbol)||!nzchar(human$symbol))human$symbol<-human$id
  pair<-gene_map$pairs[gene_map$pairs$human_id==human$id,,drop=FALSE]
  mouse<-if(nrow(pair)==1L)list(id=pair$mouse_id[[1]],symbol=pair$mouse_symbol[[1]]) else NULL
  if(!is.null(mouse)&&(is.na(mouse$symbol)||!nzchar(mouse$symbol)))mouse$symbol<-mouse$id
  list(human=human,mouse=mouse,status=if(is.null(mouse))'No unambiguous 1:1 mouse orthologue is available.' else 'Resolved')
}

cerebellum_read_gene<-function(asset,gene) {
  if(is.null(gene))return(list(status='No unambiguous 1:1 mouse orthologue is available.'))
  row<-match(gene$id,asset$genes$gene_id)
  if(is.na(row))return(list(status='Gene absent from this species\' exonic dataset.',gene=gene))
  idx<-asset$genes[row,];values<-rep(NA_real_,nrow(asset$cells))
  values[asset$measured_cells]<-0
  con<-file(asset$binary,'rb');on.exit(close(con))
  seek(con,where=idx$offset,origin='start')
  cells<-readBin(con,'integer',n=idx$n,size=4,endian='little')
  counts<-readBin(con,'integer',n=idx$n,size=4,endian='little')
  if(length(cells)!=idx$n||length(counts)!=idx$n||any(cells<1L|cells>length(values))||
    any(counts<=0L)||anyDuplicated(cells)||anyNA(counts))stop('Invalid sparse gene record')
  values[cells]<-counts
  list(status='complete',gene=gene,counts=values,n_measured=sum(!is.na(values)),n_expressing=sum(values>0,na.rm=TRUE))
}

get_cerebellum_profile<-function(cfg,gene) {
  assets<-cerebellum_assets(cfg);identity<-cerebellum_resolve(gene,assets$genes,assets$species$human$genes$gene_id)
  species<-lapply(c('human','mouse'),function(s){
    if(is.null(identity$human))return(list(status=identity$status))
    tryCatch(cerebellum_read_gene(assets$species[[s]],identity[[s]]),error=function(e){
      message('[Cerebellum ',s,'] ',conditionMessage(e));list(status='Expression data are currently unavailable for this species.')})
  });names(species)<-c('human','mouse')
  list(query=gene,identity=identity,species=species,manifest=assets$manifest)
}

cerebellum_colors<-function(assets,field) {
  if(!field %in% unname(CEREBELLUM_ANNOTATIONS))stop('Unsupported annotation')
  labels<-sort(unique(unlist(lapply(assets$species,function(x)x$cells[[field]]))))
  colors<-setNames(grDevices::hcl.colors(length(labels),'Dark 3'),labels)
  colors['Unclassified']<-'#a7adb5';colors
}

cerebellum_plot<-function(asset,field='broad_lineage',expression=NULL,colors=NULL) {
  cells<-asset$cells
  # Static, deterministic stratified overview sample; all positive cells are added
  # for expression plots. Exact published coordinates and full data remain on disk.
  pick<-asset$display_cells
  p<-plotly::plot_ly()
  if(is.null(expression)) {
    for(label in sort(unique(cells[[field]][pick]))) {
      i<-pick[cells[[field]][pick]==label]
      p<-plotly::add_trace(p,x=cells$UMAP1[i],y=cells$UMAP2[i],type='scattergl',mode='markers',
        name=label,text=label,hoverinfo='text',marker=list(color=unname(colors[label]),size=2,opacity=.7))
    }
  } else {
    if(!identical(expression$status,'complete'))stop(expression$status)
    values<-expression$counts
    for(kind in c('unavailable','zero')) {
      i<-pick[if(kind=='zero')!is.na(values[pick])&values[pick]==0 else is.na(values[pick])]
      if(length(i))p<-plotly::add_trace(p,x=cells$UMAP1[i],y=cells$UMAP2[i],type='scattergl',mode='markers',
        name=if(kind=='zero')'0 exonic UMI' else 'No exonic measurement',hoverinfo='name',
        marker=list(color=if(kind=='zero')'#c5c9cf' else '#eceef1',size=2,opacity=.6))
    }
    i<-which(!is.na(values)&values>0);i<-i[order(values[i])]
    if(length(i)) {
      maximum<-max(values[i]);ticks<-unique(round(expm1(seq(0,log1p(maximum),length.out=5))))
      p<-plotly::add_trace(p,x=cells$UMAP1[i],y=cells$UMAP2[i],type='scattergl',mode='markers',
        name='Expressing nuclei',showlegend=FALSE,text=paste(values[i],'exonic UMI'),hoverinfo='text',
        marker=list(size=3,opacity=.8,color=log1p(values[i]),cmin=0,cmax=log1p(maximum),
          colorscale=list(c(0,'#EDBDFE'),c(.5,'#BB36EA'),c(1,'#630C65')),showscale=TRUE,
          colorbar=list(title=list(text='Exonic UMI'),tickvals=log1p(ticks),ticktext=as.character(ticks),len=.7)))
    }
  }
  p<-plotly::layout(p,xaxis=list(title='UMAP 1',showgrid=FALSE,zeroline=FALSE),
    yaxis=list(title='UMAP 2',showgrid=FALSE,zeroline=FALSE,scaleanchor='x',scaleratio=1),
    margin=list(l=45,r=65,b=110,t=10),legend=list(orientation='h',font=list(size=13),
      itemsizing='constant',itemwidth=38,y=-.3,yanchor='top'),
    dragmode='pan',paper_bgcolor='white',plot_bgcolor='white')
  plotly::config(p,displaylogo=FALSE,scrollZoom=TRUE,modeBarButtonsToRemove=c('select2d','lasso2d'))
}

cerebellum_export_plot<-function(asset,expression,title,field='broad_lineage',colors=NULL) {
  d<-asset$cells
  if(is.null(expression)) {
    d$value<-d[[field]]
    p<-ggplot2::ggplot(d,ggplot2::aes(UMAP1,UMAP2,color=value))+ggplot2::geom_point(size=.15)+
      ggplot2::scale_color_manual(values=colors,name=names(CEREBELLUM_ANNOTATIONS)[match(field,CEREBELLUM_ANNOTATIONS)])+
      ggplot2::guides(color=ggplot2::guide_legend(override.aes=list(size=4,alpha=1)))
  } else {
    d$value<-expression$counts;d<-d[order(d$value,na.last=FALSE),]
    p<-ggplot2::ggplot(d,ggplot2::aes(UMAP1,UMAP2))+
      ggplot2::geom_point(data=d[is.na(d$value),],color='#eceef1',size=.15)+
      ggplot2::geom_point(data=d[!is.na(d$value)&d$value==0,],color='#c5c9cf',size=.15)+
      ggplot2::geom_point(data=d[!is.na(d$value)&d$value>0,],ggplot2::aes(color=value),size=.2)+
      ggplot2::scale_color_gradientn(colors=c('#EDBDFE','#BB36EA','#630C65'),transform='log1p',name='Exonic UMI')
  }
  p+ggplot2::coord_equal()+ggplot2::theme_classic(base_size=12)+ggplot2::labs(title=title,
    caption='Sepp, Leiss et al., Nature (2023). Published UMAP; exonic counts. CC BY 4.0.')
}
