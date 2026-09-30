# Browser-only dimensions; export functions keep their own publication dimensions.
explorer_plot_height <- function(size="standard") {
  sizes<-c(small="280px",standard="420px",large="560px",pathways="650px")
  unname(sizes[[size]])
}
explorer_web_plot <- function(plot,id,session) {
  if(is.null(plot))return(plot)
  width<-session$clientData[[paste0("output_",id,"_width")]]
  if(is.null(width)||!is.finite(width)||width<100)width<-900
  wrap<-function(x,n) {
    if(!is.character(x)||length(x)!=1L)return(x)
    paste(vapply(strsplit(x,"\n",fixed=TRUE)[[1]],function(line)paste(strwrap(line,width=n),collapse="\n"),character(1)),collapse="\n")
  }
  labels<-plot$labels
  if(is.character(labels$y)&&length(labels$y)==1L)labels$y<-wrap(labels$y,22)
  for(k in intersect(c("title","subtitle","caption","x"),names(labels)))
    labels[[k]]<-wrap(labels[[k]],max(22L,floor((width-90)/if(k %in% c("title","subtitle","x"))11 else 8)))
  plot<-plot+do.call(ggplot2::labs,labels[intersect(c("title","subtitle","caption","x","y"),names(labels))])+ggplot2::theme(
    plot.title=ggplot2::element_text(size=13),plot.subtitle=ggplot2::element_text(size=11),
    plot.caption=ggplot2::element_text(size=9,hjust=0),plot.margin=ggplot2::margin(8,12,8,8),
    plot.title.position="plot",plot.caption.position="plot",axis.title.y=ggplot2::element_text(size=10))
  if(width<850) {
    plot<-plot+ggplot2::theme(legend.position="bottom",legend.box="vertical",legend.title.position="top",
      legend.title=ggplot2::element_text(size=10),legend.text=ggplot2::element_text(size=10),
      legend.margin=ggplot2::margin(2,0,2,0),legend.spacing.y=grid::unit(.1,"cm"))
    if(id!="enrichment_plot")plot<-plot+ggplot2::guides(color=ggplot2::guide_legend(nrow=2),fill=ggplot2::guide_legend(nrow=2))
  }
  if(id=="subtype_plot" && width<850) {
    label_map<-plot$scales$get_scales("x")$labels
    if(is.character(label_map))plot<-suppressMessages(plot+ggplot2::scale_x_discrete(labels=gsub("\n"," ",label_map,fixed=TRUE)))
    plot<-plot+ggplot2::coord_flip()+ggplot2::labs(x="R2 subtype")+ggplot2::theme(axis.text.y=ggplot2::element_text(angle=0,hjust=1),axis.text.x=ggplot2::element_text(angle=0,hjust=.5))
  }
  if(id %in% c("subtype_plot","pfister_plot","subgroup_plot"))plot<-plot+ggplot2::guides(fill="none")
  plot
}
explorer_style <- function() shiny::tags$style(shiny::HTML("
  html {font-size:15px}
  body {line-height:1.45}
  h1 {font-size:1.9rem} h2 {font-size:1.55rem} h3 {font-size:1.3rem} h4 {font-size:1.1rem}
  .nav-link {font-size:.95rem} .card {--bs-card-spacer-y:.85rem;--bs-card-spacer-x:1rem}
  .card p {margin-bottom:.65rem} .shiny-input-container {margin-bottom:.8rem}
  .table {font-size:.92rem} .table>:not(caption)>*>* {padding:.38rem .5rem}
  .shiny-html-output:has(>table) {overflow-x:auto}
  .shiny-plot-output {width:100%;max-width:100%}
  .tab-pane {min-width:0} .alert {padding:.7rem 1rem}
  .nav-tabs {flex-wrap:wrap} .card-header a {overflow-wrap:anywhere}
  .explorer-methods {margin-top:1rem;color:#45545d;font-size:.9rem}
"))
