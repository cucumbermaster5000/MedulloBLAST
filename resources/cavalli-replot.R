# Standalone: only local CSV/JSON files plus ggplot2 and jsonlite are required.
cavalli_replot <- function(directory='.') {
  for(package in c('ggplot2','jsonlite'))if(!requireNamespace(package,quietly=TRUE))stop('Install the ',package,' package first.')
  m<-jsonlite::read_json(file.path(directory,'metadata.json'),simplifyVector=TRUE)
  read_table<-function(name)utils::read.csv(file.path(directory,paste0(name,'.csv')),check.names=FALSE,stringsAsFactors=FALSE,na.strings='NA',colClasses='character',fileEncoding='UTF-8')
  s<-m$plotting;groups<-as.character(unlist(s$groups));colors<-setNames(as.character(unlist(s$colors)),groups)
  labels<-setNames(as.character(unlist(s$group_labels)),groups)
  if(m$plot=='survival') {
    d<-read_table('curves');for(k in c('time','survival','lower','upper','censored'))d[[k]]<-as.numeric(d[[k]])
    p<-ggplot2::ggplot(d,ggplot2::aes(time,survival,color=group,group=group))+
      ggplot2::geom_step(linewidth=.8)+
      ggplot2::geom_step(ggplot2::aes(y=lower),linetype=3,alpha=.4,na.rm=TRUE)+
      ggplot2::geom_step(ggplot2::aes(y=upper),linetype=3,alpha=.4,na.rm=TRUE)+
      ggplot2::geom_point(data=d[d$censored>0,],shape=3,size=1.8)+
      ggplot2::scale_color_manual(values=colors,breaks=groups,labels=labels)+
      ggplot2::coord_cartesian(ylim=c(0,1))+ggplot2::theme_classic(base_size=12)
  } else {
    d<-read_table('plot_data');d$expression<-as.numeric(d$expression)
    field<-switch(m$plot,subgroups='subgroup',subtypes='plot_subtype',metastasis='metastasis_label',stop('Unknown plot type'))
    d$plot_group<-factor(d[[field]],levels=groups)
    p<-ggplot2::ggplot(d,ggplot2::aes(plot_group,expression,fill=plot_group))
    if(m$plot=='subgroups')p<-p+ggplot2::geom_boxplot(width=.5,outlier.shape=NA,alpha=.35,linewidth=.45)+
      ggplot2::geom_point(position=ggplot2::position_jitter(width=.16,height=0,seed=85217),size=1.1,alpha=.5,shape=16)
    if(m$plot=='subtypes')p<-p+ggplot2::geom_boxplot(width=.6,outlier.shape=NA,alpha=.45)+
      ggplot2::geom_point(position=ggplot2::position_jitter(width=.18,height=0,seed=85217),size=.9,alpha=.5)
    if(m$plot=='metastasis')p<-p+ggplot2::geom_boxplot(outlier.shape=NA,alpha=.4)+
      ggplot2::geom_point(position=ggplot2::position_jitter(width=.15,seed=85217),alpha=.5,size=1)
    p<-p+ggplot2::scale_x_discrete(drop=m$plot!='subtypes',labels=labels)+
      ggplot2::scale_fill_manual(values=colors,na.value='#777777')+ggplot2::guides(fill='none')+ggplot2::theme_classic(base_size=12)
    if(m$plot=='subtypes')p<-p+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1),legend.position='top',plot.caption=ggplot2::element_text(hjust=0,size=9))
    if(m$plot=='subgroups')p<-p+ggplot2::theme(plot.caption=ggplot2::element_text(hjust=0,size=9))
  }
  do.call(ggplot2::labs,s$labels)->lab
  p+lab
}
if(sys.nframe()==0L) {
  arg<-grep('^--file=',commandArgs(),value=TRUE)
  directory<-dirname(normalizePath(sub('^--file=','',arg[1]),mustWork=TRUE))
  args<-commandArgs(trailingOnly=TRUE);out<-if(length(args))args[1]else directory
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  p<-cavalli_replot(directory);m<-jsonlite::read_json(file.path(directory,'metadata.json'),simplifyVector=TRUE)
  ggplot2::ggsave(file.path(out,'reproduced_plot.png'),p,width=m$plotting$width,height=m$plotting$height,dpi=300)
  ggplot2::ggsave(file.path(out,'reproduced_plot.pdf'),p,width=m$plotting$width,height=m$plotting$height,device=grDevices::pdf)
  message('Recreated plot from local package files: ',normalizePath(out))
}
