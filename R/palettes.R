# Biological identity only; independent of factor order and of UI accent colours.
MB_GROUP_COLORS <- c(wnt="#5266AD",shh="#E23650",group3="#FBDB14",group4="#368450")
CAVALLI_SUBTYPE_COLORS <- c(wnt="#5266AD",shh_beta="#823C34",shh_gamma="#B78CBA",
  shh_alpha="#CE5248",shh_delta="#E8A4A5",group3_beta="#DBB840",group3_gamma="#FCB333",
  group3_alpha="#FBFA49",group4_alpha="#22F84B",group4_gamma="#1C822C",group4_beta="#D9ECD0")
HPA_PROFILE_PALETTE <- c("#EDBDFE","#D698FD","#BB36EA","#7D117F","#630C65")
mb_color_key <- function(labels) {
  x<-tolower(trimws(as.character(labels)))
  for(pair in list(c("\u03b1","alpha"),c("\u03b2","beta"),c("\u03b3","gamma"),c("\u03b4","delta")))x<-gsub(pair[1],pair[2],x,fixed=TRUE)
  x<-gsub("[^a-z0-9]","",x)
  x<-sub("^(medulloblastoma|mb)","",x)
  x<-sub("^g([34])","group\\1",x)
  x<-sub("like$","",x)
  x<-sub("(alpha|beta|gamma|delta)$","_\\1",x)
  x
}
category_colors <- function(labels,cavalli=FALSE) {
  labels<-unique(as.character(labels));keys<-mb_color_key(labels)
  # Stable label-derived hues do not change when categories are absent.
  hues<-vapply(labels,function(s)sum(utf8ToInt(tolower(s))*seq_along(utf8ToInt(tolower(s))))%%360,numeric(1))
  colors<-grDevices::hcl(h=hues,c=48,l=62)
  colors[keys %in% c("","medulloblastoma")]<-"#34465E"
  broad<-match(keys,names(MB_GROUP_COLORS));ok<-!is.na(broad);colors[ok]<-MB_GROUP_COLORS[broad[ok]]
  if(cavalli) {
    keys[keys %in% c("wnt_alpha","wnt_beta")]<-"wnt"
    i<-match(keys,names(CAVALLI_SUBTYPE_COLORS));ok<-!is.na(i);colors[ok]<-CAVALLI_SUBTYPE_COLORS[i[ok]]
  }
  setNames(colors,labels)
}
