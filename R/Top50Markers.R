Top50Markers=function(mrkrs){

cm=mrkrs
cm$cluster=as.character(cm$cluster)

#  cut to topNum for each cluster
TopNum=50

nclusters=length(unique(cm$cluster))
MperC=as.vector(table(cm$cluster))   #  Markers per Cluster
MperC=ifelse(MperC >TopNum,TopNum,MperC)

#  Accumulate the genes (features) 
genesFile=NULL
for (i in 1:(nclusters-1))
{
  a=cm[1:MperC[i],] # TopNum genes for cluster i
  c=a$cluster     #  capture cluster identifier
  genesFile=rbind.data.frame(genesFile,a,stringsAsFactors = F)
  #  remove data for cluster just processed
  cm=cm[! (cm$cluster %in% c),]
}
#genesFile$cluster=as.numeric(genesFile$cluster)
genesFile$cluster=as.character(genesFile$cluster)
genesFile    
}
