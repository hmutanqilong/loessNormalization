# Rima's project normalization.
# date 4/26/2017
source("utils.R")
# ---- read data ---- #
pacman::p_load(readxl, data.table,ggplot2)
# positive mode or negative mode
POSorNEG = "POS"
data = read_excel("MX211830_Rima_Lipidomics_ReProcessed.xlsx",sheet = ifelse(POSorNEG=="POS",2,3))
f = data.table(data[3:nrow(data),1]); colnames(f) = "compound"
p = data.table(t(data[c(1,2),])); colnames(p) = as.character(p[1,]); p = p[-1]; p$TimeStamp = as.numeric(p$TimeStamp);
p$index = 1:nrow(p)
order = order(p$TimeStamp)
p = p[order,]
p[,QC.index:=grepl("QC",SampleLabel)]
# define batches according to time.
timeInterval = which(diff(p[,TimeStamp])%in%sort(diff(p[,TimeStamp]),decreasing = T)[1:3])
batch = rep(LETTERS[1:4],diff(c(0,timeInterval,nrow(p))))
p[,batches:=batch]
e = as.matrix(data[3:nrow(data),2:ncol(data)]); e=apply(e,2,as.numeric)
e = e[,order]
batches=DefineBatches(e,f,p)

# loess normalization.
loess = loessNormalization(e,f,p)
e_loess = loess$e

# save normalized data.
fwrite(data.table(e_loess[,order(p$index)]),"loessNormalized.csv")

# RSD.
QC_raw = apply(e,1,function(x,QC.index=p[,QC.index]){
  sd(x[QC.index])/mean(x[QC.index])
})
QC_loess = apply(e_loess,1,function(x,QC.index=p[,QC.index]){
  sd(x[QC.index])/mean(x[QC.index])
})
sample_raw = apply(e,1,function(x,QC.index=p[,QC.index]){
  sd(x[!QC.index])/mean(x[!QC.index])
})
sample_loess = apply(e_loess,1,function(x,QC.index=p[,QC.index]){
  sd(x[!QC.index])/mean(x[!QC.index])
})
fwrite(data.table(QC_raw, QC_loess, sample_raw, sample_loess),"performance.csv")
# PCA performance.
pca = prcomp(t(e),center = T,scale. = T)
pca.d = data.frame(pca$x[,c(1,2)],
                   batch=p$batches,
                   InjectionOrder = 1:nrow(p),
                   type = ifelse(p[,QC.index],"QC"," Sample"),
                   label = 1:nrow(p))
variance = pca$sdev^2/sum(pca$sdev^2)
pca_score = ggplot(pca.d, aes(x=PC1, y=PC2)) + 
  geom_point(aes(colour=batch,size=type,alpha=type))+
  ggtitle(label = paste("raw ",POSorNEG))+
  scale_alpha_discrete(range = c(0.3, 1))+
  # geom_text(aes(label=label), size=3)+
  labs(x = paste0("PC1: ",signif(variance[1]*100,3),"%"), y = paste0("PC2: ",signif(variance[2]*100,3),"%"))
pca_score +scale_colour_Publication()+ theme_Publication()

pca = prcomp(t(e_loess),center = T,scale. = T)
pca.d = data.frame(pca$x[,c(1,2)],
                   batch=p$batches,
                   InjectionOrder = 1:nrow(p),
                   type = ifelse(p[,QC.index],"QC"," Sample"),
                   label = 1:nrow(p))
variance = pca$sdev^2/sum(pca$sdev^2)
pca_score = ggplot(pca.d, aes(x=PC1, y=PC2)) + 
  geom_point(aes(colour=batch,size=type,alpha=type))+
  ggtitle(label = paste("loess Normalization ",POSorNEG))+
  scale_alpha_discrete(range = c(0.3, 1))+
  # geom_text(aes(label=label), size=3)+
  labs(x = paste0("PC1: ",signif(variance[1]*100,3),"%"), y = paste0("PC2: ",signif(variance[2]*100,3),"%"))
pca_score +scale_colour_Publication()+ theme_Publication()




























