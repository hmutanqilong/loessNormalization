# Sili Fan UC Davis May 2017
source("utils.R")
# ---- read data ---- #
pacman::p_load(readxl, data.table,ggplot2,gridExtra,ggthemes)
# positive mode or negative mode
POSorNEG = "NEG"
data = read_excel("ADNI_April_2017_reprocessing_scientific_data.xlsx",sheet = ifelse(POSorNEG=="POS",2,3))
f = data.table(data[3:nrow(data),c(1,2,3)]); colnames(f) = c("compound","RT","m.z")
p = data.table(t(data[c(1,2),])); colnames(p) = c("SampleLabel","TimeStamp"); p = p[-c(1,2,3)];
p$index = 1:nrow(p)
order = order(p$TimeStamp)
p$TimeStamp = as.numeric(p$TimeStamp)
p = p[order,]
p[,QC.index:=grepl("QC",SampleLabel)]
# define batches according to time.
timeInterval = which(diff(p[,TimeStamp])%in%sort(diff(p[,TimeStamp]),decreasing = T)[1:3])
batch = rep(LETTERS[1:4],diff(c(0,timeInterval,nrow(p))))
p[,batches:=batch]
e = as.matrix(data[3:nrow(data),4:ncol(data)]); e=apply(e,2,as.numeric)
e = e[,order]
batches=DefineBatches(e,f,p,auto.merge.batch = F)



QC_raw.train  = QC_loess.train = QC_raw.test = QC_loess.test = sample_raw= sample_loess = list()

for(i in 1:10){
  print(paste0("WORKING on ", i, 'th round'))
  
  QC.index = which(p[,QC.index])
  QC.index.train.temp = sample(QC.index,round(length(QC.index)*.8))
  QC.index.test.temp = QC.index[!QC.index%in%QC.index.train.temp]
  QC.index.train = rep(F,nrow(p))
  QC.index.test = rep(F,nrow(p))
  QC.index.train[QC.index.train.temp] = T
  QC.index.test[QC.index.test.temp] = T
  p[,QC.index.train := QC.index.train]
  p[,QC.index.test := QC.index.test]
  
  
  loess = loessNormalization(e,f,p,QC.index = p[,QC.index.train])
  e_loess = loess$e
  
  # # RSD.
  QC_raw.train[[i]] = apply(e,1,function(x,QC.index=p[,QC.index.train]){
    sd(x[QC.index])/mean(x[QC.index])
  })
  QC_loess.train[[i]] = apply(e_loess,1,function(x,QC.index=p[,QC.index.train]){
    sd(x[QC.index])/mean(x[QC.index])
  })
  
  QC_raw.test[[i]] = apply(e,1,function(x,QC.index=p[,QC.index.test]){
    sd(x[QC.index])/mean(x[QC.index])
  })
  QC_loess.test[[i]] = apply(e_loess,1,function(x,QC.index=p[,QC.index.test]){
    sd(x[QC.index])/mean(x[QC.index])
  })
  
  sample_raw[[i]] = apply(e,1,function(x,QC.index=p[,QC.index]){
    sd(x[!QC.index])/mean(x[!QC.index])
  })
  sample_loess[[i]] = apply(e_loess,1,function(x,QC.index=p[,QC.index]){
    sd(x[!QC.index])/mean(x[!QC.index])
  })
}
QC_raw.train. = do.call("cbind",QC_raw.train)
QC_loess.train. = do.call("cbind",QC_loess.train)
QC_raw.test. = do.call("cbind",QC_raw.test)
QC_loess.test. = do.call("cbind",QC_loess.test)
sample_raw. = do.call("cbind",sample_raw)
sample_loess. = do.call("cbind",sample_loess)

QC_raw.train.. = apply(QC_raw.train.,1,mean)
QC_loess.train.. = apply(QC_loess.train.,1,mean)
QC_raw.test.. = apply(QC_raw.test.,1,mean)
QC_loess.test.. = apply(QC_loess.test.,1,mean,na.rm=T)
sample_raw.. = apply(sample_raw.,1,mean)
sample_loess.. = apply(sample_loess.,1,mean)

fwrite(data.table(QC_raw.train.., QC_loess.train..,
                  QC_raw.test..,QC_loess.test..,
                  sample_raw..,sample_loess..),"performance_NEG.csv")

# loess normalization.
loess = loessNormalization(e,f,p)
e_loess = loess$e

# save normalized data.
fwrite(data.table(e_loess[,order(p$index)]),"loessNormalized.csv")

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




























