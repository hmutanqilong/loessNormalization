# functions.
# ---- define batches ---- #
# Batches first defined as time intervals. Five batches defined. For each compound, if a t test shows a non-significant 
# result between neighbor batches, then merge these two batches.
DefineBatches = function(e,f,p,
                         time.intval.batch = TRUE, time = 'TimeStamp',number_of_batches=4,
                         batch = NULL,
                         QC.index = p$QC.index,
                         auto.merge.batch = T, p.adjust.method = 'none'){
  # first define 5 batches according to time interval. We define 5 batches.
  if(time.intval.batch){
    time.diff = diff(p[[time]])
    batch.index = which(time.diff%in%sort(time.diff,decreasing=T)[1:(number_of_batches-1)])
    batch = rep(LETTERS[1:number_of_batches],times = diff(c(0,batch.index,nrow(p))))
  }
  # then if auto.merge.batch is true, for each compound, we use t test to merge two neighbor batches. Rule: if p > 0.05, we merge.
  # this makes differnt compounds have different batch.
  batches = matrix(NA,nrow=nrow(p),ncol=nrow(f))
  if(auto.merge.batch){
    for(i in 1:nrow(f)){
      test=pairwise.wilcox.test(as.numeric(e[i,QC.index]), batch[QC.index],p.adjust.method = p.adjust.method)
      test.p = diag(test$p.value)
      batch.opt = LETTERS[1:(length(which(test.p<0.05))+1)]
      batches[,i] = rep(batch.opt,times = diff(c(0,batch.index[which(test.p<0.05)],nrow(p))))
    }
  }else{
    batches = matrix(rep(batch,ncol(batches)),ncol = ncol(batches))
  }
  return(t(batches))
}
  
get_loess_para = function(x,y,loess.span.limit = 0.5){ # use leave-one-out to select the best span parameter.
  j = 0
  error = rep(0, length(c(seq(loess.span.limit,1.5,0.1),1.75,2,2.25,2.5)))
  if(class(x)=='numeric'|class(x)=='integer'){
    for(par in c(seq(loess.span.limit,1.5,0.1),1.75,2,2.25,2.5)){
      j = j+1
      for(i in 2:(length(y)-1)){ # if i from 1 or to length(y), the prediction would be NA
        o = loess(y[-i]~x[-i],span = par)
        if(sum(is.na(o))){
          error[j] = Inf
        }else{
          err = abs(predict(o,newdata = x[i]) - y[i])
          error[j] = sum(error[j],err,na.rm = T)
        }
      }
    }
  }else if(class(x)=='data.frame'){
    d = data.frame(y,x)
    for(par in c(seq(loess.span.limit,1.5,0.1),1.75,2,2.25,2.5)){
      j = j+1
      for(i in 2:(length(y)-1)){ # if i from 1 or to length(y), the prediction would be NA
        o = loess(y~.,d[-i,],span=par)
        # o = loess(y[-i]~x[-i,],span = par)
        if(sum(is.na(o))){
          error[j] = Inf
        }else{
          err = abs(predict(o,newdata = d[i,colnames(x)]) - y[i])
          error[j] = sum(error[j],err,na.rm = T)
        }
      }
    }
  }else{
    stop("Error: The class of x should be either numeric or matrix")
  }
  
  
  return(c(seq(loess.span.limit,1.5,0.1),1.75,2,2.25,2.5)[which.min(error)])
}
remove_outlier = function(v){ # sometimes, the outlier of QC would be disaster when fitting loess. So we have to remove them.
  out = boxplot(v,plot=F)$out
  return(list(value = v[!v%in%out],index = which(v%in%out)))
}
shiftData = function(ori,norm){
  # to make sure that after normalization the data still at the same scale as raw data.
  ori.min = apply(ori,1,min,na.rm=T)
  norm.min = apply(norm,1,min,na.rm=T)
  return(norm - c(norm.min - ori.min))
}
loessNormalization = function(e,f,p,
         batch = DefineBatches(e,f,p),
         QC.index = p$QC.index,
         time = "TimeStamp",
         span.para = 'auto',loess.span.limit=0.5,
         divide = F){
  # batch is a matrix. columns are samples and rows are compounds.
  # batch = define_batch(e=e,p=p,f=f)
  # for each compound, calculate the loess line.
  library(parallel)
  cl = makeCluster(10)
  e = data.matrix(e)
  norms = parSapply(cl, X = 1:nrow(e), function(i,e,f,p,QC.index,batch,time,remove_outlier,span.para,get_loess_para,
                                                loess.span.limit){
    models = by(data.frame(v=e[i,QC.index],t=p[[time]][QC.index]),
                batch[i,QC.index],function(x){
                  # x = data.frame(v=e[i,QC.index],t=p[[time]][QC.index])[batch[i,QC.index]=="B",]
                  if(length(remove_outlier(x$v)[[2]])>0){# if outlier exists.
                    span = ifelse(span.para=='auto',
                                  get_loess_para(x=x$t[-remove_outlier(x$v)[[2]]],y=remove_outlier(x$v)[[1]],
                                                 loess.span.limit = loess.span.limit),span.para) # find a proper span.
                  }else{
                    span = ifelse(span.para=='auto',
                                  get_loess_para(x=x$t,y=x$v,
                                                 loess.span.limit = loess.span.limit),span.para) # find a proper span.
                    
                  }
                  if(length(remove_outlier(x$v)[[2]])>0){
                    loess(v~t,data=x[-remove_outlier(x$v)[[2]],],span=span)
                  }else{
                    loess(v~t,data=x,span=span)
                  }
                })
    # predict using the models.
    norm = mapply(function(u,v){
      o = tryCatch({
        predict(u,newdata = v)
      },
      error = function(e){
        print(e)
        v
      })
    },models,by(p[[time]],batch[i,],function(x){x}))
    norm = unlist(norm)
    
    # replace NA with the closest value.
    if(length(which(is.na(norm)))>0){
      for(j in which(is.na(norm))){
        time_notNA = p[[time]][-which(is.na(norm))]
        closest_time = time_notNA[which.min(abs(time_notNA - p[[time]][j]))]
        norm[j] = norm[which(p[[time]]==closest_time)]
      }
    }
    return(norm)
  },e,f,p,QC.index,batch,time,remove_outlier,span.para,get_loess_para,loess.span.limit)
  norms = t(norms)
  e_norm = matrix(NA,nrow=nrow(f),ncol=nrow(p))
  if(divide){
    for(i in 1:nrow(f)){
      e_norm[i,] = e[i,] / (norms[i,] / median(e[i,]))
    }
  }else{
    for(i in 1:nrow(f)){
      e_norm[i,] = e[i,] - (norms[i,] - median(e[i,]))
    }
  }
  
  rownames(e_norm) = rownames(e)
  colnames(e_norm) = colnames(e)
  stopCluster(cl)
  return(list(e = shiftData(ori=e,norm = e_norm),f=f,p=p,normalize_line = norms))
}
  
theme_Publication <- function(base_size=14, base_family="helvetica",title.size = 2.4) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(title.size), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}
scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
} 
  
