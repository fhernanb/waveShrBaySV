errorF<-function(y,vec){
  SMem<-matrix(NA,nrow=ncol(vec),ncol=1)
  rownames(SMem)<-names(vec)
  colnames(SMem)<-"RMSE"
  Serr<-NULL
  for(j in 1:ncol(vec)){
    for(i in 1:length(y)){
      Serr[i]<-(vec[i,j]-y[i])^2
    }
    SMem[j,]<-sqrt(sum(Serr)/length(y))
  }
  return(SMem)
}
