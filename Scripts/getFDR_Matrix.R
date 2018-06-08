getFDR_Matrix <- function(p1, p0, alpha=0.1, z=NULL, subset=NULL){
  if(is.null(z)){
    a=0
    for(itr in 1:10){
      a=getFDR_Matrix(p1, p0, alpha, rev(a+0:100/100^itr), subset)
    }
    a
  }else{
    if(!is.null(subset)){
      p1=p1[subset]
      p0=p0[subset]
    }
    p1=p1[!is.na(p1)]
    p0=p0[complete.cases(p0),]
    x=NULL;
    for(i in z){
      x=c(x, min(colSums(p0<i))/nrow(p0)/(sum(p1<i)/length(p1)))
    };
    max(c(0,z[x<alpha]),na.rm=T)
  }
}