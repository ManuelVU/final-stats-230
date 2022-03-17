# Compute TIP
# @Input: X is a array with dim=3, e.g. X=array(dim=c(n_cow,days,n_Pen))
TIP<-function(X){
  n_pen<-dim(X)[3]
  n_cow<-dim(X)[1]
  TIP<-0
  for(i in 1:n_pen){
    for(j in 1:n_cow){
      x<-X[j,,i][!is.na(X[j,,i])]
      TIP<-TIP+sum(x)
    }
  }
  return(TIP)
}