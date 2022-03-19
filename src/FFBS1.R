# Forward filtering backward sampling  within a pen - consider missing data 
# cow: index of cow to be updated
# s: the day that a cow from the pen withdrawn
# a: external infection rate
# b: with-in pen infection rate
# m0: mean infectious period - 1
# nu: initial probability of infection
# ObservedR, ObservedF: Observed test results for RAMS and faecal respectively
# thetaR, thetaF : sensitivity of the RAMS and faecal test respectively
# @OUtput of function: generated hidden states
FFBS1<-function(a, b,m0, nu, cow, s, X,thetaR, thetaF, ObservedR, ObservedF,id_pen){
  #forward filtering
  I <- apply(X[,-cow, id_pen], 1, function(x) sum(x, na.rm = T))
  s<-length(X[,cow, id_pen])-sum(is.na(X[, cow, id_pen]))
  predictive.prob_X<-matrix(0,nrow = s,ncol = 2) #Compute the one-step ahead modified conditional predictive probabilities
  filtered.prob_X<-matrix(0,nrow = s,ncol = 2) #compute the modified conditional filtered prob
  predictive.prob_X[1,]<-c(nu,1-nu)
  prob1<-predictive.prob_X[1,1] # Xt=1
  #emission matrix:
  f<-matrix(0,nrow = 2,ncol = 4)
  f[1,]<-c(thetaR*thetaF, thetaR*(1-thetaF), (1-thetaR)*thetaF, (1-thetaR)*(1-thetaF))
  f[2,]<-c(0,0,0,1)  # y: (1,1),(1,0),(0,1),(0,0)
  
  if(is.na(ObservedR[1, cow, id_pen])==FALSE){
    if((ObservedR[1, cow, id_pen]==0)&(ObservedF[1, cow, id_pen]==0)){
      no1<-prob1*f[1,4]  # x1=1，y1=(0，0)
      no2<-(1-prob1)*f[2,4] # x1=0,y1=(0,0)
    }else if((ObservedR[1, cow, id_pen]==1)&(ObservedF[1, cow, id_pen]==0)){
      no1<-prob1*f[1,2]    # x1=1，y1=(1,0)
      no2<-(1-prob1)*f[2,2]  # x1=0，y1=(1,0)
    }else if((ObservedR[1, cow, id_pen]==0)&(ObservedF[1, cow, id_pen]==1)){
      no1<-prob1*f[1,3]   # x1=1，y1=(0,1)
      no2<-(1-prob1)*f[2,3]  # x1=0，y1=(0,1)
    }else{
      no1<-prob1*f[1,1]  # x1=1，y1=(1,1)
      no2<-(1-prob1)*f[2,1] # x1=0,y1=(1,1)
    }
    filtered.prob_X[1,]<-c(no1/(no1+no2),no2/(no1+no2))
  }else{
    filtered.prob_X[1,]<-predictive.prob_X[1,]  # observation is empty, f_xi(y_t)=1
  }
  # t = 2 to s
  for(i in 2:s){
    #transition prob:
    p00<-exp(-a-b*I[i-1])
    p11<-1/(m0+1)
    # compute predictive prob
    predictive.prob_X[i,1]<-filtered.prob_X[i-1,1]*(1-p11)+filtered.prob_X[i-1,2]*(1-p00)
    predictive.prob_X[i,2]<-filtered.prob_X[i-1,1]*p11+filtered.prob_X[i-1,2]*p00
    # Compute filtered prob
    if(is.na(ObservedR[i])==FALSE){
      if((ObservedR[1, cow, id_pen]==0)&(ObservedF[1, cow, id_pen]==0)){
        no1<-predictive.prob_X[i,1]*f[1,4]  # x1=1，y1=(0，0)
        no2<-predictive.prob_X[i,2]*f[2,4] # x1=0,y1=(0,0)
      }else if((ObservedR[1, cow, id_pen]==1)&(ObservedF[1, cow, id_pen]==0)){
        no1<-predictive.prob_X[i,1]*f[1,2]    # x1=1，y1=(1,0)
        no2<-predictive.prob_X[i,2]*f[2,2]  # x1=0，y1=(1,0)
      }else if((ObservedR[1, cow, id_pen]==0)&(ObservedF[1, cow, id_pen]==1)){
        no1<-predictive.prob_X[i,1]*f[1,3]   # x1=1，y1=(0,1)
        no2<-predictive.prob_X[i,2]*f[2,3]  # x1=0，y1=(0,1)
      }else{
        no1<-predictive.prob_X[i,1]*f[1,1]  # x1=1，y1=(1,1)
        no2<-predictive.prob_X[i,2]*f[2,1] # x1=0,y1=(1,1)
      }
      filtered.prob_X[i,]<-c(no1/(no1+no2),no2/(no1+no2))
    }else{
      filtered.prob_X[i,]<-predictive.prob_X[i,]/sum(predictive.prob_X[i,])
    }
  }
  ## Backward sampling
  Hstate<-c()
  # compute proposal distribution
  proposal.prob_X<-matrix(0,nrow = s,ncol = 2)
  # Likelihood
  likelihood<-c()
  u<-runif(1)
  if(u<=filtered.prob_X[s,1]){
    Hstate[s]<-1
    likelihood[s]<-filtered.prob_X[s,1]
  }else{
    Hstate[s]<-0
    likelihood[s]<-filtered.prob_X[s,2]}
  for(i in (s-1):1){
    #transition prob 
    p00<-exp(-a-b*I[i])
    p11<-1/(m0+1)
    
    if(Hstate[i+1]==1){    # X_{t+1} = 1
      num1<-filtered.prob_X[i,1]*p11
      num2<-filtered.prob_X[i,2]*(1-p00)
      proposal.prob_X[i,1]<-num1/(num1+num2) # Xt = 1
      proposal.prob_X[i,2]<-1-proposal.prob_X[i,1] # Xt = 0
    }else{  # X_{t+1} = 0
      num1<-filtered.prob_X[i,1]*(1-p11)
      num2<-filtered.prob_X[i,2]*p00
      proposal.prob_X[i,1]<-num1/(num1+num2) # Xt = 1
      proposal.prob_X[i,2]<-1-proposal.prob_X[i,1] # Xt = 0
    }
    u1<-runif(1)
    if(u1<=proposal.prob_X[i,1]){
      Hstate[i]<-1
      likelihood[i]<-proposal.prob_X[i,1]
    }else{
      Hstate[i]<-0
      likelihood[i]<-proposal.prob_X[i,2]}
  }
  result<-list(Hstate,sum(log(likelihood)))
  return(result)
}
