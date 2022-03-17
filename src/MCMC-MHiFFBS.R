### MCMC - MH sampler
MCMC_MHiFFBS<-function(N,burning,ObservedR,ObservedF,n_Pen,n_cow,ini_a,ini_b,ini_m,ini_thetaR,ini_thetaF,ini_mu){
  # Initialization of parameters:
  a<-ini_a
  b<-ini_b
  m<-ini_m
  mu<-ini_mu
  thetaR<-ini_thetaR
  thetaF<-ini_thetaF
  # generate hidden states:
  likelihood_x<-matrix(nrow = n_cow,ncol =n_Pen ) 
  for(i in 1:n_Pen){
    for(j in 1:n_cow){
      X[j,,i]<-sis_sim_hidenstates(dim(R)[3],nc,a, b,m,p_initial)  ##??
      likelihood1[j,i]<-pi_traget()
    }
  }
  #TIP
  tip<-c()
  #MH sampling
  for(i in 1:N){
    for(j in 1:n_Pen){
      for(k in 1:n_cow){
        #(proposal_pinfected, pen_states, id_number, id_pen, 
         #external_inf_rate, within_inf_rate, 
         #infectious_period)
        sampling<-FFBS1()
        X_prop<-sampling[[1]]
        likeli_prop<-sampling[[2]]
        pi_cur<-1
        pi_prop<-2
        # Compute acceptance ratio
        accept<-min(1,(likelihood1[k,j]*pi_prop)/(likeli_prop*pi_cur))
        #accept-reject step
        u1<-runif(1)
        if(u1<accept){
          X[k,,j]<-X_prop
        }else{X[k,,j]<-X[k,,j]}
        
      }
      
      # HMC update  a,b,m
      cur_theta<-hamiltonian_mc(,30,0.03) 
      a<-cur_theta[1]
      b<-cur_theta[2]
      m<-cur_theta[3]
      # update thetaR thetaF
      b_thetaR<-1;c_thetaR<-1  # prior params for thetaR
      b_thetaF<-1;c_thetaF<-1  # prior params for thetaR
      sum_r1<-0;sum_r2<-0 # true positive test of R; false positive test of R
      sum_f1<-0;sum_f2<-0 # true positive test of F; false positive test of F
      for(i in 1:n_Pen){
        for(j in 1:n_cow){
          t<-sum(!is.na(X[j,,i]))
          for(k in 1:t){
            if(X[j,k,i]==1){
              if(ObservedR[j,k,i]==1){
                sum_r1<-sum_r1+1
              }else{sum_r2<-sum_r2+1}
              if(ObservedF[j,k,i]==1){
                sum_f1<-sum_f1+1
              }else{sum_f2<-sum_f2+1}
            }
          }
        }
      }
      thetaR_new<-rbeta(sum_r1+b_thetaR,sum_r2+c_thetaR)
      thetaF_new<-rbeta(sum_f1+b_thetaF,sum_f2+c_thetaF)
    }
    tip[i]<-TIP(X[[i]])
  }
  ACF_TIP<-acf(tip)$acf
  theta<-c(a,b,m,thetaR_new,thetaF_new)
}