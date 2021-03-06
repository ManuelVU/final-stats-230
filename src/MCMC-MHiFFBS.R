### MCMC - MH sampler
MCMC_MHiFFBS<-function(N,burning,ObservedR,ObservedF,n_Pen,n_cow,ini_a,ini_b,ini_m,ini_thetaR,ini_thetaF,ini_mu){
  # load functions
  source(file = "src/simulation-hidden-states.R")
  source(file = "src/metropolis-hastings-hidden-states.R")
  #source(file = "src/hamiltonian-mc.R")
  source(file = "src/TIP.R")
  source(file = "src/gradient-log-joint.R")
  source(file = "src/log-joint-posterior.R")
  # Initialization of parameters:
  a1<-ini_a
  b1<-ini_b
  m1<-ini_m # mean of infected periods 
  mu<-ini_mu
  thetaR<-ini_thetaR
  thetaF<-ini_thetaF
  ini_theta<-c(a1,b1,m1-1,mu,thetaR,thetaF)
  n_day<-dim(ObservedR)[1]
  #Samples
  a_s<-c(); b_s<-c(); m_s<-c()
  thetaR_s<-c();thetaF_s<-c(); mu_s<-c()
  # generate hidden states:
  X<-array(dim = c(n_day,n_cow,n_Pen))
  for(i in 1:n_Pen){
      X[,,i]<-sis_sim_hidenstates(n_day,n_cow,a1, b1,m1-1,mu) 
  }
  #TIP
  tip<-c()
  #MH sampling
  # Acceptance rate:
  accept<-array(dim=c(n_cow,n_Pen,N))
  # Compute number of infected cows
  n_infected<-matrix(nrow = N,ncol = n_day)
  theta1<-c(a1,b1,m1-1,mu,thetaR,thetaF)
  for(i in 1:N){
    for(j in 1:n_Pen){
      for(k in 1:n_cow){
        update<-mh_chain(X,ObservedR,ObservedF,k,j,theta1,log_q_current = NULL, log_pi_current = NULL)
        X[,k,j]<-update$mhchain
        accept[k,j,i]<-update$accepted
      }
    }
    # HMC update  a,b,m
    theta_tilde<-c(log(a1),log(b1),m1-1)
    initial_theta<-matrix(c(1,1,1,1,1,10),nrow = 3,ncol = 2,byrow = TRUE)
    updated_theta<-hamiltonian_mc(u=log_joint_post,grad_u=gradient_log_joint_post,theta_tilde,initial_theta,X,30,0.001) 
    a_s[i]<-a1<-exp(updated_theta[1])
    b_s[i]<-b1<-exp(updated_theta[2])
    m_s[i]<-m1<-updated_theta[3]+1
    # update thetaR thetaF
    b_thetaR<-1;c_thetaR<-1  # prior params for thetaR
    b_thetaF<-1;c_thetaF<-1  # prior params for thetaR
    sum_r1<-0;sum_r2<-0 # true positive test of R; false positive test of R
    sum_f1<-0;sum_f2<-0 # true positive test of F; false positive test of F
    for(i1 in 1:n_Pen){
      for(j in 1:n_cow){
        t<-sum(!is.na(X[,j,i1]))
        for(k in 1:t){
          if(X[k,j,i1]==1){
            if(ObservedR[k,j,i1]==1){
              sum_r1<-sum_r1+1
            }else{sum_r2<-sum_r2+1}
            if(ObservedF[k,j,i1]==1){
              sum_f1<-sum_f1+1
            }else{sum_f2<-sum_f2+1}
          }
        }
      }
    }
    thetaR_s[i]<-thetaR<-rbeta(1,sum_r1+b_thetaR,sum_r2+c_thetaR)
    thetaF_s[i]<-thetaF<-rbeta(1,sum_f1+b_thetaF,sum_f2+c_thetaF)
    # Update mu
    b_v<-1;c_v<-1;
    sum_mu1<-sum(X[1,,])
    sum_mu2<-sum(1-X[1,,])
    mu_s[i]<-mu<-rbeta(1,sum_mu1+b_v,sum_mu2+c_v)
    #number of infected cows
    for(l in 1:n_day){
      n_infected[i,l]<-sum(X[l,,])
    }
    # Compute TIP
    tip[i]<-TIP(X)
  }
  effective<-c((burning+1):N)
  acceptance<-mean(accept[effective])
  ACF_TIP<-acf(tip[effective],lag.max = 30)$acf
  
  theta0<-matrix(nrow = 6,ncol = 5) #a1,b1,m1,thetaR,thetaF,mu
  theta0[1,]<-quantile(a_s[effective])
  theta0[2,]<-quantile(b_s[effective])
  theta0[3,]<-quantile(m_s[effective])
  theta0[4,]<-quantile(thetaR_s[effective])
  theta0[5,]<-quantile(thetaF_s[effective])
  theta0[6,]<-quantile(mu_s[effective])
  rownames(theta0)<-c("alpha","beta","m","thetaR","thetaF","mu")
  n_infected_est<-colMeans(n_infected[effective,])
  CI_lower<-c();CI_upper<-c()
  for(i in 1:n_day){
    CI_lower[i]<-quantile(n_infected[effective,i],0.025,na.rm = TRUE)
    CI_upper[i]<-quantile(n_infected[effective,i],0.975,na.rm = TRUE)
  }
  N_infected<-cbind(n_infected_est,CI_lower,CI_upper)
  result<-list(theta0,acceptance,N_infected,tip,ACF_TIP)
  names(result)<-c("post_theta","acceptance","number of infected cows","TIP","ACF of TIP")
  return(result)
}

output<-MCMC_MHiFFBS(50,20,sis_data$resp$rams,sis_data$resp$fec,20,8,0.02,0.03,3,0.9,0.95,0.3)
save(output, file = "data/output.Rdata")



