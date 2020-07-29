library(gsDesign)


################# STATISTIC TESTS (2-TAILED) REAL ################# 



######   FREQUENTIST TEST  ######

frequentist_test <- function (visit_Start,visit_End,mean_Start,mean_End,var_Start,var_End,alpha,confidence,assessed_uplift){
  
  #print(c(visit_End,mean_End,var_End,assessed_uplift))
  
  # 2-tailed Welch's t-test
  SE <- sqrt((var_Start/visit_Start) + (var_End/visit_End))
  df <- (((var_Start/visit_Start) + (var_End/visit_End))^2)/((((var_Start/visit_Start)^2)/(visit_Start-1)) + ((var_End/visit_End)^2/(visit_End-1)))
  t <- (mean_End - mean_Start)/SE
  
  #print(c(SE,df,t))
  "if (t>0){
    pvalfreq <- 1-pt(t,df)
  }else{
    pvalfreq <- pt(t,df)
  }
  significant <- pvalfreq<=alpha
  assessed_output <- ifelse(significant,sign(t),0) #significant mean difference ?"
  
  pval_freq_sup <- (1-pt(t,df))
  pval_freq_inf <- pt(t,df)
  
  #print(c(pval_freq_sup,pval_freq_inf))
  if (pval_freq_sup<=alpha){
    assessed_output <- 1
    pvalfreq <- pval_freq_sup
  }else if (pval_freq_inf<=alpha){
    assessed_output <- -1
    pvalfreq <- pval_freq_inf
  }else{
    assessed_output <- 0
    pvalfreq <- 1
  }
  
  output <- data.frame(confidence,pvalfreq,assessed_output,assessed_uplift)
  names(output) <- c("confidence","parameters","assessed_output","assessed_uplift")
  return(output)
  
}





######   BAYESIAN TEST  ######

bayesian_test <- function(conv_Start,conv_End,fail_Start,fail_End,confidence,assessed_uplift){
  
  #Beta function calculation  
  endbestthanstart <- 0
  for (i in 0:conv_End) {
    endbestthanstart <- endbestthanstart + exp(lbeta(conv_Start+1+i,fail_End+1+fail_Start+1) - log(fail_End+1+i) - lbeta(1+i,fail_End+1) - lbeta(conv_Start+1,fail_Start+1))
  }
  
  startbestthanend <- 0
  for (i in 0:conv_Start) {
    startbestthanend <- startbestthanend + exp(lbeta(conv_End+1+i,fail_End+1+fail_Start+1) - log(fail_Start+1+i) - lbeta(1+i,fail_Start+1) - lbeta(conv_End+1,fail_End+1))
  }
  
  #print(c(endbestthanstart,startbestthanend))
  
  # Significant ?
  #if((endbestthanstart>=confidence) || (endbestthanstart<=1-confidence)){
  
  if(startbestthanend>=confidence){
    #print("Start")
    output <- data.frame(confidence,startbestthanend,-1,assessed_uplift)
    names(output) <- c("confidence","parameters","assessed_output","assessed_uplift")
  }else if(endbestthanstart>=confidence){
    #print("End")
    output <- data.frame(confidence,endbestthanstart,1,assessed_uplift)
    names(output) <- c("confidence","parameters","assessed_output","assessed_uplift")
  }else{
    #print("No")
    output <- data.frame(confidence,endbestthanstart,0,assessed_uplift)
    names(output) <- c("confidence","parameters","assessed_output","assessed_uplift")
  }
  #.print(c(conv_Start,conv_End,fail_Start,fail_End,confidence,assessed_uplift,endbestthanstart,startbestthanend))
  return(output)
}

######   CHI-2 TEST  ######


chi2_test <- function(visit_Start,visit_End,conv_Start,conv_End,fail_Start,fail_End,mean_Start,mean_End,tot_vis,tot_conv,tot_fail,alpha,assessed_uplift){
  
  confidence<-1-alpha
  
  #print(visit_Start,visit_End,conv_Start,conv_End,fail_Start,fail_End,mean_Start,mean_End,tot_vis,tot_conv,tot_fail)
  
  EA1<-visit_Start*(tot_conv/tot_vis)
  A1<-(conv_Start-EA1)^2/EA1
  
  #print(c(conv_Start,EA1))
  
  EA2<-visit_Start*(tot_fail/tot_vis)
  A2<-(fail_Start-EA2)^2/EA2
  
  EB1<-visit_End*(tot_conv/tot_vis)
  B1<-(conv_End-EB1)^2/EB1
  
  EB2<-visit_Start*(tot_fail/tot_vis)
  B2<-(fail_End-EB2)^2/EB2
  
  #print(c(A1,A2,B1,B2))
  
  chi2 = A1+A2+B1+B2
  df = 1
  
  pvalchi2<- pchisq(chi2, df=1, lower.tail=FALSE)
  
  #print(pvalchi2)
  
  significant <- pvalchi2<=alpha
  assessed_output <- ifelse(significant,sign(mean_End-mean_Start),0) #significant mean difference ?
  
  output <- data.frame(confidence,pvalchi2,assessed_output,assessed_uplift)
  names(output) <- c("confidence","parameters","assessed_output","assessed_uplift") 
  
  return(output)
  
}


######  BOUNDS TEST  ######


bounds_test <- function (visit_Start,visit_End,mean_Start,mean_End,var_Start,var_End,day_nb,alpha,confidence,assessed_uplift,Bounds,N_max_bounds,N_days_max) {
  
  
  #if (bounds_finish==TRUE){
  # return("Done !") #WHAT TO RETURN (LISTE VIDE ????)
  #}
  if (is.null(Bounds)){
    Bounds <- gsDesign(test.type=3, alpha = alpha, beta = 0.2, k=N_max_bounds, sfu=sfPower, sfupar=3, sfl=sfPower, sflpar=2)
  }
  #plot(Bounds)
  
  Upper<-Bounds$upper$bound
  #print(Upper)
  Lower<-Bounds$lower$bound
  
  bounds_day <- ceiling((day_nb*N_max_bounds)/N_days_max)
  
  # 2-tailed Welch's t-test
  SE <- sqrt((var_Start/visit_Start) + (var_End/visit_End))
  #df <- (((var_Start/visit_Start) + (var_End/visit_End))^2)/((((var_Start/visit_Start)^2)/(visit_Start-1)) + ((var_End/visit_End)^2/(visit_End-1)))
  t <- (mean_End - mean_Start)/SE
  #pvalfreq <- (1 - pt(abs(t), df)) * 2
  if (t>=Upper[[bounds_day]]) {
    assessed_output<-1
    output<-data.frame(confidence,t,assessed_output,assessed_uplift)
    names(output)<-c("confidence","parameters","assessed_output","assessed_uplift") 
    #bounds_finish<-TRUE
  }else if (t<=Lower[[bounds_day]]) {
    assessed_output<- -1
    output<-data.frame(confidence,t,assessed_output,assessed_uplift)
    names(output)<-c("confidence","parameters","assessed_output","assessed_uplift") 
    #bounds_finish<-TRUE
  }else{
    assessed_output <- 0
    output<-data.frame(confidence,t,assessed_output,assessed_uplift)
    names(output)<-c("confidence","parameters","assessed_output","assessed_uplift") 
  }
  
  #if (N_week==N_week_max){
  # output <- frequentist_test(visit_Start,visit_End,mean_Start,mean_End,var_Start,var_End,alpha,confidence,assessed_uplift) 
  #}
  return(output)
}


######   FREQUENTIST + BOUNDS TEST  ######


freq_bounds_test <- function (visit_Start,visit_End,mean_Start,mean_End,var_Start,var_End,day_nb,alpha,confidence,assessed_uplift,Bounds,N_max_bounds,N_days_max) {
  
  
  #if (bounds_finish==TRUE){
  # return("Done !") #WHAT TO RETURN (LISTE VIDE ????)
  #}
  if (is.null(Bounds)){
    Bounds <- gsDesign(test.type=3, alpha = alpha, beta = 0.2, k=N_max_bounds, sfu=sfPower, sfupar=3, sfl=sfPower, sflpar=2)
  }
  #plot(Bounds)
  
  Upper<-Bounds$upper$bound
  Lower<-Bounds$lower$bound
  
  bounds_day <- ceiling((day_nb*N_max_bounds)/N_days_max)
  
  # 2-tailed Welch's t-test
  SE <- sqrt((var_Start/visit_Start) + (var_End/visit_End))
  #df <- (((var_Start/visit_Start) + (var_End/visit_End))^2)/((((var_Start/visit_Start)^2)/(visit_Start-1)) + ((var_End/visit_End)^2/(visit_End-1)))
  t <- (mean_End - mean_Start)/SE
  #pvalfreq <- (1 - pt(abs(t), df)) * 2
  if (t>=Upper[[bounds_day]]) {
    assessed_output<-1
    output<-data.frame(confidence,t,assessed_output,assessed_uplift)
    names(output)<-c("confidence","parameters","assessed_output","assessed_uplift") 
    #bounds_finish<-TRUE
  }else if (t<=Lower[[bounds_day]]) {
    assessed_output<- -1
    output<-data.frame(confidence,t,assessed_output,assessed_uplift)
    names(output)<-c("confidence","parameters","assessed_output","assessed_uplift") 
    #bounds_finish<-TRUE
  }else{
    assessed_output <- 0
    output<-data.frame(confidence,t,assessed_output,assessed_uplift)
    names(output)<-c("confidence","parameters","assessed_output","assessed_uplift") 
  }
  
  if (bounds_day==N_max_bounds){
    output <- frequentist_test(visit_Start,visit_End,mean_Start,mean_End,var_Start,var_End,alpha,confidence,assessed_uplift) 
  }
  return(output)
}


######   PERMUTATION TEST (RESAMPLING TEST 1)  ######

## Permutation function ##

perm_func <- function (Total_sample,nS,nE){
  
  n<-nS+nE
  idx_S <- sample(1:n,nS,replace=F) #nS random indexes in [1,n] WITHOUT REPLACEMENT
  idx_E <- setdiff(1:n,idx_S) #the rest of the indexes in [1,n]
  mean_diff <- mean(Total_sample[idx_E])-mean(Total_sample[idx_S])
  return(mean_diff)
}

perm_func2 <- function (Total_sample,nS,nE){
  
  n<-nS+nE
  X<-list(Total_sample)
  Random_sample<-X[[1]][sample(1:n)]
  Random_S <- Random_sample[1:(n/2)]
  Random_E <- Random_sample[((n/2)+1):n]
  mean_diff <- mean(Random_E)-mean(Random_S)
  return(mean_diff)
}

perm_func3 <- function(total_conv,visit_Start,visit_End){
  total_vis0 <- visit_Start+visit_End
  total_vis <- total_vis0
  randS <- 0
  randE <- 0
  while (total_vis>0){
    proba <- total_conv/total_vis
    a <- runif(1,0,1)
    b <- runif(1,0,1)
    if (a<proba){
      randS <- randS+1
      total_conv <- total_conv-1
    }
    proba <- total_conv/total_vis
    if (b<proba){
      randE <- randE+1
      total_conv <- total_conv-1
    }
    total_vis <- total_vis-2
  }
  mean_diff<-(randE/visit_End)-(randS/visit_Start)
  return(mean_diff)
  #rand<-min(sample(1:total_conv, 1),visit_Start)
  #mean_diff<-((total_conv-rand)/visit_End)-(rand/visit_Start)
}



## TEST FUNCTION ##

permutation_test <- function(Start_sample,End_sample,alpha,assessed_uplift) {
  
  visit_Start<-length(Start_sample)
  visit_End<-length(End_sample)
  mean_Start<-mean(Start_sample)
  mean_End<-mean(End_sample)
  delta_mean<-mean_End-mean_Start
  
  R<-100
  Total_sample <- append(Start_sample,End_sample)
  perm_diff <- rep(0,R)
  for (i in 1:R){
    #print(i)
    perm_diff[i]<-perm_func2(Total_sample,visit_Start,visit_End)
  }
  
  "if (delta_mean>0){
    pval_perm <- mean(perm_diff>delta_mean)
  }else{
    pval_perm <- mean(perm_diff<delta_mean)
  }
  significance <- pval_perm<=alpha
  assessed_output <- ifelse(significance,sign(delta_mean),0)"
  
  pval_perm_sup <- mean(perm_diff>=delta_mean)
  pval_perm_inf <- mean(perm_diff<=delta_mean)
  if (pval_perm_sup<=alpha){
    assessed_output <- 1
    pval_perm <- pval_perm_sup
  }else if (pval_perm_inf<=alpha){
    assessed_output <- -1
    pval_perm <- pval_perm_inf
  }else{
    assessed_output <- 0
    pval_perm <- 1
  }
  
  output<-data.frame(1-alpha,pval_perm,assessed_output,assessed_uplift)
  names(output)<-c("confidence","parameters","assessed_output","assessed_uplift") 
  
  return(output)
  
}


permutation_test_real<-function(conv_Start,conv_End,visit_Start,visit_End,alpha,assessed_uplift){
  mean_Start<-conv_Start/visit_Start
  mean_End<-conv_End/visit_End
  delta_mean<-mean_End-mean_Start
  total_conv<-conv_Start+conv_End
  R<-100
  #Total_sample <- append(Start_sample,End_sample)
  perm_diff <- rep(0,R)
  for (i in 1:R){
    #print(i)
    perm_diff[i]<-perm_func3(total_conv,visit_Start,visit_End)
  }
  
  pval_perm <- mean(perm_diff>(delta_mean))
  significance <- pval_perm<=alpha
  assessed_output <- ifelse(significance,1,0)
  
  output<-data.frame(1-alpha,pval_perm,assessed_output,assessed_uplift)
  names(output)<-c("confidence","parameters","assessed_output","assessed_uplift") 
  
  return(output)
}



######   BOOTSTRAP TEST (RESAMPLING TEST 2)  ######

## Bootstrap function ##

boot_func <- function (Total_sample,nS,nE){
  
  n<-nS+nE
  idx_S <- sample(1:n,nS,replace=T) #nS random indexes in [1,n] WITHOUT REPLACEMENT
  idx_E <- sample(1:n,nE,replace=T) #the rest of the indexes in [1,n]
  mean_diff <- mean(Total_sample[idx_E])-mean(Total_sample[idx_S])
  return(mean_diff)
}

boot_func3 <- function(total_conv,visit_Start,visit_End){
  total_vis <- visit_Start+visit_End
  proba <- total_conv/total_vis
  randS <- 0
  randE <- 0
  for (i in 1:visit_Start){
    a <- runif(1,0,1)
    if (a<proba){
      randS <- randS+1
    }
  }
  for (j in 1:visit_End){
    b <- runif(1,0,1)
    if (b<proba){
      randE <- randE+1
    }
  }
  #rand_S<-min(sample(0:total_conv, 1),visit_Start)
  #rand_E<-min(sample(0:total_conv, 1),visit_End)
  mean_diff<-(randE/visit_End)-(randS/visit_Start)
  return(mean_diff)
}




## TEST FUNCTION ##

bootstrap_test <- function(Start_sample,End_sample,alpha,assessed_uplift) {
  
  visit_Start<-length(Start_sample)
  visit_End<-length(End_sample)
  mean_Start<-mean(Start_sample)
  mean_End<-mean(End_sample)
  delta_mean<-mean_End-mean_Start
  
  #print(c(mean_Start,mean_End))
  
  R<-100
  Total_sample <- append(Start_sample,End_sample)
  boot_diff <- rep(0,R)
  for (i in 1:R){
    boot_diff[i]<-boot_func(Total_sample,visit_Start,visit_End)
    #print(boot_diff[i])
  }
  #hist(boot_diff,xlab="Click Through Rate Difference")
  #abline(v=(mean_End-mean_Start))
  
  ##Pvalue version 1##
  
  "if (delta_mean>0){
    pval_boot <- mean(boot_diff>delta_mean)
  }else{
    pval_boot <- mean(boot_diff<delta_mean)
  }
  
  #pval_boot <- mean(boot_diff>abs(delta_mean))+mean(boot_diff<(-abs(delta_mean)))
  significance <- pval_boot<=alpha
  assessed_output <- ifelse(significance,sign(delta_mean),0)
  #print(mean_End-mean_Start)"
  
  
  pval_boot_sup <- mean(boot_diff>=delta_mean)
  pval_boot_inf <- mean(boot_diff<=delta_mean)
  if (pval_boot_sup<=alpha){
    assessed_output <- 1
    pval_boot <- pval_boot_sup
  }else if (pval_boot_inf<=alpha){
    assessed_output <- -1
    pval_boot <- pval_boot_inf
  }else{
    assessed_output <- 0
    pval_boot <- 1
  }
  
  output<-data.frame(1-alpha,pval_boot,assessed_output,assessed_uplift)
  names(output)<-c("confidence","parameters","assessed_output","assessed_uplift") 
  
  return(output)
}

bootstrap_test_real <-function(conv_Start,conv_End,visit_Start,visit_End,alpha,assessed_uplift){
  mean_Start<-conv_Start/visit_Start
  mean_End<-conv_End/visit_End
  delta_mean<-mean_End-mean_Start
  total_conv<-conv_Start+conv_End
  #print(c(mean_Start,mean_End))
  R<-100
  #Total_sample <- append(Start_sample,End_sample)
  boot_diff <- rep(0,R)
  for (i in 1:R){
    boot_diff[i]<-boot_func3(total_conv,visit_Start,visit_End)
  }
  
  pval_boot <- mean(boot_diff>(delta_mean))
  significance <- pval_boot<=alpha
  assessed_output <- ifelse(significance,1,0)
  
  output<-data.frame(1-alpha,pval_boot,assessed_output,assessed_uplift)
  names(output)<-c("confidence","parameters","assessed_output","assessed_uplift") 
  
  return(output)
}

"
time1 <- Sys.time()
print(Sys.time())
print(boot_func3(54000,100000,100000))
time2 <- Sys.time()-time1
print(time2)
print(Sys.time())
out <- bootstrap_test_real(26000,28000,100000,100000,0.2,0.01)
time3 <- Sys.time()-time2
print(Sys.time())


### Test code ###

ConvS <- sample(1:1000,1000)
ConvE <- sample(1:1000,1000)

AO1 <- 0
AO2 <- 0
TA <- 0

time1 <- Sys.time()

for (i in 1:1000){
  print(i)
  moy1 <- 0
  moy2 <- 0
  for (j in 1:100){ 
  out1 <- boot_func(c(rep(1,ConvS[i]+ConvE[i]),rep(0,200-ConvS[i]-ConvE[i])),100,100)
  out2 <- boot_func3(ConvS[i]+ConvE[i],100,100)
  moy1 <- moy1+out1
  moy2 <- moy2+out2
  #print(c(out1,out2))
  }
  #print(c((moy1/100),(moy2/100)))
  #print(c(ConvS[i],ConvE[i]))
  output1 <- bootstrap_test(c(rep(1,ConvS[i]),rep(0,10000-ConvS[i])),c(rep(1,ConvE[i]),rep(0,10000-ConvE[i])),0.2,0.01)
  #output2 <- bootstrap_test_real(ConvS[i],ConvE[i],1000,1000,0.2,0.01)
  AO1 <- AO1 + output1$assessed_output
  AO2 <- AO2 + output2$assessed_output
  TA <- ifelse(ConvE[i]>ConvS[i],TA +1,TA)
  #print(c(i,AO1,output1$parameters,TA))#AO2,output2$parameters,
}

print(Sys.time()-time1)
print(c(AO1,TA)) #AO2
print(ConvS)
print(ConvE)

#print(bootstrap_test_real(10,9,100-10,100-9,0.2,0.01))
#print(bootstrap_test(c(rep(1,10),rep(0,9)),c(rep(1,100-10),rep(0,100-9)),0.2,0.01))
"
