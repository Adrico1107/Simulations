source(file="StatsSim.R")


######   BIG TEST FUNCTION  ######


testing_distribution <- function(Start_distrib, End_distrib, methodology_used, alpha , N_week,Bounds,esti_obs,N_week_max) {
  
  #### INITIALISE PARAMETERS ####  
  
  confidence <- 1-alpha
  
  conv_Start <- sum(Start_distrib) #number of conversions on start page
  conv_End <- sum(End_distrib) #number of conversions on end page
  
  visit_Start <- length(Start_distrib) #number of visits on start page
  visit_End <- length(End_distrib)
  
  
  tot_vis<-visit_Start+visit_End 
  
  fail_Start <- visit_Start - conv_Start #number of failures on start page
  fail_End <- visit_Start - conv_End
  
  tot_conv<-conv_Start+conv_End
  tot_fail<-fail_Start+fail_End
  
  mean_Start <- mean(Start_distrib) #click through rate (CTR)
  mean_End <- mean(End_distrib)
  
  var_Start <- mean_Start*(1-mean_Start)
  var_End <- mean_End*(1-mean_End)
  
  if (conv_Start==0){
    assessed_uplift<-0
  }else{
    assessed_uplift<-(conv_End-conv_Start)/conv_Start
  }
  #assessed_value<-assessed_uplift*100*10^6
  
  
  ### APPLY DIFFERENT STATISTIC TESTS ### 
  
  if (methodology_used == "frequentist") {
    output <- frequentist_test(visit_Start,visit_End,mean_Start,mean_End,var_Start,var_End,alpha,confidence,assessed_uplift)
    return(output)
    
  }else if (methodology_used == "bayesian") {
    output <- bayesian_test(conv_Start,conv_End,fail_Start,fail_End,confidence,assessed_uplift)
    return(output)
    
  }else if (methodology_used == "chi2"){
    output <- chi2_test(visit_Start,visit_End,conv_Start,conv_End,fail_Start,fail_End,mean_Start,mean_End,tot_vis,tot_conv,tot_fail,alpha,assessed_uplift)
    return(output)
    
  }else if (methodology_used == "bounds"){ 
    N_week<-N_week
    output <- bounds_test(visit_Start,visit_End,mean_Start,mean_End,var_Start,var_End,N_week,alpha,confidence,assessed_uplift,Bounds = Bounds,esti_obs,N_week_max)
    return(output)
    
  }else if (methodology_used == "freq/bounds"){ 
    N_week<-N_week
    output <- freq_bounds_test(visit_Start,visit_End,mean_Start,mean_End,var_Start,var_End,N_week,alpha,confidence,assessed_uplift,Bounds = Bounds,esti_obs,N_week_max)
    return(output)
    
  }else if (methodology_used=="permutation"){
    Start_sample<-Start_distrib
    End_sample<-End_distrib
    output<-permutation_test(Start_sample,End_sample,alpha,assessed_uplift)
    return(output)
  }else
    Start_sample<-Start_distrib
  End_sample<-End_distrib
  output<-bootstrap_test(Start_sample,End_sample,alpha,assessed_uplift) 
  return(output)
  
}



###### Theoretical sample size (Assuming a normal distribution) ######

sample_size <- function (m1,m2,alpha,beta){
  q1<-1-m1
  q2<-1-m2
  M<-(m1+m2)/2
  Q<-1-M
  Zalpha <- qnorm(1-(alpha/2), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  Zbeta <- -qnorm(1-beta, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  
  N <- ((Zalpha*((2*M*Q)^0.5)+Zbeta*((m1*q1+m2*q2)^0.5))^2)/((m2-m1)^2)
  return(N)
}


######   MAIN LOOP  ######

main_loop <- function (test_list,output_results_empty) {
  
  ## LOOP THROUGH ALL THE SIMULATIONS ##  
  
  ## FOR DIFFERENT ALPHA VALUES ##  
  
  for (alpha in c(0.01, 0.05,0.1,0.15,0.2,0.25)){
    Bounds <- gsDesign(test.type=3, alpha = alpha, beta = 0.2, k=10, sfu=sfPower, sfupar=3, sfl=sfPower, sflpar=2)
    
    for (test_nb in 1:(length(test_list)-1)){
      test<-test_list[test_nb,] #select the characteristics of a specific test
      #print(test_nb)
      #print(test)
      #print(is.data.frame(test))
      #print(test[1,3])
      
      
      ## REPEAT EACH SIMULATION 10 TIMES ##  
      
      for (version in c(1:10)) { #repeat the simulation 10 times
        Start_obs<- test[1,3]
        End_obs<-test[1,4]
        Start_distrib <- rbinom(Start_obs,1,test[1,1]) #full ditributions
        End_distrib <- rbinom(End_obs,1,test[1,2])
        #bounds_finish<-FALSE #tells us whether the Bounds algo is already done or not
        
        
        ## APPLY TESTS TO INCREASING SIZES OF THE SIMULATION ##  
        
        for (N_tile in (1:10)) {#assess increasing samples of the Start and End distributions
          Start_sample <- Start_distrib[1:(N_tile*Start_obs/10)] 
          End_sample <- End_distrib[1:(N_tile*End_obs/10)]
          
          
          ## APPLY ALL THE DIFFERENT STATISTIC TESTS ##  
          
          for (methodo in c("frequentist", "bayesian","chi2","freq/bounds")) {
            output <- testing_distribution(Start_sample, End_sample, methodology_used = methodo, alpha = alpha, N_tile,Bounds=Bounds) 
            theo_sample_size<-sample_size(test[1,1],test[1,2],alpha,0.8)
            output_df <-  cbind(test_nb,version,N_tile,methodo,alpha,output,theo_sample_size)
            print(output_df)
            output_results_empty <- rbind(output_results_empty,output_df)
          }
          
        }
        
      }
      
      
    }
    
  }
  
  return (output_results_empty)
  
  
}



main_table <- function (){
  
  ######   INITIALISE DATAFRAMES  ######
  
  output_results_empty <- data.frame(Test_id = c()
                                     , version = c() 
                                     , decile =c()
                                     , methodology = c()
                                     , alpha=c()
                                     , confidence = c()
                                     , parameters=c()
                                     , assessed_output = c()
                                     , assessed_uplift = c())
  
  test_list <- data.frame(Start_mean = c(0.01, 0.02, 0.03)
                          , End_mean = c(0.0103, 0.0203, 0.0305)
                          , Start_obs = 1000000
                          , End_obs = 1000000)
  
  output_results <- main_loop(test_list, output_results_empty)
  
  N_rows<-dim(output_results)[1]
  
  N_obs <- rep(test_list[1,3],N_rows)
  output_results[,"Nb Obs"] <- N_obs
  
  A_mean <- test_list[,1]
  B_mean <- test_list[,2]
  
  Theo_uplift <- rep(c((B_mean[1]-A_mean[1])/A_mean[1],(B_mean[2]-A_mean[2])/A_mean[2],(B_mean[3]-A_mean[3])/A_mean[3]),each=(N_rows/length(A_mean)))
  
  Theo_output <- rep(c(as.numeric(B_mean[1]>A_mean[1]),as.numeric(B_mean[2]>A_mean[2]),as.numeric(B_mean[3]>A_mean[3])),each=(N_rows/length(A_mean)))
  output_results[,"Theo output"]<-Theo_output
  
  S_mean <- rep(c(A_mean[1],A_mean[2],A_mean[3]),eagch=(N_rows/length(A_mean)))
  output_results[,"Start CR"]<-S_mean
  E_mean <- rep(c(B_mean[1],B_mean[2],B_mean[3]),each=(N_rows/length(A_mean)))
  output_results[,"End CR"]<-E_mean
  
  output_results[,"Theo uplift"]<-Theo_uplift
  
  return(output_results)
}
#View(output_results)

#output<-main_table()

#write.csv(output,"AB Test simulation 2.csv",row.names = FALSE)




