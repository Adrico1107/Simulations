source(file="StatsSim.R")
source(file="MainSim.R")

library(ggplot2)

type <- "Simulation"

Obs_week <- 10000/2

#THE RATIO MAX_OBS/inc_size HAS TO BE BETWEEN 2 AND 19 (2 AND 13 FOR ALPHA=0.01)

#### WE NEED TO HAVE MAX_OBS/inc_size <19 OTHERWISE FREQ/BOUNDS CANNOT HANDLE IT

###### ESTIMATED NUMBER OF OBSERVATIONS NEEDED ######

esti_N_obs <- function(Start_mean,theo_uplift,alpha){
  esti_uplift <- (abs(theo_uplift)*runif(1,0.5,1.5)) #estimated uplift by the Pierre et Vacances teams
  End_mean <- Start_mean*(1 + esti_uplift) #MUST BE SMALLER THAN 1 !!!!!
  #print(c(Start_mean,End_mean,alpha))
  return(sample_size(Start_mean,End_mean,alpha,0.8))
}

###### DECISION MAKING ######

roll_out <- function(conv_Start,conv_End,methodology,current_obs,alpha,assessed_output,esti_obs){
  #determines whether we can kill the test and its outcome
  #esti_obs <- esti_N_obs(Start_mean,theo_uplift,alpha)
  
  if ((current_obs<=Obs_week)|(conv_Start<100)|(conv_End<100)){ #Basic rules for a usual AB test
    return("No decision")
  }else{
    if (assessed_output==0){
      return("No decision")
    }else if (assessed_output==1){
      return("Roll out")
    }else{
      return("No roll out")
    }
  }
}



###### CUMULATED RESULTING UPLIFT ######

cumulated_uplift <- function (output_results_empty,boxplot_empty,time){
 
    ## FOR DIFFERENT ALPHA VALUES ##  

  
  for (version in c(1:1)) {  
    
    ### Initialize the A and B CTR distributions ###
    
    uplift_list<-runif(30,0.94,1.06) #theoretical uplift
    
    A<-(runif(30,0.005,0.015))
    B<-(A*uplift_list) #random list of 100 AB tests
    
    distrib_list_Start <- vector("list", 30)
    distrib_list_End <- vector("list", 30)
    for (j in 1:length(A)){
      Start_mean <- A[j]
      End_mean <- B[j]
      theo_uplift <- ((End_mean-Start_mean)/Start_mean)
      inc <- floor(5000+(j-1)*((500000-5000)/(length(A)-1))) #nombre de visites par semaine augmente à pas constant entre 10k et 1M en tout (idem/2 pour chaque page)
      #inc<-(1400000/2)
      max_obs <- inc * 12 #12 semaines d'observation au max
      #alpha <- 0.05 #On prend le alpha le plus petit pour avoir les distributions de plus grande taille possible
      
      #esti_obs_dist <- min(floor(sample_size(Start_mean,Start_mean*(1+1.5*abs(theo_uplift)),alpha,0.8)),max_obs) 
      esti_obs_dist <- max_obs
      #upper bound of estimated sample size 

      Start_distrib <- rbinom(esti_obs_dist,1,Start_mean) #full ditributions
      End_distrib <- rbinom(esti_obs_dist,1,End_mean)
      
      distrib_list_Start[[j]] <- Start_distrib
      distrib_list_End[[j]] <- End_distrib
    }
    
   for (alpha in c(0.05,0.1,0.15,0.2,0.3,0.5)){
          
     
     ### Initialize the test list ###
     
          Esti_obs <- vector()
          Inc_size <- vector()
          for (i in 1:length(A)){
            Start_mean <- A[i]
            End_mean <- B[i]
            theo_uplift <- ((End_mean-Start_mean)/Start_mean)
            inc <- floor(5000+(i-1)*((500000-5000)/(length(A)-1))) #nombre de visites par semaine augmente à pas constant entre 10k et 1M en tout (idem/2 pour chaque page)
            #inc<-(1400000/2)
            max_obs <- inc * 12 #12 semaines d'observation au max
            Esti_obs <- c(Esti_obs,min(floor(esti_N_obs(Start_mean,theo_uplift,alpha)),max_obs)) #estimated sample size cannot be greater than the nb of elements in the sample (according to an estimated uplift)
            Inc_size <- c(Inc_size,inc)
          }
          

          rd_test_list <- data.frame(A,B,Esti_obs,Inc_size)
                      
          
          #print(rd_test_list)
          
          for (methodology in c("full_bayes","full_freq","frequentist","bayesian","chi2","bounds","permutation","bootstrap")){ #

            cumul_ideal_uplift<-0
            cumul_assessed_uplift <- 0  
            cumul_real_uplift<-0 #theo_uplift when roll out
            cumul_good_output <- 0
            cumul_obs<-0
            cumul_esti_obs<-0
            cumul_weeks<-0
            cumul_esti_weeks<-0
            
            print(c(version,alpha,methodology))
          
            for (test_id in 1:(nrow(rd_test_list))){
            
              
              
            test<-rd_test_list[test_id,] #select the characteristics of a specific test
          
            Start_mean<-test[1,1]
            End_mean<-test[1,2]
            esti_obs <- test[1,3]
            inc_size <- test[1,4]
            
          
            
            N_week_max <- ceiling(esti_obs/inc_size)
            
            theo_uplift <- ((End_mean-Start_mean)/Start_mean)
            theo_output<-sign(theo_uplift)
            
            final_decision <- roll_out(1000,1000,"frequentist",10000,0.05,theo_output,esti_obs)
            
            
            if (theo_output==1){
              cumul_ideal_uplift<-(cumul_ideal_uplift+theo_uplift) #real uplift if we had rolled out only the ones with a positive uplift (ideal world with no mistakes)
            }
            
            
            cumul_esti_obs<-(cumul_esti_obs+esti_obs)
            cumul_esti_weeks<-ceiling(cumul_esti_weeks+(esti_obs/inc_size))#we check at most once every week hence the ceiling function
            #print(c(theo_uplift,Start_mean,End_mean,esti_obs))
            
            Start_distrib <- distrib_list_Start[[test_id]] #full ditributions
            End_distrib <- distrib_list_End[[test_id]]
            Bounds<-0
            
            
            if ((methodology=="full_freq")||(methodology=="full_bayes")){ #We compare to the results of a traditional frequentist A/B test
              
              Start_sample <- Start_distrib[1:esti_obs]
              End_sample <- End_distrib[1:esti_obs]
              
              outcome <- TRUE
              conv_Start<-sum(Start_sample)
              conv_End<-sum(End_sample)
              
              current_obs<-esti_obs
              
              N_week<-0
              if (methodology=="full_freq"){
                output <- testing_distribution(Start_sample, End_sample, methodology_used="frequentist", alpha , N_week=0,Bounds=0,esti_obs,N_week_max=0)
                
              }else{
                output <- testing_distribution(Start_sample, End_sample, methodology_used="bayesian", alpha , N_week=0,Bounds=0,esti_obs,N_week_max=0)
                
              }
              assessed_output<-output$assessed_output
              decision<-roll_out(conv_Start,conv_End,methodology,current_obs,alpha,assessed_output,esti_obs)
              
              #print(c(version,alpha,methodology,test_id,assessed_output,theo_output,decision,conv_Start,conv_End))
              
              if (decision=="Roll out"){
                cumul_assessed_uplift<-(cumul_assessed_uplift+output$assessed_uplift)#assessed uplifts of roll outs (increases over time)
                cumul_real_uplift<-(cumul_real_uplift+theo_uplift) #real uplift of tests rolled out (can increase if the test is correct or decrease if not over time)
              }
              
              
              if (assessed_output==theo_output){
                cumul_good_output <- (cumul_good_output+1)
              }
              
              
              cumul_weeks<-cumul_esti_weeks
              cumul_obs<-cumul_esti_obs
              
              output_df <-  cbind(type,alpha,version,methodology,test_id,N_week,current_obs,inc_size,Start_mean,End_mean,esti_obs,esti_up=0,conv_Start
                                  ,conv_End,output,theo_output,theo_uplift,decision,cumul_good_output,cumul_obs,cumul_weeks,cumul_esti_obs
                                  ,cumul_esti_weeks
                                  ,cumul_ideal_uplift,cumul_assessed_uplift,cumul_real_uplift,outcome
                                  , day_nb=N_week*7, cumul_days=cumul_weeks*7
                                  , current_obs_S=current_obs, current_obs_E=current_obs, cumul_obs_S=cumul_obs, cumul_obs_E=cumul_obs
                                  , current_conv_S=conv_Start, current_conv_E=conv_End, final_decision)

              output_results_empty <- rbind(output_results_empty,output_df)
              
              
              }
            
            if ((methodology=="freq/bounds")|(methodology=="bounds")){
              N_week_max<-max(2,min(19,ceiling(esti_obs/inc_size))) #HAS TO BE BETWEEN 2 AND 19 FOR ALPHA>0.05 AND 2 AND 13 FOR ALPHA=0.01
              Bounds <- gsDesign(test.type=3, alpha = alpha, beta = 0.2, k=N_week_max, sfu=sfPower, sfupar=3, sfl=sfPower, sflpar=2)
            }
            
            outcome<-FALSE
            N_week<-1
            current_obs<-floor(min(inc_size,esti_obs)) #WE HAVE TO CHECK AT LEAST ONCE
            delta_obs<-current_obs
            
            
            
            while ((outcome==FALSE) && (N_week<=12) && (methodology!="full_freq") && (methodology!="full_bayes")){#assess increasing samples of the Start and End distributions until a decision is made

              Start_sample <- Start_distrib[1:min(current_obs,esti_obs)] #we go to the very end if the test sample if needed
              End_sample <- End_distrib[1:min(current_obs,esti_obs)]
              
              #print(c(length(Start_distrib),current_obs,esti_obs))
            
              conv_Start<-sum(Start_sample)
              conv_End<-sum(End_sample)
              
              #print(c(conv_Start,conv_End))
            
            ## APPLY ALL THE DIFFERENT STATISTIC TESTS ##  

              output <- testing_distribution(Start_sample, End_sample, methodology_used = methodology, alpha = alpha, N_week,Bounds=Bounds
                                             ,esti_obs,N_week_max=0)
              assessed_output<-output$assessed_output
            
              decision<-roll_out(conv_Start,conv_End,methodology,current_obs,alpha,assessed_output)

              cumul_obs<-(cumul_obs+delta_obs)
              cumul_weeks<-(cumul_weeks+1)
              

              
              if ((assessed_output==theo_output)&&((decision!="No decision")|((decision=="No decision")&&(current_obs==esti_obs)))){
                cumul_good_output <- (cumul_good_output+1)
              }
              
              #print(c(N_week,current_obs,assessed_output,theo_output,decision))
                

              if (decision=="Roll out"){
                  cumul_assessed_uplift<-(cumul_assessed_uplift+output$assessed_uplift)#the uplift only really changes when we roll out a change
                  cumul_real_uplift<-(cumul_real_uplift+theo_uplift)
              }
              
              if ((decision != "No decision")|(current_obs==esti_obs)){ #We have either reached a conclusive decision or the end of the test
                outcome <- TRUE
              }
              
              

              output_df <-  cbind(type,alpha,version,methodology,test_id,N_week,current_obs,inc_size,Start_mean,End_mean,esti_obs,esti_up=0,conv_Start
                                  ,conv_End,output,theo_output,theo_uplift,decision,cumul_good_output,cumul_obs,cumul_weeks,cumul_esti_obs
                                  ,cumul_esti_weeks
                                  ,cumul_ideal_uplift,cumul_assessed_uplift,cumul_real_uplift,outcome, day_nb=N_week*7, cumul_days=cumul_weeks*7
                                  , current_obs_S=current_obs, current_obs_E=current_obs, cumul_obs_S=cumul_obs, cumul_obs_E=cumul_obs
                                  , current_conv_S=conv_Start, current_conv_E=conv_End, final_decision)
              
              output_results_empty <- rbind(output_results_empty,output_df)

            
                
                N_week<-N_week+1
                current_obs_prec<-current_obs
                current_obs<-min(current_obs+inc_size,esti_obs) #we can't observe more than what was estimated
                delta_obs<-ifelse(current_obs==esti_obs,esti_obs-current_obs_prec,inc_size)
                
            

            }
            
          
            }
            print(c(cumul_weeks,cumul_esti_weeks,cumul_assessed_uplift,cumul_real_uplift,cumul_ideal_uplift,(Sys.time()-time)))
            
            yearly_assessed_uplift <- ((52/cumul_weeks)*cumul_assessed_uplift)
            yearly_real_uplift <- ((52/cumul_weeks)*cumul_real_uplift) #uplift sur un an normalisé pour comparer entre versions
            yearly_ideal_uplift <- ((52/cumul_weeks)*cumul_ideal_uplift)

            box_output <- data.frame(type,alpha,version,methodology,cumul_obs,cumul_weeks,cumul_assessed_uplift,cumul_real_uplift
                                     ,cumul_ideal_uplift,yearly_assessed_uplift,yearly_real_uplift,yearly_ideal_uplift
                                     ,cumul_good_output, cumul_days=cumul_weeks*7, cumul_max_days=0, cumul_obs_S=cumul_obs
                                     , cumul_obs_E=cumul_obs)
            
            #names(box_output)<-c("alpha","version","methodology","total observations","total weeks","cumul assessed uplift"
                                # ,"cumul real uplift","cumul ideal uplift","yearly assessed uplift","yearly real uplift","yearly ideal uplift")
            #print(box_output)
            boxplot_empty <- rbind(boxplot_empty,box_output)
          
      }
     }
      
   } 
  big_list <- list("dataframe" = output_results_empty, "boxplot" = boxplot_empty)  
  return(big_list)
    
}



dataset <- function(){
  
  time <- Sys.time()
  output_results_empty <- data.frame( type=c(),
                                      alpha=c()
                                     , version = c() 
                                     , methodology = c()
                                     , test_id=c()
                                     , N_week=c()
                                     , current_obs =c()
                                     , inc_size=c()
                                     , Start_mean=c()
                                     , End_mean=c()
                                     , esti_obs=c()
                                     , esti_up=c()
                                     , conv_Start=c()
                                     , conv_End=c()
                                     , confidence = c()
                                     , parameters=c()
                                     , assessed_output = c()
                                     , assessed_uplift = c()
                                     , theo_output=c()                                     
                                     , theo_uplift=c()
                                     , desision=c()
                                     , cumul_good_output=c()
                                     , cumul_obs=c()
                                     , cumul_weeks=c()
                                     , cumul_esti_obs=c()
                                     , cumul_esti_weeks=c()
                                     , cumul_ideal_uplift=c()
                                     , cumul_assessed_uplift=c()
                                     , cumul_real_uplift=c()
                                     , outcome=c()
                                     , day_nb=c()
                                     , cumul_days=c()
                                     , current_obs_S=c()
                                     , current_obs_E=c()
                                     , cumul_obs_S=c()
                                     , cumul_obs_E=c()
                                     , current_conv_S=c()
                                     , current_conv_E=c()
                                     , final_decision=c())


  boxplot_empty <- data.frame(  type=c()
                              , alpha=c()
                              , version = c()
                              , methodology = c()
                              , cumul__observations=c()
                              , cumul_weeks=c()
                              , cumul_assessed_uplift=c()
                              , cumul_real_uplift=c()
                              , cumul_ideal_uplift=c()
                              , yearly_assessed_uplift=c()
                              , yearly_real_uplift=c()
                              , yearly_ideal_uplift=c()
                              , cumul_good_output=c()
                              , cumul_days=c()
                              , cumul_max_days=c()
                              , cumul_obs_S=c()
                              , cumul_obs_E=c())
  

  results <- cumulated_uplift(output_results_empty,boxplot_empty,time)
  print(Sys.time()-time)
  return(results)


}


#print(distrib_list_Start)

results1 <- dataset()

output_results_final4 <- results1$dataframe
View(output_results_final4)

#print(ncol(output_results_final1))

write.csv(output_results_final1,"D:/Stage DataMa/Final datasets/Simulated AB tests results 1 version.csv",row.names = FALSE)

#real_results <- as.data.frame(read.csv("D:/Stage DataMa/Final datasets/Real AB tests results final.csv"))
#results_final <- rbind(output_results_final1,real_results)
#write.csv (results_final,"D:/Stage DataMa/Final datasets/Final dataset Sim+Real.csv",row.names = FALSE)

boxplot_final4 <- results1$boxplot
View(boxplot_final1)
write.csv(boxplot_final1,"D:/Stage DataMa/Final datasets/Simulated AB tests boxplot 1 version.csv",row.names = FALSE)

#real_boxplot <- as.data.frame(read.csv("D:/Stage DataMa/Final datasets/Real AB tests boxplot final.csv"))
#boxplot_final <- rbind(boxplot_final1,real_boxplot)
#write.csv (boxplot_final,"D:/Stage DataMa/Final datasets/Final boxplot Sim+Real.csv",row.names = FALSE)



"
#### MIX THE 2 SIMULATIONS (+ CHANGE VERSION NB + GET RID OF FREQ/BOUNDS) ####

results_final <- output_results_final1
for (i in 1:nrow(results_final)){
  results_final[i,3] <- results_final[i,3]+3
}

#View(results_final)
results_real <- as.data.frame(read.csv(D:/Stage DataMa/Final datasets/Final dataset test.csv))
results_real <- results_real %>% filter(methodology!=freq/bounds)
#View(results_real)

results_final <- rbind (results_final,results_real)

write.csv (results_final,D:/Stage DataMa/Final datasets/Final dataset test total 1.csv,row.names = FALSE)

boxplot_final <- boxplot_final1
for (i in 1:nrow(boxplot_final)){
  boxplot_final[i,3] <- boxplot_final[i,3]+3
}
View(boxplot_final)
boxplot_real <- as.data.frame(read.csv(D:/Stage DataMa/Final datasets/Final boxplot test.csv))
boxplot_real <- boxplot_real %>% filter(methodology!=freq/bounds)
boxplot_final <- rbind(boxplot_final,boxplot_real)
#write.csv (boxplot_final,D:/Stage DataMa/Final datasets/Final boxplot test total 1.csv,row.names = FALSE)

"
