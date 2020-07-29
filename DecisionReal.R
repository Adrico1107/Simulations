source(file="StatsReal.R")
source(file="MainReal.R")
#source(file="DecisionSimFinal.R")

#library(ggplot2)

Obs_week <- 10000/2

#THE RATIO MAX_OBS/inc_size HAS TO BE BETWEEN 2 AND 19 (2 AND 13 FOR ALPHA=0.01)

#### WE NEED TO HAVE MAX_OBS/inc_size <19 OTHERWISE FREQ/BOUNDS CANNOT HANDLE IT


###### DECISION MAKING ######

roll_out <- function(conv_Start,conv_End,alpha,assessed_output,day_nb,methodology){
  #determines whether we can kill the test and its outcome
  #esti_obs <- esti_N_obs(Start_mean,theo_uplift,alpha)
  
  if ((day_nb<=7)|(conv_Start<100)|(conv_End<100)){ #Basic rules for a usual AB test
    return("No decision")
  }else{
      if (assessed_output==1){
        return("Roll out")
      }else if (assessed_output==-1){
        return("No roll out")
      }else{
        return("No decision")
      }
  }
}


###### IMPORT THE DATA ######

library(tidyverse)

real_data0 = read.csv("Real AB tests v4.txt",header=TRUE,sep="\t") #import the real data

real_data <-  ((real_data0 %>% filter(Event_Label!="D-615700-Version centrée")) %>% filter(Event_Label!="C-615699-Version claire")) %>% filter(Event_Label!="C-537119-No offer block")

#View(real_data)
##Replace with Start and End##

Event_Label_clean <- as.vector(real_data[,3])
Event_Label_clean <- gsub(".*-A$|^A-.*|^0-.*", "Start", Event_Label_clean, perl=TRUE)
Event_Label_clean <- gsub(".*-B$|^B-.*|.*-Variation 1$|^C-537119.*", "End", Event_Label_clean, perl=TRUE)
real_data[,3] <- Event_Label_clean


#real_data

##Reorder each test by date##

test_name_list <- as.vector(unique(real_data[,1]))
#test_name_list

clean_data <- data.frame(Event_Action=c()
                         , Date=c()
                         , Event_Label=c()
                         , Sessions=c()
                         , Transactions=c())

for (test_name in test_name_list){
  new_data <- (real_data %>% filter(Event_Action==test_name))%>% arrange(Date) #subset of the real_data with rows in ascending order (by date)
  
  "Date <- as.vector(new_data[,2]) #solving the date problem
  Date_clean<-vector()
  for (i in 2:length(Date_clean)){ 
    Date_clean[i] <- Date_clean[i]-baseline+1
    if (Date_clean[i]>=8870){
      Date_clean[i] <- Date_clean[i]-8869
    }
  }"
  N<-(dim(new_data))[1]/2
  new_data[,2]<-rep(c(1:N),each=2)
  clean_data <- rbind(clean_data,  new_data) #append to the new dataset
}

#View(clean_data)

###### CUMULATED RESULTING UPLIFT ######

cumulated_uplift_real <- function (real_results_empty,real_boxplot_empty,time){
  
  type <- "Real"
  
## Real test characteristics ##
  
  Test_desc <- data.frame(test_nb=c(),test_name=c(),N_days=c(),obs_S=c(),conv_S=c(),obs_E=c(),conv_E=c(),theo_uplift=c())
  test_nb <- 0
  
  for (test_name in test_name_list){
    test_nb <- test_nb +1
    obs_S <- sum((clean_data %>% filter(Event_Action==test_name,Event_Label=="Start"))$Sessions) #sum of all sessions on Start page for a certain test
    obs_E <- sum((clean_data %>% filter(Event_Action==test_name,Event_Label=="End"))$Sessions) #same for End page
    conv_S <- sum((clean_data %>% filter(Event_Action==test_name,Event_Label=="Start"))$Transactions) #sum of all transactions on Start page for a certain test
    conv_E <- sum((clean_data %>% filter(Event_Action==test_name,Event_Label=="End"))$Transactions) #same for End page
    CTR_S <- conv_S/obs_S
    CTR_E <- conv_E/obs_E
    theo_uplift <- (CTR_E-CTR_S)/CTR_S
    
    N_days <- max((clean_data %>% filter(Event_Action==test_name,Event_Label=="End"))$Date) #-min((clean_data %>% filter(Event_Action==test_name,Event_Label=="End"))$Date)+1 #number of days the test lasted
    Test_desc <- rbind(Test_desc,data.frame(test_nb,test_name,N_days,obs_S,conv_S,obs_E,conv_E,theo_uplift))
  }
  
  test_name_pos <- c("473666-[PROD][AB][PV][FR][ABT][MOB][50/50][SEARCH]Mini tiles V2"
  ,"AB Test Split | OnePage FP - MVP4"
  ,"AB Test Split | Elastic Lot 3"
  ,"AB Test Split | New Funnel CP"
  ,"420651-[ANIM][AB][CP][FR][ABT][MOB][50/50][SCK]no offer block on mobile"
  ,"490368-[PROD][AB][PV][FR][ABT][MOB][50/50][FP]Header Lightering V2")
  
  test_name_neg <- c("AB Test Split | Payment - Price in the Button (Language=nl-nl)"
                     ,"AB Test Split | New Funnel CP Lot 3.1"
                     ,"475942-[ANIM][AB][CP][FR][ABT][DESK][50/50][FP]Video FR park pages LA"
                     ,"AB Test Split | Contact Form Design V3"
                     ,"AB Test Split | SEARCH Target2Sell")

  
    for (alpha in c(0.05,0.1,0.15,0.2,0.3,0.5)){
      
      
      for (methodology in c("full_bayes","full_freq","frequentist","bayesian","chi2","bounds","permutation","bootstrap")){
        
        
        cumul_assessed_uplift <- 0 
        cumul_real_uplift <- 0
        cumul_obs_S<-0
        cumul_obs_E<-0
        cumul_days<-0
        cumul_good_output <- 0
        #cumul_week <- 0

          "Start_mean<-test[1,1]
          End_mean<-test[1,2]
          esti_obs <- test[1,3]
          inc_size <- test[1,4]"
          
          #print(c(alpha,version,methodology,test_nb))
          
          "theo_uplift <- ((End_mean-Start_mean)/Start_mean)
          theo_output<-sign(theo_uplift)
          
          Bounds<-0"
         
        
        for (test_id in test_name_list){
          
          test<-clean_data %>% filter(Event_Action==test_id) #select the data of a specific test
          final_decision <- "No decision"
          theo_uplift <- (Test_desc %>% filter(test_name==test_id))$theo_uplift
          
          if (test_id %in% test_name_pos){
            cumul_real_uplift <- cumul_real_uplift + theo_uplift
            final_decision <- "Roll out"
          }else if (test_id %in% test_name_neg){
            final_decision <- "No roll out"
          }
        
          N_days_max <- (Test_desc %>% filter(test_name==test_id))$N_days
          N_max_bounds <-0
          conv_Start <- (Test_desc %>% filter(test_name==test_id))$conv_S
          conv_End <- (Test_desc %>% filter(test_name==test_id))$conv_E
          
          #print(c(conv_Start,conv_End))
          
          #final_output <- roll_out(conv_Start,conv_End,alpha,final_uplift,N_days_max,methodology)
          
          visit_Start <- (Test_desc %>% filter(test_name==test_id))$obs_S
          visit_End <- (Test_desc %>% filter(test_name==test_id))$obs_E
          
          test_nb <- (Test_desc %>% filter(test_name==test_id))$test_nb
          
          outcome <- FALSE
          

          
          #Start_distrib <- c(rep(1,conv_Start),rep(0,visits_Start-conv_Start))
          #End_distrib <- c(rep(1,conv_End),rep(0,visits_End-conv_End))
          
          
          
          ##### CLASSIC A/B FREQUENTIST TES #####
          
          
          if ((methodology=="full_freq")|(methodology=="full_bayes")){ #We compare to the results of a traditional frequentist/bayesian A/B test
            
            outcome <- TRUE
            
            #conv_Start,conv_End,visit_Start,visit_End, methodology_used, alpha , day_nb,Bounds,N_days_max
            if (methodology=="full_freq"){
              #print(c(conv_Start,conv_End,visit_Start,visit_End, "frequentist", alpha , day_nb,0,0,N_days_max))
            output <- testing_distribution(conv_Start,conv_End,visit_Start,visit_End, methodology_used="frequentist", alpha , day_nb,Bounds=0,N_max_bounds=0,N_days_max)
            }else if (methodology=="full_bayes"){
              output <- testing_distribution(conv_Start,conv_End,visit_Start,visit_End, methodology_used="bayesian", alpha , day_nb,Bounds=0,N_max_bounds=0,N_days_max)
            }
            
            #print(output)
            
            "output <- testing_distribution(9909,8689,1982987,1419239, methodology_used=bayesian, 0.05 , 100,Bounds=0,N_max_bounds=0,100)
            print(output)
            assessed_output<-output$assessed_output
            assessed_uplift<-output$assessed_uplift
            decision<-roll_out(9909,8689,0.05,assessed_output,65,bayesian)"
            
            assessed_output<-output$assessed_output
            assessed_uplift<-output$assessed_uplift
            decision <- roll_out(conv_Start,conv_End,alpha,assessed_output,N_days_max,methodology)
            
            #print(c(assessed_output,decision))
            
            if (decision=="Roll out"){
              cumul_assessed_uplift<-(cumul_assessed_uplift+assessed_uplift)#assessed uplifts of roll outs (increases over time)
            }
            
            if (decision==final_decision){
              cumul_good_output <- cumul_good_output+1
            }
            
            day_nb <- N_days_max
            current_obs_S<-visit_Start
            current_obs_E<-visit_End
            cumul_days<-cumul_days+N_days_max
            cumul_obs_S<-cumul_obs_S+visit_Start
            cumul_obs_E<-cumul_obs_E+visit_End
            N_week<-N_days_max/7
            
            #print(c(test_nb,alpha,methodology,day_nb,N_days_max,assessed_uplift,theo_uplift,current_obs_S,cumul_obs_S))
            print(c(test_nb,methodology, alpha, day_nb,N_days_max))
            
            #print(theo_uplift)
            
            "A <- c(type,alpha,version=0,methodology,test_id,N_week=0,current_obs=0,inc_size=0,Start_mean=0,End_mean=0,esti_obs=0,esti_up=0,
                    conv_Start,conv_End,
                    output,theo_output=0,theo_uplift,decision,cumul_good_output,cumul_obs=0,cumul_weeks=0,cumul_esti_obs=0,cumul_esti_weeks=0,
                    cumul_ideal_uplift=0,cumul_assessed_uplift,cumul_real_uplift,day_nb,cumul_days,current_obs_S, current_obs_E,
                    cumul_obs_S,cumul_obs_E,current_conv_S=conv_Start,current_conv_E=conv_End,final_decision)
           
             
           
            for (x in A){
              print(x)
            }"
            
            #print(names(A))
            #print(names(real_results_empty))
            
            #print(nrow(type))
            #print(nrow(alpha))
            
            output_df <-  cbind(type,alpha,version=0,methodology,test_id,N_week=0,current_obs=0,inc_size=0,Start_mean=0,End_mean=0,esti_obs=0,esti_up=0,
                                conv_Start,conv_End,
                                output,theo_output=0,theo_uplift,decision,cumul_good_output,cumul_obs=0,cumul_weeks=0,cumul_esti_obs=0,cumul_esti_weeks=0,
                                cumul_ideal_uplift=0,cumul_assessed_uplift,cumul_real_uplift,outcome,day_nb,cumul_days,current_obs_S, current_obs_E,
                                cumul_obs_S,cumul_obs_E,current_conv_S=conv_Start,current_conv_E=conv_End,final_decision)
            
            #names(output_df) <- c("type","alpha","version","methodology","test_nb","N_week","current_obs", "inc_size", "Start_mean", "End_mean"
             #                     , "esti_obs", "esti_up", "conv_Start", "conv_End" , "confidence", "parameters"
              #                    , "assessed_output", "assessed_uplift", "theo_output" , "theo_uplift", "decision", "cumul_good_output"
               #                   , "cumul_obs", "cumul_weeks", "cumul_esti_obs" , "cumul_esti_weeks", "cumul_ideal_uplift", "cumul_assessed_uplift"
                #                  , "cumul_real_uplift", "day_nb", "cumul_days"
                 #                 , "current_obs_S" , "current_obs_E" , "cumul_obs_S" , "cumul_obs_E", "current_conv_S" 
                  #                , "current_conv_E" , "final_decision")
            
            #print("YA")
            real_results_empty <- rbind(real_results_empty,output_df)
          }
          
          
          
          ##### SHORT TESTS #####
          
          
          if ((methodology=="freq/bounds")|(methodology=="bounds")){
            N_max_bounds <- ifelse(alpha==0.5,max(2,min(14,N_days_max)),max(2,min(19,N_days_max))) #HAS TO BE BETWEEN 2 AND 19 FOR ALPHA>0.05 AND 2 AND 13 FOR ALPHA=0.01
            Bounds <- gsDesign(test.type=3, alpha = alpha, beta = 0.2, k=N_max_bounds, sfu=sfPower, sfupar=3, sfl=sfPower, sflpar=2)
            plot(Bounds)
          }
          
                    
          #outcome<-FALSE
          #N_week<-1
          day_nb<-1
          current_obs_S <- 0
          current_obs_E <- 0
          current_conv_S<-0
          current_conv_E<-0
          #delta_obs_S<-current_obs_S
          #delta_obs_E<-current_obs_E
          
          
          
          while ((outcome==FALSE) & (day_nb<=N_days_max) & (methodology!="full_freq")& (methodology!="full_bayes")){#assess increasing samples of the Start and End distributions until a decision is made
            
            
            current_conv_S<-current_conv_S + (clean_data %>% filter(Event_Action==test_id,Date==day_nb,Event_Label=="Start"))$Transactions
            current_conv_E<-current_conv_E + (clean_data %>% filter(Event_Action==test_id,Date==day_nb,Event_Label=="End"))$Transactions
            
            #print(current_conv_S)
            
            current_obs_S<- current_obs_S + (clean_data %>% filter(Event_Action==test_id,Date==day_nb,Event_Label=="Start"))$Sessions
            current_obs_E<- current_obs_E + (clean_data %>% filter(Event_Action==test_id,Date==day_nb,Event_Label=="End"))$Sessions
            
            #visit_Start <- visit_Start + (clean_data %>% filter(test_name==test_id,N_days==day_nb))$obs_S
            #visit_End <- visit_End + (clean_data %>% filter(test_name==test_id,N_days==day_nb))$obs_E
            
            ## APPLY ALL THE DIFFERENT STATISTIC TESTS ## 
            #plot(Bounds)
            
            if ((methodology=="bounds")|(methodology=="freq/bounds")){
              output <- testing_distribution(current_conv_S,current_conv_E,current_obs_S,current_obs_E, methodology_used = methodology, alpha = alpha, day_nb,Bounds,N_max_bounds,N_days_max)
            }else{
              #print(c(current_conv_S,current_conv_E,current_obs_S,current_obs_E))
              
            output <- testing_distribution(current_conv_S,current_conv_E,current_obs_S,current_obs_E, methodology_used = methodology, alpha = alpha, day_nb,0,0,N_days_max)
            }
                        
            print(c(test_nb,methodology, alpha, day_nb,N_days_max))
            
            assessed_output<-output$assessed_output
            assessed_uplift<-output$assessed_uplift
            pval<-output$parameters
            decision<-roll_out(current_conv_S,current_conv_E,alpha,assessed_output,day_nb,methodology)
            
            
            cumul_days<-cumul_days+1
            cumul_obs_S<-cumul_obs_S+(clean_data %>% filter(Event_Action==test_id,Date==day_nb,Event_Label=="Start"))$Sessions
            cumul_obs_E<-cumul_obs_E+(clean_data %>% filter(Event_Action==test_id,Date==day_nb,Event_Label=="End"))$Sessions
            
            
            if (decision=="Roll out"){
              #print("YYYYYYY")
              cumul_assessed_uplift<-(cumul_assessed_uplift+assessed_uplift)#the uplift only really changes when we roll out a change
              outcome <- TRUE
            }
            
            if (decision=="No roll out"){
              #print("NNNNNN")
              outcome <- TRUE
            }
            
            #print(outcome)
            
            if ((outcome==TRUE)&&(decision==final_decision)){ #The test is done early (we have a clear answer)
              cumul_good_output <- cumul_good_output+1
            }
            
            if ((outcome==FALSE)&&(day_nb==N_days_max)&&(decision==final_decision)){#We have reached the end of the test without a clear decision
              cumul_good_output <- cumul_good_output+1
            }
            
            #if (floor(day_nb/7)==(day_nb/7)){
            #  cumul_weeks <- cumul_weeks+1
            #}
            
            #print(c(test_nb,alpha,methodology,day_nb,N_days_max,pval,assessed_output,assessed_uplift,theo_uplift,current_conv_S
             #       ,current_conv_E,current_obs_S,current_obs_E,decision,final_decision,outcome))
            
            output_df <-  cbind(type,alpha,version=0,methodology,test_id,N_week=0,current_obs=0,inc_size=0,Start_mean=0,End_mean=0,esti_obs=0,esti_up=0,
                                conv_Start,conv_End
                                ,output,theo_output=0,theo_uplift,decision,cumul_good_output,cumul_obs=0,cumul_weeks=0
                                ,cumul_esti_obs=0,cumul_esti_weeks=0,cumul_ideal_uplift=0,cumul_assessed_uplift,cumul_real_uplift,outcome
                                ,day_nb,cumul_days,current_obs_S, current_obs_E,cumul_obs_S,cumul_obs_E,current_conv_S,current_conv_E,final_decision)
          
           # names(output_df) <- c("type","alpha","version","methodology","test_nb","N_week","current_obs", "inc_size", "Start_mean", "End_mean"
            #                      , "esti_obs", "esti_up", "conv_Start", "conv_End" , "confidence", "parameters"
             #                     , "assessed_output", "assessed_uplift", "theo_output" , "theo_uplift", "decision", "cumul_good_output"
              #                    , "cumul_obs", "cumul_weeks", "cumul_esti_obs" , "cumul_esti_weeks", "cumul_ideal_uplift", "cumul_assessed_uplift"
               #                   , "cumul_real_uplift", "day_nb", "cumul_days"
                #                  , "current_obs_S" , "current_obs_E" , "cumul_obs_S" , "cumul_obs_E", "current_conv_S" 
                 #                 , "current_conv_E" , "final_decision")
            
            real_results_empty <- rbind(real_results_empty,output_df)
            
            day_nb <- day_nb+1
            
            #if (current_obs==esti_obs){ #We have either reached a decision or the end of the test
             # outcome <- TRUE
              
            }
            
          }
          
          
        
        print(c(cumul_days,cumul_assessed_uplift,cumul_real_uplift,(Sys.time()-time)))
        
        
        cumul_max_days <- sum(Test_desc$N_days)
        yearly_assessed_uplift <- ((365/cumul_days)*cumul_assessed_uplift)
        yearly_real_uplift <- ((365/cumul_max_days)*cumul_real_uplift) #uplift sur un an normalisé pour comparer entre versions
        
        box_output <- data.frame(type,alpha,version=0,methodology,cumul_obs=0,cumul_weeks=0,cumul_assessed_uplift,cumul_real_uplift
                                 ,cumul_ideal_uplift=0,yearly_assessed_uplift,yearly_real_uplift,yearly_ideal_uplift=0,cumul_good_output,
                                 cumul_days
                                 ,cumul_max_days,cumul_obs_S,cumul_obs_E)
        
        #names(box_output)<-c("type","alpha","version","methodology","cumul__observations","cumul_weeks","cumul_assessed_uplift","cumul_real_uplift"
         #                    ,"cumul_ideal_uplift","yearly_assessed_uplift","yearly_real_uplift","yearly_ideal_uplift",
          #                   "cumul_good_output","cumul_days","cumul_max_days"
           #                  ,"cumul_obs_S","cumul_obs_E")
      


        #print(box_output)
        real_boxplot_empty <- rbind(real_boxplot_empty,box_output)
        
      }
    }  
    
  big_list <- list("dataframe" = real_results_empty, "boxplot" = real_boxplot_empty)  
  return(big_list)
  
}



datasetReal <- function(){
  
  time <- Sys.time()
  real_results_empty <- data.frame( type=c()
                                      , alpha=c()
                                      , version = c() 
                                      , methodology = c()
                                      , test_nb=c()
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
                                      , decision=c()
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
  
  
  
  
  
  real_boxplot_empty <- data.frame(  type=c()
                                     ,alpha=c()
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
  
  
  results_sim_real <- cumulated_uplift_real(real_results_empty,real_boxplot_empty,time)
  print(Sys.time()-time)
  return(results_sim_real)
  
  
}



results_sim_real <- datasetReal()

#output_results_final1 <- cbind(output_results_final1,N_days=c(),cumul_days=c(),current_obs_S=c(), current_obs_E=c(),
#                               cumul_obs_S=c(),cumul_obs_E=c())


real_results <- results_sim_real$dataframe
View(real_results)
write.csv(real_results,"D:/Stage DataMa/Final datasets/Real AB tests results final.csv",row.names = FALSE)

library(dplyr)

data <- read.csv("D:/Stage DataMa/Final datasets/Real AB tests results final.csv")
data <- data %>% filter (test_id=="AB Test Split | New Funnel CP")
View(data)

real_boxplot <- results_sim_real$boxplot
View(real_boxplot)
write.csv(real_boxplot,"D:/Stage DataMa/Final datasets/Real AB tests boxplot final.csv",row.names = FALSE)


#results1 <- as.data.frame(read.csv("D:/Stage DataMa/Final datasets/Final dataset test total 1.csv"))
#results1 <- rbind (results1,real_results)
#write.csv (results1,"D:/Stage DataMa/Final datasets/Final dataset total 2.csv",row.names = FALSE)


#boxplot1 <- as.data.frame(read.csv("D:/Stage DataMa/Final datasets/Final boxplot test.csv"))
#boxplot1 <- rbind(boxplot1,real_boxplot)
#write.csv (boxplot1,"D:/Stage DataMa/Final datasets/Final boxplot total 2.csv",row.names = FALSE)
