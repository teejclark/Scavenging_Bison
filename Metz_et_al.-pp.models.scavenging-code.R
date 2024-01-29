#Code for three wolf-prey models in Metz et al. "Scavenging stabilizes predator-prey dynamics: insights from twenty-five years of wolf-elk-bison relationships in northern Yellowstone"

#--------------------------------------------------------------------------------#
#Set working directory 
#--------------------------------------------------------------------------------#

setwd() #location where have 'data_model.variable.values_el.ltr.csv'
getwd()

#--------------------------------------------------------------------------------#
#Bring in data 
#--------------------------------------------------------------------------------#

#####-----#####-----#####-----#####-----
#estimates for parameter values
###---###---###---###---

#This can be pulled in from available .csv file
df_params <- read.csv('data_model.variable.values_el.ltr.csv', header=TRUE,sep=",", na.strings="NA")

#--------------------------------------------------------------------------------#
#Create functions for p-p models
#--------------------------------------------------------------------------------#

#####-----#####-----#####-----#####-----
#create dataframe to hold predicted populations
###---###---###---###---

#Provide length of predictions and initial abundances for wolves, elk, and bison
create.all.popn.df = function(time.length, wolf.start, elk.start, bison.start){
  
  pop.df <- data.frame(time = seq(0, time.length-1, by = 1),
                       
                       #wolf-elk population - killing only
                       elk_st.1_max.popn = rep(NA, time.length), elk_st.1_kill = rep(NA,time.length), elk_st.1_end.popn = rep(NA,time.length), elk_st.1_pred.rate = rep(NA,time.length), elk_st.1_r.popn = rep(NA,time.length),
                       wolf_st.1_max.popn = rep(NA, time.length), wolf_st.1_death = rep(NA,time.length), wolf_st.1_end.popn = rep(NA,time.length), wolf_st.1_r.popn = rep(NA,time.length),
                       
                       #wolf-elk-bison population - killing only
                       elk_st.2_max.popn = rep(NA, time.length), elk_st.2_kill = rep(NA,time.length), elk_st.2_end.popn = rep(NA,time.length), elk_st.2_pred.rate = rep(NA,time.length), elk_st.2_r.popn = rep(NA,time.length),
                       bison_st.2_max.popn = rep(NA, time.length), bison_st.2_kill = rep(NA,time.length), bison_st.2_end.popn = rep(NA,time.length), bison_st.2_pred.rate = rep(NA,time.length), bison_st.2_r.popn = rep(NA,time.length),
                       wolf_st.2_max.popn = rep(NA, time.length), wolf_st.2_death = rep(NA,time.length), wolf_st.2_end.popn = rep(NA,time.length), wolf_st.2_r.popn = rep(NA,time.length),
                       
                       #wolf-elk-bison population - killing AND scavenging
                       elk_st.3_max.popn = rep(NA, time.length), elk_st.3_kill = rep(NA,time.length), elk_st.3_end.popn = rep(NA,time.length), elk_st.3_pred.rate = rep(NA,time.length), elk_st.3_r.popn = rep(NA,time.length),
                       bison_st.3_max.popn = rep(NA, time.length), bison_st.3_kill = rep(NA,time.length), bison_st.3_end.popn = rep(NA,time.length), bison_st.3_pred.rate = rep(NA,time.length), bison_st.3_r.popn = rep(NA,time.length),
                       wolf_st.3_max.popn = rep(NA, time.length), wolf_st.3_death = rep(NA,time.length), wolf_st.3_end.popn = rep(NA,time.length),  wolf_st.3_r.popn = rep(NA,time.length))
  
  #Fill initial populations
  pop.df$elk_st.1_end.popn[1] <- elk.start
  pop.df$wolf_st.1_end.popn[1] <- wolf.start
  
  pop.df$elk_st.2_end.popn[1] <- elk.start
  pop.df$bison_st.2_end.popn[1] <- bison.start
  pop.df$wolf_st.2_end.popn[1] <- wolf.start
  
  pop.df$elk_st.3_end.popn[1] <- elk.start
  pop.df$bison_st.3_end.popn[1] <- bison.start
  pop.df$wolf_st.3_end.popn[1] <- wolf.start
  
  print(head(pop.df))
  
  print(tail(pop.df))
  
  return(pop.df)
}

#####-----#####-----#####-----#####-----
#function that predicts the populations for all 3 food webs (W-E, W-E-B, W-E-B-S)
###---###---###---###---

#yy = df that holds results, xx = list of parameter values
predict.pp.dynamics_all.food.webs = function(yy, xx){
  
  for(t in 2:(nrow(yy))){ #for each row
    
    ###
    #wolf-elk
    ###
    
    #Elk population equations
    yy$elk_st.1_max.popn[t] <- yy$elk_st.1_end.popn[t-1] + ((xx$r_elk[[1]] * yy$elk_st.1_end.popn[t-1]) * (1 - (yy$elk_st.1_end.popn[t-1] / xx$K_elk[[1]])))
    
    yy$elk_st.1_kill[t] <- yy$wolf_st.1_end.popn[t-1] *
      ((xx$a.elk[[1]] * yy$elk_st.1_max.popn[t]) /
         (1 + (xx$a.elk[[1]] * xx$h.elk[[1]] * yy$elk_st.1_max.popn[t])))
    
    yy$elk_st.1_end.popn[t] <- yy$elk_st.1_max.popn[t] - yy$elk_st.1_kill[t]
    
    #Predation rate on elk
    yy$elk_st.1_pred.rate[t] <- yy$elk_st.1_kill[t] / yy$elk_st.1_max.popn[t]
    
    #Growth rate of elk
    yy$elk_st.1_r.popn[t] <- log(yy$elk_st.1_end.popn[t] / yy$elk_st.1_end.popn[t-1])
    
    #Wolf population equations
    yy$wolf_st.1_max.popn[t] <- yy$wolf_st.1_end.popn[t-1] *
      (((xx$b.wolf[[1]] * xx$a.elk[[1]] * yy$elk_st.1_max.popn[t]) /
          (1 + (xx$a.elk[[1]] * xx$h.elk[[1]] * yy$elk_st.1_max.popn[t]))))
    
    yy$wolf_st.1_death[t] <- yy$wolf_st.1_end.popn[t-1] * (xx$d.wolf[[1]] * yy$wolf_st.1_end.popn[t-1])
    
    yy$wolf_st.1_end.popn[t] <- yy$wolf_st.1_max.popn[t] - yy$wolf_st.1_death[t]
    
    #Growth rate of wolves
    yy$wolf_st.1_r.popn[t] <- log(yy$wolf_st.1_end.popn[t] / yy$wolf_st.1_end.popn[t-1])
    
    ###
    #wolf-elk-bison (killing only)
    ###
    
    #Initial ungulate growth equations
    
    yy$elk_st.2_max.popn[t] <- yy$elk_st.2_end.popn[t-1] + ((xx$r_elk[[1]] * yy$elk_st.2_end.popn[t-1]) * (1 - (yy$elk_st.2_end.popn[t-1] / xx$K_elk[[1]])))
    
    yy$bison_st.2_max.popn[t] <- yy$bison_st.2_end.popn[t-1] + ((xx$r_bison[[1]] * yy$bison_st.2_end.popn[t-1]) * (1 - (yy$bison_st.2_end.popn[t-1] / xx$M_bison[[1]])))
    
    #Subsequent elk population equations
    
    yy$elk_st.2_kill[t] <- yy$wolf_st.2_end.popn[t-1] *
      ((xx$a.elk[[1]] * yy$elk_st.2_max.popn[t]) /
         (1 +
            #elk part (type 2)
            (xx$a.elk[[1]] * xx$h.elk[[1]] * yy$elk_st.2_max.popn[t]) +
            #bison part (type 2, with size parameter)
            xx$a.bison[[1]] * xx$h.bison[[1]] * xx$wkb_s.param[[1]] * yy$bison_st.2_max.popn[t]))
    
    yy$elk_st.2_end.popn[t] <- yy$elk_st.2_max.popn[t] - yy$elk_st.2_kill[t]
    
    #Predation rate on elk
    yy$elk_st.2_pred.rate[t] <- yy$elk_st.2_kill[t] / yy$elk_st.2_max.popn[t]
    
    #Growth rate of elk
    yy$elk_st.2_r.popn[t] <- log(yy$elk_st.2_end.popn[t] / yy$elk_st.2_end.popn[t-1])
    
    #Subsequent bison population equations
    
    yy$bison_st.2_kill[t] <- yy$wolf_st.2_end.popn[t-1] *
      ((xx$a.bison[[1]] * yy$bison_st.2_max.popn[t]) /
         (1 +
            #elk part (type 2)
            (xx$a.elk[[1]] * xx$h.elk[[1]] * yy$elk_st.2_max.popn[t]) +
            #bison part (type 2, with size parameter)
            xx$a.bison[[1]] * xx$h.bison[[1]] * xx$wkb_s.param[[1]] * yy$bison_st.2_max.popn[t]))
    
    yy$bison_st.2_end.popn[t] <- yy$bison_st.2_max.popn[t] - yy$bison_st.2_kill[t]
    
    #Predation rate on bison
    yy$bison_st.2_pred.rate[t] <- yy$bison_st.2_kill[t] / yy$bison_st.2_max.popn[t]
    
    #Growth rate of bison
    yy$bison_st.2_r.popn[t] <- log(yy$bison_st.2_end.popn[t] / yy$bison_st.2_end.popn[t-1])
    
    #Wolf population equations
    yy$wolf_st.2_max.popn[t] <- yy$wolf_st.2_end.popn[t-1] *
      (xx$b.wolf[[1]] * (xx$a.elk[[1]] * yy$elk_st.2_max.popn[t] + xx$a.bison[[1]] * xx$wkb_s.param[[1]] * yy$bison_st.2_max.popn[t])) /
      (1 + (xx$a.elk[[1]] * xx$h.elk[[1]] * yy$elk_st.2_max.popn[t] +
              xx$a.bison[[1]] * xx$h.bison[[1]] * xx$wkb_s.param[[1]] * yy$bison_st.2_max.popn[t]))
    
    yy$wolf_st.2_death[t] <- yy$wolf_st.2_end.popn[t-1] * (xx$d.wolf[[1]] * yy$wolf_st.2_end.popn[t-1])
    
    yy$wolf_st.2_end.popn[t] <- yy$wolf_st.2_max.popn[t] - yy$wolf_st.2_death[t]
    
    #Growth rate of wolves
    yy$wolf_st.2_r.popn[t] <- log(yy$wolf_st.2_end.popn[t] / yy$wolf_st.2_end.popn[t-1])
    
    ###
    #wolf-elk-bison (including scavenging)
    ###
    
    #Initial ungulate growth equations
    yy$elk_st.3_max.popn[t] <- yy$elk_st.3_end.popn[t-1] + ((xx$r_elk[[1]] * yy$elk_st.3_end.popn[t-1]) * (1 - (yy$elk_st.3_end.popn[t-1] / xx$K_elk[[1]])))
    
    yy$bison_st.3_max.popn[t] <- yy$bison_st.3_end.popn[t-1] + ((xx$r_bison[[1]] * yy$bison_st.3_end.popn[t-1]) * (1 - (yy$bison_st.3_end.popn[t-1] / xx$M_bison[[1]])))
    
    #Subsequent elk population equations
    
    yy$elk_st.3_kill[t] <- yy$wolf_st.3_end.popn[t-1] *
      ((xx$a.elk[[1]] * yy$elk_st.3_max.popn[t]) /
         (1 +
            #elk part (type 2)
            (xx$a.elk[[1]] * xx$h.elk[[1]] * yy$elk_st.3_max.popn[t]) +
            #bison part (type 2, with size parameter)
            xx$a.bison[[1]] * xx$h.bison[[1]] * xx$wkb_s.param[[1]] * yy$bison_st.3_max.popn[t] +
            #bison scavenging part (type 2, with size parameter)
            xx$a.scav.bison[[1]] * xx$h.scav.bison[[1]] * xx$wsb_s.param[[1]] * yy$bison_st.3_max.popn[t]))
    
    yy$elk_st.3_end.popn[t] <- yy$elk_st.3_max.popn[t] - yy$elk_st.3_kill[t]
    
    #Predation rate on elk
    yy$elk_st.3_pred.rate[t] <- yy$elk_st.3_kill[t] / yy$elk_st.3_max.popn[t]
    
    #Growth rate of elk
    yy$elk_st.3_r.popn[t] <- log(yy$elk_st.3_end.popn[t] / yy$elk_st.3_end.popn[t-1])
    
    #Subsequent bison population equations
    
    yy$bison_st.3_kill[t] <- yy$wolf_st.3_end.popn[t-1] *
      ((xx$a.bison[[1]] * yy$bison_st.3_max.popn[t]) /
         (1 +
            #elk part (type 2)
            (xx$a.elk[[1]] * xx$h.elk[[1]] * yy$elk_st.3_max.popn[t]) +
            #bison part (type 2, with size parameter)
            xx$a.bison[[1]] * xx$h.bison[[1]] * xx$wkb_s.param[[1]] * yy$bison_st.3_max.popn[t] +
            #bison scavenging part (type 2, with size parameter)
            xx$a.scav.bison[[1]] * xx$h.scav.bison[[1]] * xx$wsb_s.param[[1]] * yy$bison_st.3_max.popn[t]))
    
    yy$bison_st.3_end.popn[t] <- yy$bison_st.3_max.popn[t] - yy$bison_st.3_kill[t]
    
    #Predation rate on bison
    yy$bison_st.3_pred.rate[t] <- yy$bison_st.3_kill[t] / yy$bison_st.3_max.popn[t]
    
    #Growth rate of bison
    yy$bison_st.3_r.popn[t] <- log(yy$bison_st.3_end.popn[t] / yy$bison_st.3_end.popn[t-1])
    
    #Wolf population equations
    yy$wolf_st.3_max.popn[t] <- yy$wolf_st.3_end.popn[t-1] *
      (xx$b.wolf[[1]] * (xx$a.elk[[1]] * yy$elk_st.3_max.popn[t] + xx$a.bison[[1]] * xx$wkb_s.param[[1]] * yy$bison_st.3_max.popn[t] + xx$a.scav.bison[[1]] * xx$wsb_s.param[[1]] * yy$bison_st.3_max.popn[t])) /
      (1 + (xx$a.elk[[1]] * xx$h.elk[[1]] * yy$elk_st.3_max.popn[t] +
              xx$a.bison[[1]] * xx$h.bison[[1]] * xx$wkb_s.param[[1]] * yy$bison_st.3_max.popn[t] +
              xx$a.scav.bison[[1]] * xx$h.scav.bison[[1]] * xx$wsb_s.param[[1]] * yy$bison_st.3_max.popn[t]))
    
    yy$wolf_st.3_death[t] <- yy$wolf_st.3_end.popn[t-1] * (xx$d.wolf[[1]] * yy$wolf_st.3_end.popn[t-1])
    
    yy$wolf_st.3_end.popn[t] <- yy$wolf_st.3_max.popn[t] - yy$wolf_st.3_death[t]
    
    #Growth rate of wolves
    yy$wolf_st.3_r.popn[t] <- log(yy$wolf_st.3_end.popn[t] / yy$wolf_st.3_end.popn[t-1])
    
  }
  
  return(yy)
  
}

#--------------------------------------------------------------------------------#
#Predict populations
#--------------------------------------------------------------------------------#

#####-----#####-----#####-----#####-----
#Create df to hold results
###---###---###---###---

init.popn.df <- create.all.popn.df(time.length = 501, wolf.start = 70, elk.start = 10600, bison.start = 700)

#####-----#####-----#####-----#####-----
#Create list of param values
###---###---###---###---

init.param.values.list <- list(K_elk = df_params %>% dplyr::filter(variable == "K_elk") %>% dplyr::pull(value),
                               r_elk = df_params %>% dplyr::filter(variable == "r_elk") %>% dplyr::pull(value),
                               M_bison = df_params %>% dplyr::filter(variable == "M_bison") %>% dplyr::pull(value),
                               r_bison = df_params %>% dplyr::filter(variable == "r_bison") %>% dplyr::pull(value),
                               wkb_s.param = df_params %>% dplyr::filter(variable == "wkb_s.param") %>% dplyr::pull(value),
                               wsb_s.param = df_params %>% dplyr::filter(variable == "wsb_s.param") %>% dplyr::pull(value),
                               b.wolf = df_params %>% dplyr::filter(variable == "b.wolf") %>% dplyr::pull(value),
                               d.wolf = df_params %>% dplyr::filter(variable == "d.wolf") %>% dplyr::pull(value),
                               #Elk killing functional response
                               a.elk = df_params %>% dplyr::filter(variable == "a.elk") %>% dplyr::pull(value),
                               h.elk = df_params %>% dplyr::filter(variable == "h.elk") %>% dplyr::pull(value),
                               #Bison killing functional response
                               a.bison = df_params %>% dplyr::filter(variable == "a.bison") %>% dplyr::pull(value),
                               h.bison = df_params %>% dplyr::filter(variable == "h.bison") %>% dplyr::pull(value),
                               #Bison scavenging functional response
                               a.scav.bison = df_params %>% dplyr::filter(variable == "a.scav.bison") %>% dplyr::pull(value),
                               h.scav.bison = df_params %>% dplyr::filter(variable == "h.scav.bison") %>% dplyr::pull(value)) 

#####-----#####-----#####-----#####-----
#Predict the populations
###---###---###---###---

#columns with:
#1) 'species_st.1' in their name represent predictions from W-E model
#2) 'species_st.2' in their name represent predictions from W-E-B model
#3) 'species_st.3' in their name represent predictions from W-E-B-S model

df.populations <- predict.pp.dynamics_all.food.webs(yy = init.popn.df, xx = init.param.values.list)

head(df.populations)
tail(df.populations)


