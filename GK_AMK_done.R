library(tidyverse)
library(outliers)
library(drc)


setwd("C:/Users/15183/Documents/Ziva data/Adarza")
concentration <- c(10000,2000,400,80,16,3.2,0.64, 0, 10000,2000,400,80,16,3.2,0.64, 0, 10000,2000,400,80,16,3.2,0.64, 0)
concentration_2 <- c(10000,2000,400,80,16,3.2,0.64, 0)
conc_3 <- rep(concentration_2,times=c(12,12,12,12,12,12,12,12))
data <- read.csv(file = "segresults_BY_IMAGE_SEGMENT_389.csv") %>%
  filter(Exposure.MS == 9999)

outlier_deletion_mse <- function(df){
  processed_df <- df %>%
    group_by(Well.Name, Protein.Name) %>%
    filter(((0.008 > (((Raw.Thickness-mean(Raw.Thickness))^2)/(mean(Raw.Thickness)^2))))) %>%
    return(processed_df)
}
mse <- outlier_deletion_mse(data)

mse_data <- write.csv(mse, "mse_data.csv", row.names = FALSE)

#detects outliers outside 90% confidence interval and deletes them
outlier_deletion <- function(df){
  processed_df <- df %>%
    group_by(Well.Name, Protein.Name) %>%
    filter(Raw.Thickness > mean(Raw.Thickness) + qnorm(0.025) * (sd(Raw.Thickness)/sqrt(9))) %>%
    filter(Raw.Thickness < mean(Raw.Thickness) + qnorm(0.975) * (sd(Raw.Thickness)/sqrt(9))) %>%
    return(processed_df)
}
No_outliers <- outlier_deletion(data)



outlier_deletion <- function(df){
  processed_df <- df %>%
    group_by(Well.Name, Protein.Name) %>%
    filter(Raw.Thickness < quantile(Raw.Thickness, 0.75)+1.5*(IQR(Raw.Thickness))) %>%
    filter(Raw.Thickness > quantile(Raw.Thickness, 0.25)-1.5*(IQR(Raw.Thickness))) %>%
    #filter(Raw.Thickness > (mean(Raw.Thickness) - qnorm(0.05) * (sd(Raw.Thickness)/sqrt(9)))) %>%
    #filter(Raw.Thickness < (mean(Raw.Thickness) + qnorm(0.95) * (sd(Raw.Thickness)/sqrt(9)))) %>%
    return(processed_df)
}


qnorm(0.95,8)
write.csv(No_outliers,"No_outliers.csv", row.names=FALSE)

#similar to last function, but takes the average of specific probe spot but can be averaged for chips starting with the same letter
average_no_outliers_all <- function(df, chip, protein, concentration){
  df_info <- df %>% 
    filter(Well.Name %in% chip) %>%
    group_by(substr(Well.Name, 1, 1)) %>%
    filter(Protein.Name == protein) %>%
    summarize(mean_thickness = mean(Raw.Thickness),
              std_thickness = sd(Raw.Thickness),
              count = n(),
              sem_thickness = std_thickness/(sqrt(count)),
              .groups = 'drop') %>%
    as.data.frame() %>%
    mutate(dose = concentration) %>%
    return(df_info)
}
#Averages for each probe for all chips with outliers removed

all_averages <- function(df){
  df_info <- df %>%
    group_by(Protein.Name,Well.Name) %>%
    summarize(mean_thickness = mean(Raw.Thickness),
              std_thickness = sd(Raw.Thickness),
              count = n(),
              sem_thickness = std_thickness/(sqrt(count)),
              .groups = 'drop') %>%
    as.data.frame() %>%
    return(df_info)
}
all_Data <- all_averages(No_outliers)

write.csv(all_Data,"average_thickness_no_outliers.csv", row.names=FALSE)

#Group standard curve strips 2,3,and 4 together to prepare for 5PL fit

standard_averages <- function(df,conc){
  dfinfo <- df %>%
    filter(Well.Name %in% c("A2","B2","C2","D2","E2","F2","G2","H2", "A3","B3","C3","D3","E3","F3","G3","H3", "A4","B4","C4","D4","E4","F4","G4","H4"))%>%
    group_by(Protein.Name,substr(Well.Name,1,1))%>%
    summarize(standards_mean_thickness = mean(mean_thickness),
              .groups = 'drop') %>%
    as.data.frame() %>%
    mutate(dose = rep(conc,15))%>%
    return(df_info)
  
}
standards <- standard_averages(all_Data,concentration_2)
write.csv(standards,"standards_averages.csv",row.names=FALSE)


proteins <- subset(all_Data,Well.Name=="A2")$Protein.Name



backfit_blank <- function(df_train, df_test, conc){
  train_model <- drm(standards_mean_thickness ~ dose, data = subset(df_train,Protein.Name==proteins[1]), fct=LL.5())
  plot(train_model)
  dev.copy(png,'blank.png')
  dev.off()
  coef <- train_model$coefficients
  capture.output(coef, file = "blank_coef.csv") 
  well <- subset(df_test,Protein.Name==proteins[1])$Well.Name
  x_test <- coef[4]*(((((coef[3]-coef[2])/(subset(df_test,Protein.Name==proteins[1])$mean_thickness-coef[2]))^(1/coef[5]))-1)^(1/coef[1]))
  difference <- (x_test/conc[-8])
  mylist <- list("conc"=well,"x-test"=x_test, "difference"= difference)
  return(mylist)
}


blank_fit <- as.data.frame(backfit_blank(standards, all_Data, conc_3),row.names=NULL,optional=TRUE)
write.csv(blank_fit,"blank_fit.csv",row.names=FALSE)

backfit_GM_SCF <- function(df_train, df_test, conc){
  train_model <- drm(standards_mean_thickness ~ dose, data = subset(df_train,Protein.Name==proteins[2]), fct=LL.5())
  plot(train_model)
  dev.copy(png,'GM_SCF.png')
  dev.off()
  coef <- train_model$coefficients
  capture.output(coef, file = "GM_SCF_coef.csv") 
  well <- subset(df_test,Protein.Name==proteins[1])$Well.Name
  x_test <- coef[4]*(((((coef[3]-coef[2])/(subset(df_test,Protein.Name==proteins[2])$mean_thickness-coef[2]))^(1/coef[5]))-1)^(1/coef[1]))
  difference <- (x_test/conc[-8])
  mylist <- list("well"=well,"x-test"=x_test, "difference"= difference)
  return(mylist)
}


GM_SCF_fit <- as.data.frame(backfit_GM_SCF(standards, all_Data, conc_3),row.names=NULL,optional=TRUE)
write.csv(GM_SCF_fit,"GM_SCF_fit.csv",row.names=FALSE)

backfit_IFNg <- function(df_train, df_test, conc){
  train_model <- drm(standards_mean_thickness ~ dose, data = subset(df_train,Protein.Name==proteins[3]), fct=LL.5())
  plot(train_model)
  dev.copy(png,'IFNg.png')
  dev.off()
  coef <- train_model$coefficients
  capture.output(coef, file = "IFNg_coef.csv") 
  well <- subset(df_test,Protein.Name==proteins[1])$Well.Name
  x_test <- coef[4]*(((((coef[3]-coef[2])/(subset(df_test,Protein.Name==proteins[3])$mean_thickness-coef[2]))^(1/coef[5]))-1)^(1/coef[1]))
  difference <- (x_test/conc)
  mylist <- list("well"=well,"x-test"=x_test,"difference"= difference )
  return(mylist)
}


IFNg_fit <- as.data.frame(backfit_IFNg(standards, all_Data, conc_3),row.names=NULL,optional=TRUE)
write.csv(IFNg_fit,"IFNg_fit.csv",row.names=FALSE)

backfit_IL10 <- function(df_train, df_test, conc){
  train_model <- drm(standards_mean_thickness ~ dose, data = subset(df_train,Protein.Name==proteins[4]), fct=LL.5())
  plot(train_model)
  dev.copy(png,'IL10.png')
  dev.off()
  coef <- train_model$coefficients
  capture.output(coef, file = "IL10_coef.csv") 
  well <- subset(df_test,Protein.Name==proteins[1])$Well.Name
  x_test <- coef[4]*(((((coef[3]-coef[2])/(subset(df_test,Protein.Name==proteins[4])$mean_thickness-coef[2]))^(1/coef[5]))-1)^(1/coef[1]))
  difference <- (x_test/conc)
  mylist <- list("well"=well,"x-test"=x_test,"difference"= difference )
  return(mylist)
}


IL10_fit <- as.data.frame(backfit_IL10(standards, all_Data, conc_3),row.names=NULL,optional=TRUE)
write.csv(IL10_fit,"IL10_fit.csv",row.names=FALSE)

backfit_IL17A <- function(df_train, df_test, conc){
  train_model <- drm(standards_mean_thickness ~ dose, data = subset(df_train,Protein.Name==proteins[5]), fct=LL.5())
  plot(train_model)
  dev.copy(png,'IL17A.png')
  dev.off()
  coef <- train_model$coefficients
  capture.output(coef, file = "IL17A_coef.csv") 
  well <- subset(df_test,Protein.Name==proteins[1])$Well.Name
  x_test <- coef[4]*(((((coef[3]-coef[2])/(subset(df_test,Protein.Name==proteins[5])$mean_thickness-coef[2]))^(1/coef[5]))-1)^(1/coef[1]))
  difference <- (x_test/conc)
  mylist <- list("well"=well,"x-test"=x_test,"difference"= difference )
  return(mylist)
}


IL17A_fit <- as.data.frame(backfit_IL17A(standards, all_Data, conc_3),row.names=NULL,optional=TRUE)
write.csv(IL17A_fit,"IL17A_fit.csv",row.names=FALSE)


backfit_IL2 <- function(df_train, df_test, conc){
  train_model <- drm(standards_mean_thickness ~ dose, data = subset(df_train,Protein.Name==proteins[6]), fct=LL.5())
  plot(train_model)
  dev.copy(png,'IL2.png')
  dev.off()
  coef <- train_model$coefficients
  capture.output(coef, file = "IL2_coef.csv") 
  well <- subset(df_test,Protein.Name==proteins[1])$Well.Name
  x_test <- coef[4]*(((((coef[3]-coef[2])/(subset(df_test,Protein.Name==proteins[6])$mean_thickness-coef[2]))^(1/coef[5]))-1)^(1/coef[1]))
  difference <- (x_test/conc)
  mylist <- list("well"=well,"x-test"=x_test,"difference"= difference )
  return(mylist)
}


IL2_fit <- as.data.frame(backfit_IL2(standards, all_Data, conc_3),row.names=NULL,optional=TRUE)
write.csv(IL2_fit,"IL2_fit.csv",row.names=FALSE)

backfit_IL3 <- function(df_train, df_test, conc){
  train_model <- drm(standards_mean_thickness ~ dose, data = subset(df_train,Protein.Name==proteins[7]), fct=LL.5())
  plot(train_model)
  dev.copy(png,'IL3.png')
  dev.off()
  coef <- train_model$coefficients
  capture.output(coef, file = "IL3_coef.csv") 
  well <- subset(df_test,Protein.Name==proteins[1])$Well.Name
  x_test <- coef[4]*(((((coef[3]-coef[2])/(subset(df_test,Protein.Name==proteins[7])$mean_thickness-coef[2]))^(1/coef[5]))-1)^(1/coef[1]))
  difference <- (x_test/conc)
  mylist <- list("well"=well,"x-test"=x_test,"difference"= difference )
  return(mylist)
}


IL3_fit <- as.data.frame(backfit_IL3(standards, all_Data, conc_3),row.names=NULL,optional=TRUE)
write.csv(IL3_fit,"IL3_fit.csv",row.names=FALSE)

backfit_IL33 <- function(df_train, df_test, conc){
  train_model <- drm(standards_mean_thickness ~ dose, data = subset(df_train,Protein.Name==proteins[8]), fct=LL.5())
  plot(train_model)
  dev.copy(png,'IL33.png')
  dev.off()
  coef <- train_model$coefficients
  capture.output(coef, file = "IL33_coef.csv") 
  well <- subset(df_test,Protein.Name==proteins[1])$Well.Name
  x_test <- coef[4]*(((((coef[3]-coef[2])/(subset(df_test,Protein.Name==proteins[8])$mean_thickness-coef[2]))^(1/coef[5]))-1)^(1/coef[1]))
  difference <- (x_test/conc)
  mylist <- list("well"=well,"x-test"=x_test,"difference"= difference )
  return(mylist)
}


IL33_fit <- as.data.frame(backfit_IL33(standards, all_Data, conc_3),row.names=NULL,optional=TRUE)
write.csv(IL33_fit,"IL33_fit.csv",row.names=FALSE)

backfit_IL4 <- function(df_train, df_test, conc){
  train_model <- drm(standards_mean_thickness ~ dose, data = subset(df_train,Protein.Name==proteins[9]), fct=LL.5())
  plot(train_model)
  dev.copy(png,'IL4.png')
  dev.off()
  coef <- train_model$coefficients
  capture.output(coef, file = "IL4_coef.csv") 
  well <- subset(df_test,Protein.Name==proteins[1])$Well.Name
  x_test <- coef[4]*(((((coef[3]-coef[2])/(subset(df_test,Protein.Name==proteins[9])$mean_thickness-coef[2]))^(1/coef[5]))-1)^(1/coef[1]))
  difference <- (x_test/conc)
  mylist <- list("well"=well,"x-test"=x_test,"difference"= difference )
  return(mylist)
}


IL4_fit <- as.data.frame(backfit_IL4(standards, all_Data, conc_3),row.names=NULL,optional=TRUE)
write.csv(IL4_fit,"IL4_fit.csv",row.names=FALSE)

backfit_IL5 <- function(df_train, df_test, conc){
  train_model <- drm(standards_mean_thickness ~ dose, data = subset(df_train,Protein.Name==proteins[10]), fct=LL.5())
  plot(train_model)
  dev.copy(png,'IL5.png')
  dev.off()
  coef <- train_model$coefficients
  capture.output(coef, file = "IL5_coef.csv") 
  well <- subset(df_test,Protein.Name==proteins[1])$Well.Name
  x_test <- coef[4]*(((((coef[3]-coef[2])/(subset(df_test,Protein.Name==proteins[10])$mean_thickness-coef[2]))^(1/coef[5]))-1)^(1/coef[1]))
  difference <- (x_test/conc)
  mylist <- list("well"=well,"x-test"=x_test,"difference"= difference )
  return(mylist)
}


IL5_fit <- as.data.frame(backfit_IL5(standards, all_Data, conc_3),row.names=NULL,optional=TRUE)
write.csv(IL5_fit,"IL5_fit.csv",row.names=FALSE)

backfit_IL6 <- function(df_train, df_test, conc){
  train_model <- drm(standards_mean_thickness ~ dose, data = subset(df_train,Protein.Name==proteins[11]), fct=LL.5())
  plot(train_model)
  dev.copy(png,'IL6.png')
  dev.off()
  coef <- train_model$coefficients
  capture.output(coef, file = "IL6_coef.csv") 
  well <- subset(df_test,Protein.Name==proteins[1])$Well.Name
  x_test <- coef[4]*(((((coef[3]-coef[2])/(subset(df_test,Protein.Name==proteins[11])$mean_thickness-coef[2]))^(1/coef[5]))-1)^(1/coef[1]))
  difference <- (x_test/conc)
  mylist <- list("well"=well,"x-test"=x_test,"difference"= difference )
  return(mylist)
}


IL6_fit <- as.data.frame(backfit_IL6(standards, all_Data, conc_3),row.names=NULL,optional=TRUE)
write.csv(IL6_fit,"IL6_fit.csv",row.names=FALSE)


backfit_IL9 <- function(df_train, df_test, conc){
  train_model <- drm(standards_mean_thickness ~ dose, data = subset(df_train,Protein.Name==proteins[12]), fct=LL.5())
  plot(train_model)
  dev.copy(png,'IL9.png')
  dev.off()
  coef <- train_model$coefficients
  capture.output(coef, file = "IL9_coef.csv") 
  well <- subset(df_test,Protein.Name==proteins[1])$Well.Name
  x_test <- coef[4]*(((((coef[3]-coef[2])/(subset(df_test,Protein.Name==proteins[12])$mean_thickness-coef[2]))^(1/coef[5]))-1)^(1/coef[1]))
  difference <- (x_test/conc)
  mylist <- list("well"=well,"x-test"=x_test,"difference"= difference )
  return(mylist)
}


IL9_fit <- as.data.frame(backfit_IL9(standards, all_Data, conc_3),row.names=NULL,optional=TRUE)
write.csv(IL9_fit,"IL9_fit.csv",row.names=FALSE)

backfit_biotin <- function(df_train, df_test, conc){
  train_model <- drm(standards_mean_thickness ~ dose, data = subset(df_train,Protein.Name==proteins[13]), fct=LL.5())
  plot(train_model)
  dev.copy(png,'biotin.png')
  dev.off()
  coef <- train_model$coefficients
  capture.output(coef, file = "biotin_coef.csv") 
  well <- subset(df_test,Protein.Name==proteins[1])$Well.Name
  x_test <- coef[4]*(((((coef[3]-coef[2])/(subset(df_test,Protein.Name==proteins[13])$mean_thickness-coef[2]))^(1/coef[5]))-1)^(1/coef[1]))
  difference <- (x_test/conc)
  mylist <- list("well"=well,"x-test"=x_test,"difference"= difference )
  return(mylist)
}


biotin_fit <- as.data.frame(backfit_biotin(standards, all_Data, conc_3),row.names=NULL,optional=TRUE)
write.csv(biotin_fit,"biotin_fit.csv",row.names=FALSE)

backfit_HRP <- function(df_train, df_test, conc){
  train_model <- drm(standards_mean_thickness ~ dose, data = subset(df_train,Protein.Name==proteins[14]), fct=LL.5())
  plot(train_model)
  dev.copy(png,'HRP.png')
  dev.off()
  coef <- train_model$coefficients
  capture.output(coef, file = "HRP_coef.csv") 
  well <- subset(df_test,Protein.Name==proteins[1])$Well.Name
  x_test <- coef[4]*(((((coef[3]-coef[2])/(subset(df_test,Protein.Name==proteins[14])$mean_thickness-coef[2]))^(1/coef[5]))-1)^(1/coef[1]))
  difference <- (x_test/conc)
  mylist <- list("well"=well,"x-test"=x_test,"difference"= difference )
  return(mylist)
}


HRP_fit <- as.data.frame(backfit_HRP(standards, all_Data, conc_3),row.names=NULL,optional=TRUE)
write.csv(HRP_fit,"HRP_fit.csv",row.names=FALSE)

backfit_TNFa <- function(df_train, df_test, conc){
  train_model <- drm(standards_mean_thickness ~ dose, data = subset(df_train,Protein.Name==proteins[15]), fct=LL.5())
  plot(train_model)
  dev.copy(png,'TNFa.png')
  dev.off()
  coef <- train_model$coefficients
  capture.output(coef, file = "TNFa_coef.csv") 
  well <- subset(df_test,Protein.Name==proteins[1])$Well.Name
  x_test <- coef[4]*(((((coef[3]-coef[2])/(subset(df_test,Protein.Name==proteins[15])$mean_thickness-coef[2]))^(1/coef[5]))-1)^(1/coef[1]))
  difference <- (x_test/conc)
  mylist <- list("well"=well,"x-test"=x_test,"difference"= difference )
  return(mylist)
}


TNFa_fit <- as.data.frame(backfit_TNFa(standards, all_Data, conc_3),row.names=NULL,optional=TRUE)
write.csv(TNFa_fit,"TNFa_fit.csv",row.names=FALSE)