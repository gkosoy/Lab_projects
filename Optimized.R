library(tidyverse)
library(outliers)
library(drc)


setwd("C:/Users/15183/Documents/Ziva data/Adarza")
concentration <- c(10000,2000,400,80,16,3.2,0.64, 0, 10000,2000,400,80,16,3.2,0.64, 0, 10000,2000,400,80,16,3.2,0.64, 0)
concentration_2 <- c(10000,2000,400,80,16,3.2,0.64, 0)
data <- read.csv(file = "segresults_BY_IMAGE_SEGMENT_390.csv") %>%
  filter(Exposure_MS == 9999)

#detects outliers outside 90% confidence interval and deletes them
outlier_deletion <- function(df){
  processed_df <- df %>%
    group_by(Well_Name, Protein_Name) %>%
      filter(Raw_Thickness > quantile(Raw_Thickness,0.05)) %>%
      filter(Raw_Thickness < quantile(Raw_Thickness,0.95)) %>%
    return(processed_df)
}
No_outliers <- outlier_deletion(data)
#takes the average of each 9 probe area without outliers per chip name
average_no_outliers <- function(df, chip, protein, concentration){
  df_info <- df %>% 
    filter(Well_Name %in% chip) %>%
    filter(Protein_Name == protein) %>%
    summarize(mean_thickness = mean(Raw_Thickness),
              std_thickness = sd(Raw_Thickness),
              count = n(),
              sem_thickness = std_thickness/(sqrt(count)),
              .groups = 'drop') %>%
    as.data.frame() %>%
    mutate(concentration = concentration) %>%
    return(df_info)
}
IL2_data <- average_no_outliers(df = No_outliers, chip = list("A2","B2","C2","D2","E2","F2","G2","H2"), protein = "IL-2", concentration = concentration_2)
#similar to last function, but takes the average of specific probe spot but can be averaged for chips starting with the same letter
average_no_outliers_all <- function(df, chip, protein, concentration){
  df_info <- df %>% 
    filter(Well_Name %in% chip) %>%
    group_by(substr(Well_Name, 1, 1)) %>%
    filter(Protein_Name == protein) %>%
    summarize(mean_thickness = mean(Raw_Thickness),
              std_thickness = sd(Raw_Thickness),
              count = n(),
              sem_thickness = std_thickness/(sqrt(count)),
              .groups = 'drop') %>%
    as.data.frame() %>%
    mutate(dose = concentration) %>%
    return(df_info)
}
IL2_train <- average_no_outliers_all(df = No_outliers, chip = list("A2","B2","C2","D2","E2","F2","G2","H2", "A3","B3","C3","D3","E3","F3","G3","H3", "A4","B4","C4","D4","E4","F4","G4","H4"), protein = "IL-2", concentration = concentration_2)
IL3_train <- average_no_outliers_all(df = No_outliers, chip = list("A2","B2","C2","D2","E2","F2","G2","H2", "A3","B3","C3","D3","E3","F3","G3","H3", "A4","B4","C4","D4","E4","F4","G4","H4"), protein = "IL-3", concentration = concentration_2)
IL17A_train <- average_no_outliers_all(df = No_outliers, chip = list("A2","B2","C2","D2","E2","F2","G2","H2", "A3","B3","C3","D3","E3","F3","G3","H3", "A4","B4","C4","D4","E4","F4","G4","H4"), protein = "IL-17A", concentration = concentration_2)
IL17A_test_5 <- average_no_outliers(df = No_outliers, chip = list("A5","B5","C5","D5","E5","F5","G5","H5"), protein = "IL-17A", concentration = concentration_2)
IL17A_test_6 <- average_no_outliers(df = No_outliers, chip = list("A6","B6","C6","D6","E6","F6","G6","H6"), protein = "IL-17A", concentration = concentration_2)
IL17A_test_7 <- average_no_outliers(df = No_outliers, chip = list("A7","B7","C7","D7","E7","F7","G7","H7"), protein = "IL-17A", concentration = concentration_2)

#fits standard curve, fits test strip, takes difference
backfit <- function(df_train, df_test, x){
  train_model <- drm(mean_thickness ~ dose, data = df_train, fct = LL.5())
  coef <- train_model$coefficients
  test_model <- drm(mean_thickness ~ concentration, data = df_test, fct = LL.5())
  coef_2 <- test_model$coefficients
  y_train <- coef[3]+((coef[1]-coef[3])/(abs(1+((x/coef[4])^coef[2]))^coef[5]))
  y_test <- coef_2[3]+((coef_2[1]-coef_2[3])/(abs(1+((x/coef_2[4])^coef_2[2]))^coef_2[5]))
  difference <- y_test/y_train
  return(difference)
}

backfit(IL17A_train, y, concentration_2)

backfit <- function(df_train, df_test, x){
  train_model <- drm(mean_thickness ~ dose, data = df_train, fct = LL.5())
  coef <- train_model$coefficients
  x_test <- coef[4]*(((((coef[1]-coef[3])/(y-coef[3]))^(1/coef[5]))-1)^(1/coef[2]))
  difference <- x_test/concentration_2
  return(difference)
}


