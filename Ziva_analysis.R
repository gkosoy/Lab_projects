library(tidyverse)
library(drc)
library(readxl)

setwd("C:/Users/15183/Documents/Ziva data/Feb_8_2021_141_Day0,28") #where your working files are
Ziva_data_Dec2 <- read.csv(file="segresults_BY_IMAGE_SEGMENT_419.csv") %>%  
  filter(Exposure_MS == 60)  #open Ziva data file, and choose 60 or 140ms exposure data
dilutions_1 <- read_xlsx("doses_141_day0.xlsx") #make excel file of your dilution/concentration information and read here
dilutions_2 <- read_xlsx("doses_141_day28.xlsx")

#Enter your data file (df), chips to analyze (well_names), protein to analyze (protein_name)
filtered_data <- function(df, well_names, protein_name){
  processed_df <- df %>% 
    filter(Well_Name %in% well_names) %>%
    filter(Protein_Name == protein_name) %>%
    
    return(processed_df)
  
}

#example, analyze SARS-CoV-2 RBD protein of day 0 and 28 serum chips. chip "A" is my blank chip. 1's and 2's are Day 0 replicates, 3's and 4's are Day 28 replicates
RBD_Day0 <- filtered_data(Ziva_data_Dec2, list("B1", "C1", "D1", "E1","F1","G1","H1","A1", "B2", "C2", "D2", "E2", "F2", "G2", "H2", "A2"), "SARS-CoV-2 RBD")
RBD_Day28 <- filtered_data(Ziva_data_Dec2, list("B3","C3","D3","E3","F3", "G3","H3","A3","B4", "C4","E4","F4","G4","H4","A4"), "SARS-CoV-2 RBD")
CytC_Day0 <- filtered_data(Ziva_data_Dec2, list("B1", "C1", "D1", "E1","F1","G1","H1","A1", "B2", "C2", "D2", "E2", "F2", "G2", "H2", "A2"), "cytC")
CytC_Day28 <- filtered_data(Ziva_data_Dec2, list("B3","C3","D3","E3","F3", "G3","H3","A3","B4", "C4","E4","F4","G4","H4","A4"), "cytC")
#Chip D4 had defects in running
#You will see that there are five RBD probes on each chip. CytC can be used to control for non-specific binding

# Now we will get mean, and SEM information to make our graphs
info_by_protein <- function(df, dilution) {
  df_info <- df %>% 
    group_by(Well_Name) %>%
    summarize(mean_thickness = mean(Raw_Thickness),
              std_thickness = sd(Raw_Thickness),
              count = n(),
              sem_thickness = std_thickness/(sqrt(count)),
              .groups = 'drop') %>%
    as.data.frame() %>%
    mutate(dose = dilution) %>%
    return(df_info)
}
graph_RBD_1 <- info_by_protein(RBD_Day0, dilutions_1$`SARS2 RBD`)  
graph_RBD_2 <- info_by_protein(RBD_Day28, dilutions_2$`SARS2 RBD`[-8])  #-8 is to excluse chip D4
info_cytC_1 <- info_by_protein(CytC_Day0, dilutions_1$`SARS2 RBD`)
info_cytC_2 <- info_by_protein(CytC_Day28, dilutions_2$`SARS2 RBD`[-8])
#Put your CytC corrected values in new column
graph_RBD_1[,7] <- graph_RBD_1$mean_thickness-info_cytC_1$mean_thickness
graph_RBD_1[,8] <- sqrt(((graph_RBD_1$std_thickness)^2)+((info_cytC_1$std_thickness)^2))
graph_RBD_2[,7] <- graph_RBD_2$mean_thickness-info_cytC_2$mean_thickness
graph_RBD_2[,8] <- sqrt(((graph_RBD_2$std_thickness)^2)+((info_cytC_2$std_thickness)^2))
# If points are on the left side of the curve you can correct those points 
graph_RBD_1$V7[1:2] <- graph_RBD_1$V7[1:2]*-1
graph_RBD_2$V7[1:2] <- graph_RBD_2$V7[1:2]*-1
graph_RBD_1$V7[7:16] <- graph_RBD_1$V7[7:16]*-1
graph_RBD_2$V7[12:15] <- graph_RBD_2$V7[12:15]*-1


#get parameters of a four parameter logistic
RBD.m1 <- drm(V7 ~ dose, data = graph_RBD_1, fct = LL.4(fixed = c(NA,NA,NA,NA)))
RBD.m2 <- drm(V7 ~ dose, data = graph_RBD_2, fct = LL.4(fixed = c(NA,NA,NA,NA)))
parametersRBD_1 <- RBD.m1$coefficients; b = round(parametersRBD_1[1], digits = 2); c = round(parametersRBD_1[2], digits = 2); d = round(parametersRBD_1[3], digits = 2); e = round(parametersRBD_1[4], digits = 2)
parametersRBD_2 <- RBD.m2$coefficients; b2 = round(parametersRBD_2[1], digits = 2); c2 = round(parametersRBD_2[2], digits = 2); d2 = round(parametersRBD_2[3], digits = 2); e2 = round(parametersRBD_2[4], digits = 2)

RBD.m1 <- drm(V7 ~ dose, data = graph_RBD_1, fct = LL.4(fixed = c(NA,NA,45,NA)))
RBD.m2 <- drm(V7 ~ dose, data = graph_RBD_2, fct = LL.4(fixed = c(NA,NA,45,NA)))
parametersRBD_1 <- RBD.m1$coefficients; b = round(parametersRBD_1[1], digits = 2); c = round(parametersRBD_1[2], 2); d = 45; e = round(parametersRBD_1[3], digits = 2)
parametersRBD_2 <- RBD.m2$coefficients; b2 = round(parametersRBD_2[1], digits = 2); c2 = round(parametersRBD_2[2], 2); d2 = 45; e2 = round(parametersRBD_2[3], digits = 2)

#graph CytC corrected mean thickness of serum dilutions against SARS-CoV-2 RBD for Day 0 and Day 28
ggplot(data = NULL) +
  ggtitle(paste("Sample 141")) + xlab("Concentration of IgG (nM)") + ylab("Thickness change expression(ring(A))") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title=element_text(size=16)) +
  geom_pointrange(size = 0.5, data = graph_RBD_1, 
                  mapping=aes(color = "Day 0", x=dilutions_1$`SARS2 RBD`, y = V7, ymin = V7-V8, ymax = V7+V8)) +
  geom_pointrange(size = 0.5, data=graph_RBD_2, 
                  mapping=aes(color = "Day 28", x=dilutions_2$`SARS2 RBD`[-8], y = V7, ymin = V7-V8, ymax = V7+V8)) +
  scale_x_log10() +
  theme(legend.title = element_blank()) + 
  theme(axis.title.x = element_text(colour = "royalblue4", size =15)) + theme(axis.title.y = element_text(colour = "royalblue4", size =15)) +
  theme(axis.text.x= element_text(colour = "grey1", size = 14)) + theme(axis.text.y = element_text(colour = "grey1", size = 14)) +
  geom_function(fun = function(x) c2+((d2-c2)/abs(1+((x/e2)^b2))), colour = "#00BFC4", size = 1) +
  geom_function(fun = function(x) c+((d-c)/abs(1+((x/e)^b))), colour = "#F8766D", size = 1) +
  annotate(colour = "#00BFC4", "text", size = 5, x=0.035, y = 12, label = paste("Hill-", -b2, "\n", "Min-", c2, "\n", "Max-", d2, "\n", "Inflection-", e2))+
  annotate(colour = "#F8766D","text", size = 5, x=0.035, y = 5, label = paste("Hill-", -b, "\n", "Min-", c, "\n", "Max-", d, "\n", "Inflection-", e ))

