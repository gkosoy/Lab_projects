library(tidyverse)

# extract dates from sample info
#### start ####
extract_date <- function(string){
  string %>% 
    stringr::str_extract(., regex("\\d{1,2}[/]\\d{1,2}[/]\\d{2}"))
}

clean_date_digits <- function(string){
  string %>% 
    stringr::str_replace(., "[/]20$", "/2020") %>% 
    stringr::str_replace(., "[/]21$", "/2021")
}

extract_n_clean_dates <- function(string){
  string %>% 
    extract_date() %>% 
    clean_date_digits()
}

ex1 <- "asdf 12/2/23"
ex2 <- "asdf 12/20/23"
ex3 <- "asdf 1/20/23"

extract_date(ex3)
#### end ####

raw_covid_data_2 <- read.csv("vaccine sample info.csv")

# create clean data set of all info with days post vaccines
####start####
clean_covid_data_2 <- raw_covid_data %>% 
  janitor::clean_names() %>% 
  rename(x1_months_post_2nd_dose = x1_month_post_2nd_dose,
         x4_months_post_2nd_dose =  x4_months_post_2nd_dose2) %>% 
  mutate(predraw_date = extract_n_clean_dates(predraw),
         x_1_wk_post_first_dose_date = extract_n_clean_dates(x_1_wk_post_first_dose),
         x_2_wk_post_first_dose_date = extract_n_clean_dates(x_2_wk_post_first_dose),
         x_3_wk_post_first_dose_date = extract_n_clean_dates(x_3_wk_post_vaccine),
         x_1_wk_post_second_dose_date = extract_n_clean_dates(post_2nd_dose),
         x1_2_days_post_2nd_dose_date = extract_n_clean_dates(x1_2_days_post_2nd_dose),
         post_2nd_dose_date = extract_n_clean_dates(post_2nd_dose),
         
         x1_months_post_2nd_dose_date = extract_n_clean_dates(x1_months_post_2nd_dose),
         x2_months_post_2nd_dose_date = extract_n_clean_dates(x2_months_post_2nd_dose),
         x3_months_post_2nd_dose_date = extract_n_clean_dates(x3_months_post_2nd_dose),
         x4_months_post_2nd_dose_date = extract_n_clean_dates(x4_months_post_2nd_dose),
         x5_months_post_2nd_dose_date = extract_n_clean_dates(x5_months_post_2nd_dose),
         x6_months_post_2nd_dose_date = extract_n_clean_dates(x6_months_post_2nd_dose),
         x7_months_post_2nd_dose_date = extract_n_clean_dates(x7_months_post_2nd_dose),
         pre_booster_date = extract_n_clean_dates(pre_booster),
         x1_wk_post_booster_date = extract_n_clean_dates(x1_wk_post_booster)
  )  %>% 
  mutate_at(vars(contains("date")), 
            ~as.Date(., format = "%m/%d/%Y")) %>% 
  mutate(w1_since_first = difftime(x_1_wk_post_first_dose_date, vaccine_1st_dose_date, units = "days"),
         w2_since_first = difftime(x_2_wk_post_first_dose_date, vaccine_1st_dose_date, units = "days"),
         w3_since_first = difftime(x_3_wk_post_first_dose_date, vaccine_1st_dose_date, units = "days"),
         d1_from_second = difftime(x1_2_days_post_2nd_dose_date, vaccine_2nd_dose_date, units = "days"),
         w1_from_second = difftime(post_2nd_dose_date, vaccine_2nd_dose_date, units = "days"),
         m1_days_from_second = difftime(x1_months_post_2nd_dose_date, vaccine_2nd_dose_date, units = "days"),
         m2_days_from_second = difftime(x2_months_post_2nd_dose_date, vaccine_2nd_dose_date, units = "days"),
         m3_days_from_second = difftime(x3_months_post_2nd_dose_date, vaccine_2nd_dose_date, units = "days"),
         m4_days_from_second = difftime(x4_months_post_2nd_dose_date, vaccine_2nd_dose_date, units = "days"),
         m5_days_from_second = difftime(x5_months_post_2nd_dose_date, vaccine_2nd_dose_date, units = "days"),
         m6_days_from_second = difftime(x6_months_post_2nd_dose_date, vaccine_2nd_dose_date, units = "days"),
         m7_days_from_second = difftime(x7_months_post_2nd_dose_date, vaccine_2nd_dose_date, units = "days"),
         pre_booster_from_second = difftime(pre_booster_date, vaccine_2nd_dose_date, units = "days"),
         boost_from_second= difftime(x1_wk_post_booster_date, vaccine_2nd_dose_date, units = "days")
         
  ) %>%
  mutate(s_00 = stringr::str_extract_all(predraw, "^SN\\d{3}"),
         s_01 = stringr::str_extract_all(x_1_wk_post_first_dose, "^SN\\d{3}"),
         s_02 = stringr::str_extract_all(x_2_wk_post_first_dose, "^SN\\d{3}"),
         s_03 = stringr::str_extract_all(x_3_wk_post_vaccine, "^SN\\d{3}"),
         s_04 = stringr::str_extract_all(x1_2_days_post_2nd_dose, "^SN\\d{3}"),
         s_05 = stringr::str_extract_all(post_2nd_dose, "^SN\\d{3}"),
         s_06 = stringr::str_extract_all(x1_months_post_2nd_dose, "^SN\\d{3}"),
         s_07 = stringr::str_extract_all(x2_months_post_2nd_dose, "^SN\\d{3}"),
         s_08 = stringr::str_extract_all(x3_months_post_2nd_dose, "^SN\\d{3}"),
         s_09 = stringr::str_extract_all(x4_months_post_2nd_dose, "^SN\\d{3}"),
         s_10 = stringr::str_extract_all(x5_months_post_2nd_dose, "^SN\\d{3}"),
         s_11 = stringr::str_extract_all(x6_months_post_2nd_dose, "^SN\\d{3}"),
         s_12 = stringr::str_extract_all(x7_months_post_2nd_dose, "^SN\\d{3}"),
         s_13 = stringr::str_extract_all(pre_booster, "^SN\\d{3}"),
         s_14 = stringr::str_extract_all(x1_wk_post_booster, "^SN\\d{3}")
  ) %>%
  mutate(ID = stringr::str_extract_all(subject_id_for_vaccine_data_purposes, "\\d{2}")
  )
####end####

#old labeling regime to create final_df
####start####
df_SN_2 <- clean_covid_data_2 %>% 
  select(c(s_p, s_1wk, s_2wk, s_3wk, s_1_2_d, s_d, s_1mon, s_2mon, 
           s_3mon, s_4mon, s_5mon, s_6mon, s_7mon, s_p_b, s_b, ID)) %>%
  pivot_longer(cols = s_p:s_b, names_to = "time_draw") %>%
  as.data.frame() %>%
  na_if("character(0)")

df_with_dates_2<- clean_covid_data_2 %>%
  select(c(ID, w1_since_first, w2_since_first, w3_since_first, d1_from_second, w1_from_second,
           m1_days_from_second, m2_days_from_second, m3_days_from_second,m4_days_from_second, 
           m5_days_from_second, m6_days_from_second, m7_days_from_second, pre_booster_from_second, 
           boost_from_second)) %>%
  mutate(s_1wk = w1_since_first, s_2wk = w2_since_first, s_3wk = w3_since_first, s_1_2_d = d1_from_second,
         s_d = w1_from_second,s_1mon = m1_days_from_second, s_2mon = m2_days_from_second, 
         s_3mon = m3_days_from_second, s_4mon = m4_days_from_second, s_5mon = m5_days_from_second, 
         s_6mon = m6_days_from_second, s_7mon = m7_days_from_second, s_p_b = pre_booster_from_second,
         s_b = boost_from_second) %>%
  pivot_longer(cols = s_1wk:s_b, names_to = "time_draw") %>%
  select(!c(w1_since_first, w2_since_first, w3_since_first, d1_from_second, w1_from_second,
            m1_days_from_second, m2_days_from_second, m3_days_from_second,m4_days_from_second, 
            m5_days_from_second, m6_days_from_second, m7_days_from_second, pre_booster_from_second, 
            boost_from_second)) 

output_df_with_dates_2 <- df_with_dates_2 %>% mutate(ID = unlist(ID))
output_df_SN_2 <- df_SN_2 %>% 
  mutate(ID = unlist(ID), value = unlist(value))

final_df_2 <- merge(x=output_df_with_dates_2, y=output_df_SN_2, by=c("ID", "time_draw"), all=TRUE)
colnames(final_df_2) <- c("ID", "time_draw", "day_draw", "SN_num")
####end####

#new labeling regime to create final_df_2
####start####
df_SN_2 <- clean_covid_data_2 %>% 
  select(c(s_00, s_01, s_02, s_03, s_04, s_05, s_06, s_07, 
           s_08, s_09, s_10, s_11, s_12, s_13, s_14, ID)) %>%
  pivot_longer(cols = s_00:s_14, names_to = "time_draw") %>%
  as.data.frame() %>%
  na_if("character(0)")

df_with_dates_2<- clean_covid_data_2 %>%
  select(c(ID, w1_since_first, w2_since_first, w3_since_first, d1_from_second, w1_from_second,
           m1_days_from_second, m2_days_from_second, m3_days_from_second,m4_days_from_second, 
           m5_days_from_second, m6_days_from_second, m7_days_from_second, pre_booster_from_second, 
           boost_from_second)) %>%
  mutate(s_01 = w1_since_first, s_02 = w2_since_first, s_03 = w3_since_first, s_04 = d1_from_second,
         s_05 = w1_from_second, s_06 = m1_days_from_second, s_07 = m2_days_from_second, 
         s_08 = m3_days_from_second, s_09 = m4_days_from_second, s_10 = m5_days_from_second, 
         s_11 = m6_days_from_second, s_12 = m7_days_from_second, s_13 = pre_booster_from_second,
         s_14 = boost_from_second) %>%
  pivot_longer(cols = s_01:s_14, names_to = "time_draw") %>%
  select(!c(w1_since_first, w2_since_first, w3_since_first, d1_from_second, w1_from_second,
            m1_days_from_second, m2_days_from_second, m3_days_from_second,m4_days_from_second, 
            m5_days_from_second, m6_days_from_second, m7_days_from_second, pre_booster_from_second, 
            boost_from_second)) 


output_df_with_dates_2 <- df_with_dates_2 %>% mutate(ID = unlist(ID))
output_df_SN_2 <- df_SN_2 %>% 
  mutate(ID = unlist(ID), value = unlist(value))

final_df_2 <- merge(x=output_df_with_dates_2, y=output_df_SN_2, by=c("ID", "time_draw"), all=TRUE)
colnames(final_df_2) <- c("ID", "time_draw", "day_draw", "SN_num")
####end####

df_with_drawnums <- final_df_2 %>%
  group_by(ID) %>%  
  filter(!is.na(SN_num)) %>%
  mutate(COUNTER = 1:n()) %>% 
  ungroup()
