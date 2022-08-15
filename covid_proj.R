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

raw_covid_data <- read.csv("vaccine sample info.csv")

# create clean data set of all info with days post vaccines
####start####
clean_covid_data <- raw_covid_data %>% 
  janitor::clean_names() %>% 
  rename(x1_months_post_2nd_dose = x1_month_post_2nd_dose,
         x4_months_post_2nd_dose =  x4_months_post_2nd_dose2) %>% 
  mutate(predraw_date = extract_n_clean_dates(predraw),
         x_1_wk_post_first_dose_date = extract_n_clean_dates(x_1_wk_post_first_dose),
         x1_2_days_post_2nd_dose_date = extract_n_clean_dates(x1_2_days_post_2nd_dose),
         post_second_dose_date = extract_n_clean_dates(post_2nd_dose),
         
         x1_months_post_2nd_dose_date = extract_n_clean_dates(x1_months_post_2nd_dose),
         x2_months_post_2nd_dose_date = extract_n_clean_dates(x2_months_post_2nd_dose),
         x3_months_post_2nd_dose_date = extract_n_clean_dates(x3_months_post_2nd_dose),
         x4_months_post_2nd_dose_date = extract_n_clean_dates(x4_months_post_2nd_dose),
         x5_months_post_2nd_dose_date = extract_n_clean_dates(x5_months_post_2nd_dose),
         x6_months_post_2nd_dose_date = extract_n_clean_dates(x6_months_post_2nd_dose),
         x7_months_post_2nd_dose_date = extract_n_clean_dates(x7_months_post_2nd_dose)
         )  %>% 
  mutate_at(vars(contains("date")), 
            ~as.Date(., format = "%m/%d/%Y")) %>% 
  mutate(days_since_first = difftime(x_1_wk_post_first_dose_date, vaccine_1st_dose_date, units = "days"),
         m1_days_from_second = difftime(x1_months_post_2nd_dose_date, vaccine_2nd_dose_date, units = "days"),
         m2_days_from_second = difftime(x2_months_post_2nd_dose_date, vaccine_2nd_dose_date, units = "days"),
         m3_days_from_second = difftime(x3_months_post_2nd_dose_date, vaccine_2nd_dose_date, units = "days"),
         m4_days_from_second = difftime(x4_months_post_2nd_dose_date, vaccine_2nd_dose_date, units = "days"),
         m5_days_from_second = difftime(x5_months_post_2nd_dose_date, vaccine_2nd_dose_date, units = "days"),
         m6_days_from_second = difftime(x6_months_post_2nd_dose_date, vaccine_2nd_dose_date, units = "days"),
         m7_days_from_second = difftime(x7_months_post_2nd_dose_date, vaccine_2nd_dose_date, units = "days")
         ) %>%
  mutate(s_p = stringr::str_extract_all(predraw, "^SN\\d{3}"),
         s_f = stringr::str_extract_all(x_1_wk_post_first_dose, "^SN\\d{3}"),
         s_1mon = stringr::str_extract_all(x1_months_post_2nd_dose, "^SN\\d{3}"),
         s_2mon = stringr::str_extract_all(x2_months_post_2nd_dose, "^SN\\d{3}"),
         s_3mon = stringr::str_extract_all(x3_months_post_2nd_dose, "^SN\\d{3}"),
         s_4mon = stringr::str_extract_all(x4_months_post_2nd_dose, "^SN\\d{3}"),
         s_5mon = stringr::str_extract_all(x5_months_post_2nd_dose, "^SN\\d{3}"),
         s_6mon = stringr::str_extract_all(x6_months_post_2nd_dose, "^SN\\d{3}"),
         s_7mon = stringr::str_extract_all(x7_months_post_2nd_dose, "^SN\\d{3}")
         ) %>%
  mutate(ID = stringr::str_extract_all(subject_id_for_vaccine_data_purposes, "\\d{2}")
         )
####end####

df_SN <- clean_covid_data %>% 
  select(c(s_p,s_f, s_1mon, s_2mon, s_3mon, s_4mon, s_5mon, s_6mon, s_7mon, ID)) %>%
  pivot_longer(cols = s_p:s_7mon, names_to = "time_draw") %>%
  as.data.frame() %>%
  na_if("character(0)")

df_with_dates <- clean_covid_data %>%
  select(c(ID, days_since_first, m1_days_from_second, m2_days_from_second, m3_days_from_second, 
           m4_days_from_second, m5_days_from_second, m6_days_from_second, m7_days_from_second)) %>%
  mutate(s_f = days_since_first, s_1mon = m1_days_from_second, s_2mon = m2_days_from_second, 
         s_3mon = m3_days_from_second, s_4mon = m4_days_from_second, s_5mon = m5_days_from_second, 
         s_6mon = m6_days_from_second, s_7mon = m7_days_from_second) %>%
  pivot_longer(cols = s_f:s_7mon, names_to = "time_draw") %>%
  select(!c(days_since_first, m1_days_from_second, m2_days_from_second, m3_days_from_second, m4_days_from_second,
           m5_days_from_second, m6_days_from_second, m7_days_from_second)) 

output_df_with_dates <- df_with_dates %>% mutate(ID = unlist(ID))
output_df_SN <- df_SN %>% 
  mutate(ID = unlist(ID), value = unlist(value))

final_df <- merge(x=output_df_with_dates, y=output_df_SN, by=c("ID", "time_draw"), all=TRUE)
colnames(final_df) <- c("ID", "time_draw", "day_draw", "SN_num")

write.csv(final_df, "ID_to_SN_day.csv")


