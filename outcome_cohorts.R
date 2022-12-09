# ============================================================================ #
#                               Outcome cohorts                                #
#                                Núria Mercadé                                 #
#                                 01-11-2022                                   #
# ============================================================================ #


## Connect to database + load packages
source("./AURUM_CDM_connection.R")

# Cohort's ID
covidDiagnosis_id <- 1   # Broad + PC + EXCL 
AZvaccine_id      <- 2   # Astrazeneca vaccine
PFvaccine_id      <- 3   # Pfizer vaccine
covid_id          <- 4   # Covid diagnosis - washout window of 6 weeks between events
VTE_id            <- 5   # Venous Thromboembolism - washout window of 90 days between events
DVT_id            <- 6   # Deep Vein Thrombosis - washout window of 90 days between events
PE_id             <- 7   # Pulmonary Embolism - washout window of 90 days between events
HS_id             <- 8   # Haemorrhagic Stroke - washout window of 90 days between events
IS_id             <- 9   # Ischemic Stroke - washout window of 90 days between events
TIA_id            <- 10  # Transient Ischemic Attack or Transient Cerebral Ischemia - washout window of 90 days between events
MI_id             <- 11  # Myocardial Infarction - washout window of 90 days between events
MP_id             <- 12  # Myocarditis Pericarditis - washout window of 90 days between events
VACA_id           <- 13  # Ventricular Arrhythmia or Cardiac Arrest - washout window of 90 days between events
HF_id             <- 14  # Heart Failure - washout window of 90 days between events

# Vaccination period of interest
library(lubridate)
priorOUT <- 180 # wash out window for outcome event prior to covid

# Time windows between covid-19 and disease to classify pacs:
# First month after infection
date11 <- "0"
date12 <- "30"
# 31 to 90 after infection
date21 <- "31"
date22 <- "90"
# 91 to 180 after infection
date31 <- "91"
date32 <- "180"
# 181 and up to a year after infection
date41 <- "181"
date42 <- "365"


# COHORTS DB
covid_db <- cohorts_db %>% 
  filter(cohort_definition_id == covid_id) %>% 
  select(person_id = subject_id, covid_date = cohort_start_date) %>%
  compute()

# Vaccine AZ and PF
vaccine_db      <- cohorts_db %>% 
  inner_join(tibble(cohort_definition_id = c(AZvaccine_id,PFvaccine_id)), by = "cohort_definition_id", copy = TRUE) %>%
  select(person_id = subject_id, vaccine_date = cohort_start_date) %>%
  compute()

# Cohort table to fill
PASC_cohort_table <- covid_db %>% 
  transmute( cohort_definition_id = person_id, subject_id = person_id, cohort_start_date = covid_date, cohort_end_date = covid_date) %>%
  filter(subject_id == 1)


for (j in VTE_id:HF_id)
{
  # j refers to the cohort ID of the disease cohort; 
  # jj will be used to set a cohort_definition_id for the outcome cohorts
  jj <- j-4 
  
  # Disease db
  outcome_db  <- cohorts_db %>% 
    filter(cohort_definition_id == j) %>% 
    select(person_id = subject_id, outcome_date = cohort_start_date) %>% 
    compute() 
  
  # For each disease we get 16 outcome cohorts = 4 time windows after infection * 4 censoring types
  
  # ---------------------------- NO CENSORING ----------------------------- #
  # PACS db: covid + disease 
  outcome1_db <- outcome_db %>%
    inner_join(covid_db, by = "person_id") %>% 
    compute() 
  
  #  COVIDS to exclude:  have a complication event 180 days before
  covidToExclude <- outcome1_db %>% 
    group_by(person_id,covid_date) %>%
    filter(covid_date > outcome_date & covid_date < outcome_date + days(priorOUT)) %>%
    ungroup() %>% 
    compute()
  
  # Exclude covids above + ensure that disease happens within 1 year after covid
  outcomePacs11_db <- outcome1_db %>%
    anti_join(covidToExclude, by = c("person_id", "outcome_date", "covid_date")) %>% 
    group_by(person_id) %>%
    filter(covid_date < outcome_date) %>% 
    filter(covid_date >= outcome_date - days(date42)) %>% 
    ungroup() %>%  
    compute()
  
  # Add time column to classify in time intervals + delete outcome_date (disease date) column to distinct
  outcomePacs1_db <- outcomePacs11_db %>%
    mutate(distancePACS = outcome_date-covid_date) %>%
    mutate(cohort_definition_id = if_else(distancePACS >= date11 & distancePACS <= date12, 4*jj-3 + HF_id,
                                          if_else(distancePACS >= date21 & distancePACS <= date22, 4*jj-2 + HF_id,
                                                  if_else(distancePACS >= date31 & distancePACS <= date32, 4*jj-1 + HF_id,
                                                          if_else(distancePACS >= date41 & distancePACS <= date42, 4*jj + HF_id,NA))))) %>%
    select(-outcome_date,-distancePACS) %>%
    distinct() %>% 
    filter(!is.na(cohort_definition_id)) %>%
    compute()
  
  
  # Fill cohort table: cohort_id, person_id, cohort_start_date (covid), cohort_end_date (min date disease)
  PASC_cohort_table <- PASC_cohort_table %>%
    full_join(outcomePacs1_db %>%
                left_join(outcomePacs11_db %>% 
                            select(person_id, covid_date, outcome_date),by = c("person_id","covid_date")) %>%
                group_by(cohort_definition_id,person_id,covid_date) %>%
                filter(outcome_date == min(outcome_date)) %>%
                ungroup() %>% 
                rename(subject_id = person_id, cohort_start_date = covid_date, cohort_end_date = outcome_date)) %>%
    compute()
  
  
  # --------------------------- COVID CENSORING --------------------------- #
  # Associate disease outcome to the nearest prior covid
  outcomePacs21_db <- outcomePacs11_db %>%
    group_by(person_id, outcome_date) %>%
    filter(covid_date == max(covid_date)) %>%
    ungroup() %>%
    compute()
  
  # Classify by time window after covid
  outcomePacs2_db <- outcomePacs21_db %>%
    mutate(distancePACS = outcome_date-covid_date) %>%
    mutate(cohort_definition_id = if_else(distancePACS >= date11 & distancePACS <= date12, 4*jj-3 + 40 + HF_id,
                                          if_else(distancePACS >= date21 & distancePACS <= date22, 4*jj-2 + 40 + HF_id,
                                                  if_else(distancePACS >= date31 & distancePACS <= date32, 4*jj-1 + 40 + HF_id,
                                                          if_else(distancePACS >= date41 & distancePACS <= date42, 4*jj + 40 + HF_id,NA))))) %>%
    select(-outcome_date,-distancePACS) %>%
    distinct() %>% 
    filter(!is.na(cohort_definition_id)) %>%
    compute()
  
  # Fill cohort table: cohort_id, person_id, cohort_start_date (covid), cohort_end_date (min date disease)
  PASC_cohort_table <- PASC_cohort_table %>%
    full_join(outcomePacs2_db %>%
                left_join(outcomePacs21_db %>% 
                            select(person_id, covid_date, outcome_date),by = c("person_id","covid_date")) %>%
                group_by(cohort_definition_id,person_id,covid_date) %>%
                filter(outcome_date == min(outcome_date)) %>%
                ungroup() %>% 
                rename(subject_id = person_id, cohort_start_date = covid_date, cohort_end_date = outcome_date)) %>%
    compute()
  
  
  # -------------------------- VACCINE CENSORING -------------------------- #
  # excludePACSvax: covid vaccine between covid and disease
  excludePACSvax <- outcomePacs11_db %>%
    left_join(vaccine_db) %>%
    filter(covid_date < vaccine_date) %>%
    filter(vaccine_date < outcome_date) %>%
    compute()
  
  outcomePacs31_db <- outcomePacs11_db %>%
    anti_join(excludePACSvax, by = c("person_id", "outcome_date", "covid_date")) %>% compute()
  
  # Classify by time window after covid
  outcomePacs3_db <- outcomePacs31_db %>%
    mutate(distancePACS = outcome_date-covid_date) %>%
    mutate(cohort_definition_id = if_else(distancePACS >= date11 & distancePACS <= date12, 4*jj-3 + 80 + HF_id,
                                          if_else(distancePACS >= date21 & distancePACS <= date22, 4*jj-2 + 80 + HF_id,
                                                  if_else(distancePACS >= date31 & distancePACS <= date32, 4*jj-1 + 80 + HF_id,
                                                          if_else(distancePACS >= date41 & distancePACS <= date42, 4*jj + 80 + HF_id,NA))))) %>%
    select(-outcome_date,-distancePACS) %>%
    distinct() %>% 
    filter(!is.na(cohort_definition_id)) %>%
    compute()
  
  # Fill cohort table: cohort_id, person_id, cohort_start_date (covid), cohort_end_date (min date disease)
  PASC_cohort_table <- PASC_cohort_table %>%
    full_join(outcomePacs3_db %>%
                left_join(outcomePacs31_db %>% 
                            select(person_id, covid_date, outcome_date),by = c("person_id","covid_date")) %>%
                group_by(cohort_definition_id,person_id,covid_date) %>%
                filter(outcome_date == min(outcome_date)) %>%
                ungroup() %>% 
                rename(subject_id = person_id, cohort_start_date = covid_date, cohort_end_date = outcome_date)) %>%
    compute()
  
  
  # ---------------------- COVID + VACCINE CENSORING ---------------------- #
  # The last two cohorts together
  outcomePacs41_db <- outcomePacs21_db %>% inner_join(outcomePacs31_db) %>% compute()
  
  # Classify by time window after covid
  outcomePacs4_db <- outcomePacs41_db %>%
    mutate(distancePACS = outcome_date-covid_date) %>%
    mutate(cohort_definition_id = if_else(distancePACS >= date11 & distancePACS <= date12, 4*jj-3 + 120 + HF_id,
                                          if_else(distancePACS >= date21 & distancePACS <= date22, 4*jj-2 + 120 + HF_id,
                                                  if_else(distancePACS >= date31 & distancePACS <= date32, 4*jj-1 + 120 + HF_id,
                                                          if_else(distancePACS >= date41 & distancePACS <= date42, 4*jj + 120 + HF_id,NA))))) %>%
    select(-outcome_date,-distancePACS) %>%
    distinct() %>% 
    filter(!is.na(cohort_definition_id)) %>%
    compute()
  
  # Fill cohort table: cohort_id, person_id, cohort_start_date (covid), cohort_end_date (min date disease)
  PASC_cohort_table <- PASC_cohort_table %>%
    full_join(outcomePacs4_db %>%
                left_join(outcomePacs41_db %>% 
                            select(person_id, covid_date, outcome_date),by = c("person_id","covid_date")) %>%
                group_by(cohort_definition_id,person_id,covid_date) %>%
                filter(outcome_date == min(outcome_date)) %>%
                ungroup() %>% 
                rename(subject_id = person_id, cohort_start_date = covid_date, cohort_end_date = outcome_date)) %>%
    compute()
}




# ------------------------ Add cohorts to the database ----------------------- #

# sql_query <- glue::glue("INSERT INTO {results_database_schema}.{cohort_tab_name}\n",
#                         " SELECT * FROM (\n",
#                         dbplyr::sql_render(cohort_table),
#                         "\n) AS from_table")
# 
# DBI::dbExecute(con, as.character(sql_query))