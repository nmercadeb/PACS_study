# ============================================================================ #
#                                Data Analysis                                 #
#                                Núria Mercadé                                 #
#                                 28-11-2022                                   #
# ============================================================================ #

## Connect to database
source("./AURUM_CDM_connection.R")

## Packages 
library("survival")
library("cmprsk")
library("ggplot2")
library("ggpubr")
library("patchwork")

load(file = "./RData/s01_individuals_weight.Rdata")
load(file = "./RData/s01_individuals_AZ.Rdata")
load(file = "./RData/s01_individuals_PF.Rdata")

individuals_weight <- individuals_weight[individuals_weight$weight != 0, ] 

# Outcome cohorts depending on censoring
nocens    <- 1:40 + 14
covcens   <- 41:80 + 14
vaxcens   <- 81:120 + 14
covaxcens <- 121:160 + 14

# Type of censoring for the unvaccinated cohort: 
# type 1 meaning complete censoring, type 2 meaning not cens for vax between covid and disease

index      <- as_tibble(read.csv(file = "./Data/estudi_pacs_03.csv")[,2:5])
groupNames <- c("vaccinated", "unvaccinated")

# Functions ---
studyTable_unvaxCorrection <- function(pop, cohorts_id, groups, cohorts_table, type2_data) {
  # In this funtion the unvax cohort is joined with an output cohort with vax censoring
  outcome_unvax <- cohorts_table %>% inner_join(tibble(cohort_definition_id = cohorts_id), by = "cohort_definition_id", copy = TRUE) %>%
    select(person_id = subject_id, covid_date = cohort_start_date, disease_date = cohort_end_date) %>% collect() 
  pop_unvax <- pop[pop$group == "unvaccinated",]
  
  data_unvax <- pop_unvax %>% left_join(outcome_unvax)
  minimTime_unvax <- pmin(data_unvax$death_date, data_unvax$leave_date,
                          data_unvax$covid_date, data_unvax$next_vaccine, na.rm = TRUE)
  data_unvax <- data_unvax %>% mutate(time = as.integer(minimTime_unvax - index_date))
  data_unvax$event <- 0
  data_unvax$event[minimTime_unvax == data_unvax$covid_date] <- 1 # Event of interest
  data_unvax$event[minimTime_unvax == data_unvax$death_date] <- 2 # Competing event
  
  data_unvax$event <- factor(data_unvax$event, levels = c(0,1,2), labels = c("Censored","PACS","Death"))
  data_unvax$group <- as.factor(data_unvax$group) %>% relevel(group, ref = groups[2])
  
  stopifnot(data_unvax$event[data_unvax$disease_date > data_unvax$next_vaccine & !is.na(data_unvax$disease_date) & !is.na(data_unvax$next_vaccine)] == "Censored")
  
  # data_unvax <- data_unvax %>% select(person_id, group, time, event, weight)
  
  outcome_table <- data_unvax %>% union_all(type2_data[type2_data$group == "vaccinated",]) %>% distinct()
  
  return(outcome_table)
}
studyTable_same <- function(pop, cohorts_id, groups, cohorts_table, censortype) {
  # In this function vax and unvax are joined with the same output cohort
  outcome_table <- cohorts_table %>% inner_join(tibble(cohort_definition_id = cohorts_id), by = "cohort_definition_id", copy = TRUE) %>%
    select(person_id = subject_id, covid_date = cohort_start_date, disease_date = cohort_end_date) %>% collect() 
  
  outcome_table <- pop %>% left_join(outcome_table)
  
  if ({{censortype}} == "none" | {{censortype}} == "covid") {
    data_vax   <- outcome_table[outcome_table$group == "vaccinated",]
    data_unvax <- outcome_table[outcome_table$group == "unvaccinated",]
    
    minimTime_vax   <- pmin(data_vax$death_date, data_vax$leave_date,
                            data_vax$covid_date, na.rm = TRUE)
    minimTime_unvax <- pmin(data_unvax$death_date, data_unvax$leave_date,
                            data_unvax$covid_date, data_unvax$next_vaccine, na.rm = TRUE)
    
    data_vax <- data_vax %>% mutate(time = as.integer(minimTime_vax - index_date))
    data_vax$event <- 0
    data_vax$event[minimTime_vax == data_vax$covid_date] <- 1 # Event of interest
    data_vax$event[minimTime_vax == data_vax$death_date] <- 2 # Competing event
    
    data_unvax <- data_unvax %>% mutate(time = as.integer(minimTime_unvax - index_date))
    data_unvax$event <- 0
    data_unvax$event[minimTime_unvax == data_unvax$covid_date] <- 1 # Event of interest
    data_unvax$event[minimTime_unvax == data_unvax$death_date] <- 2 # Competing event
    
    outcome_table <- data_vax %>% union_all(data_unvax)
    
  } else if ({{censortype}} == "vax" | {{censortype}} == "both") {
    minimTime   <- pmin(outcome_table$death_date, outcome_table$leave_date,
                        outcome_table$covid_date, outcome_table$next_vaccine, na.rm = TRUE)
    outcome_table <- outcome_table %>% mutate(time = as.integer(minimTime - index_date))
    outcome_table$event <- 0
    outcome_table$event[minimTime == outcome_table$covid_date] <- 1 # Event of interest
    outcome_table$event[minimTime == outcome_table$death_date] <- 2 # Competing event
  } 
  
  if ({{censortype}} == "none") {} else if ({{censortype}} == "covid") {
  } else if ({{censortype}} == "vax") {
    # CENSOR VACCINE BEFORE DISEASE FOR BOTH GROUPS
    stopifnot(outcome_table$event[outcome_table$disease_date > outcome_table$next_vaccine & !is.na(outcome_table$disease_date) & !is.na(outcome_table$next_vaccine)] == 0)
  } else if ({{censortype}} == "both") {
    # CENSOR VACCINE BEFORE DISEASE FOR BOTH GROUPS + COVID BEFORE DISEASE BOTH GROUPS
    stopifnot(outcome_table$event[outcome_table$disease_date > outcome_table$next_vaccine & !is.na(outcome_table$disease_date) & !is.na(outcome_table$next_vaccine)] == 0)
  }
  
  # outcome_table <- outcome_table %>% select(person_id, group, time, event, weight)
  
  outcome_table$event <- factor(outcome_table$event, levels = c(0,1,2), labels = c("Censored","PACS","Death"))
  outcome_table$group <- as.factor(outcome_table$group) %>% relevel(group, ref = groups[2])
  outcome_table <- outcome_table %>% distinct()
  
  return(outcome_table)
}
fg           <- function(table) {
  fg_data <- finegray(Surv(time, event) ~ ., data=table, weights = weight)
  fg_regression <- coxph(Surv(fgstart, fgstop, fgstatus) ~ group, weight=fgwt, data=fg_data)
  coef <-  c(summary(fg_regression)$conf.int, summary(fg_regression)$coefficients[3])
  return(coef)
}   

# Nomenclatura ---
## Censoring in unvaccinated covid:
#  - full: complete censoring for vaccine
#  - partial: partial censoring for vaccine
## Censoring overall:
#  - 1: default (no cens)
#  - 2: covid censoring
#  - 3: vaccine censoring
#  - 4: covid and vaccine censoring

for (poblacio in c("", "AZ", "PF")) { # '' refers to any vaccine, AZ refers to only AZ vaccine and same for PF
  if(poblacio == "") {
    individuals_cru <- individuals_weight %>% mutate(weight = 1)
  } else if (poblacio == "AZ") {
    individuals_cru <- individuals_AZ %>% mutate(weight = 1)
    individuals_weight <- individuals_AZ
    individuals_weight <- individuals_weight[individuals_weight$weight != 0, ] 
    individuals_cru <- individuals_weight %>% mutate(weight = 1)
  } else if (poblacio == "PF") {
    individuals_weight <- individuals_PF
    individuals_weight <- individuals_weight[individuals_weight$weight != 0, ] 
    individuals_cru <- individuals_weight %>% mutate(weight = 1)
  }
  
  for (analisis in c("cru", "adjusted")) { # 'cru' és l'analisi amb weights = 1; 
                                           # quan es fa 'adjusted' s'extreuen dos sHR, el dels weights i el dels weights + calibration
    
    if (analisis == "cru") {population <- individuals_cru} else {population <- individuals_weight}
    
    for (kk in 1:4) {  # Tenim 4 different outcome cohorts depenent del censoring, a més, 
                       # tenim en 2 d'ells (default i covid) dos sensitivity analysis per la cohort de no vacunats
                       # kk fa referencia als 4 outcome cohorts, mentre que les dataframes per omplir amb els resultats seran 2 si
                       # estem en deafult i covid (kk = 1 i 2) o una si estem en vaccine i vaccine + covid (kk = 3 i 4)
      
      # Create tables (type 1 and 2) to fill for each overall censoring
      if (kk == 1) { cens = "none"; cohortIDs = nocens} else if (kk == 2) { cens = "covid"; cohortIDs = covcens } else if (kk == 3) { cens = "vax"; cohortIDs = vaxcens } else if (kk == 4) { cens = "both"; cohortIDs = covaxcens }
      
      outcomeHR_full <- index[cohortIDs,c(1,2)]
      outcomeHR_full <- outcomeHR_full %>% rbind(tibble(cohortId = NA, cohortName = gsub("VTE","ATE",outcomeHR_full$cohortName[1:4])))
      outcomeHR_full <- outcomeHR_full %>% cbind(tibble(HR = NA, HRinvers = NA, low95 = NA, upper95 = NA, SE = NA, vacPACS = NA, unvacPACS = NA))
      
      if (kk %in% 1:2) {outcomeHR_partial <- outcomeHR_full}
      
      for (pos in 1:nrow(outcomeHR_full)) {
        print(paste0('Status: poblacio = ', poblacio, ', analisis = ', analisis, ', k = ', kk, ', pos = ', pos))
        # Get cohort outcome ID for type 2 
        if (pos <= 40 & pos >= 5) { id <- cohortIDs[pos]
        } else if (pos <= 4) {id <- c(outcomeHR_full$cohortId[gsub("VTE","DVT",outcomeHR_full$cohortName[pos])  == outcomeHR_full$cohortName],
                                      outcomeHR_full$cohortId[gsub("VTE","PE",outcomeHR_full$cohortName[pos])   == outcomeHR_full$cohortName])
        } else if (pos >= 41) {id <- c(outcomeHR_full$cohortId[gsub("ATE","IS",outcomeHR_full$cohortName[pos])  == outcomeHR_full$cohortName],
                                       outcomeHR_full$cohortId[gsub("ATE","TIA",outcomeHR_full$cohortName[pos]) == outcomeHR_full$cohortName],
                                       outcomeHR_full$cohortId[gsub("ATE","MI",outcomeHR_full$cohortName[pos])  == outcomeHR_full$cohortName])}
        
        # Divide if we need unvax sensitibity analisis or not (1&2 we need)
        if (kk %in% 1:2) {
          studyTable_partial <- studyTable_same(population, id, groupNames, cohorts_db, cens)
          
          if (sum(studyTable_partial$event == "PACS") >= 5) {
            outcomeHR_partial$vacPACS[pos] <- sum(studyTable_partial$event == "PACS" & studyTable_partial$group == "vaccinated")
            outcomeHR_partial$unvacPACS[pos] <- sum(studyTable_partial$event == "PACS" & studyTable_partial$group == "unvaccinated")
            fineGrey_2       <- fg(studyTable_partial)
            outcomeHR_partial[pos, 3:7]  <- fineGrey_2 
            
            # Full vax censor for unvax:
            if (kk == 1) {names <- gsub("no censoring","vax censoring",outcomeHR_full$cohortName[outcomeHR_full$cohortId %in% id])
            } else {names <- gsub("covid censoring","covid+vax censoring",outcomeHR_full$cohortName[outcomeHR_full$cohortId %in% id])}
            id_censor_cohorts <- index$cohortId[index$cohortName %in% names]
            studyTable_full <- studyTable_unvaxCorrection(population, id_censor_cohorts, groupNames, cohorts_db, studyTable_partial)
            
            # See if there's enough cases:
            if (sum(studyTable_full$event == "PACS") >= 5) {
              outcomeHR_full$vacPACS[pos] <- sum(studyTable_full$event == "PACS" & studyTable_full$group == "vaccinated")
              outcomeHR_full$unvacPACS[pos] <- sum(studyTable_full$event == "PACS" & studyTable_full$group == "unvaccinated")
              fineGrey_2       <- fg(studyTable_full)
              outcomeHR_full[pos, 3:7]  <- fineGrey_2}}
          
        } else {
          studyTable <- studyTable_same(population, id, groupNames, cohorts_db, cens)
          if (sum(studyTable$event == "PACS") >= 5) {
            outcomeHR_full$vacPACS[pos] <- sum(studyTable$event == "PACS" & studyTable$group == "vaccinated")
            outcomeHR_full$unvacPACS[pos] <- sum(studyTable$event == "PACS" & studyTable$group == "unvaccinated")
            fineGrey       <- fg(studyTable)
            outcomeHR_full[pos, 3:7]  <- fineGrey}
        }
      }
      
      
      # Check if we are exploring unvax sensitivity
      fileName    <- paste0('./s01_RiskEstimates/', analisis, poblacio, '/', analisis, '_', cens)
      fileNameCal <- paste0('./s01_RiskEstimates/calibrated', poblacio, '/calibrated_', cens)
      
      if (kk %in% 1:2) {
        write.csv(outcomeHR_full, file = paste0(fileName, '_full.csv'))
        write.csv(outcomeHR_partial, file = paste0(fileName, '_partial.csv'))
        
        # Check if we need calibration
        if (analisis == "adjusted") {
          load("./RData/FE_NCO_after.RData")
          model <- fitSystematicErrorModel(log(toPlot_NCO_after$RelativeRisk), toPlot_NCO_after$SD, rep(0,nrow(toPlot_NCO_after)))
          
          result_full <- calibrateConfidenceInterval(log(outcomeHR_full$HR), outcomeHR_full$SE, model,ciWidth = 0.95)
          outcomeHR_full_cal <- outcomeHR_full[,1:2] %>% cbind(exp(result_full))
          colnames(outcomeHR_full_cal)[3:6] <- c("HR","low95", "upper95","SE")
          write.csv(outcomeHR_full_cal, file = paste0(fileNameCal, '_full.csv'))
          
          result_partial <- calibrateConfidenceInterval(log(outcomeHR_partial$HR), outcomeHR_partial$SE, model,ciWidth = 0.95)
          outcomeHR_partial_cal <- outcomeHR_partial[,1:2] %>% cbind(exp(result_partial))
          colnames(outcomeHR_partial_cal)[3:6] <- c("HR","low95", "upper95","SE")
          write.csv(outcomeHR_partial_cal, file = paste0(fileNameCal, '_partial.csv')) }
      } else {
        write.csv(outcomeHR_full, file = paste0(fileName, '.csv'))
        
        # Check if we need calibration
        if (analisis == "adjusted") {
          result_full <- calibrateConfidenceInterval(log(outcomeHR_full$HR), outcomeHR_full$SE, model,ciWidth = 0.95)
          outcomeHR_full_cal <- outcomeHR_full[,1:2] %>% cbind(exp(result_full))
          colnames(outcomeHR_full_cal)[3:6] <- c("HR","low95", "upper95","SE")
          write.csv(outcomeHR_full_cal, file = paste0(fileNameCal, '.csv'))}
      }
    }
  }
}