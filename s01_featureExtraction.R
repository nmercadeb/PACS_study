# ============================================================================ #
#                         Feature Extraction - SMD - PS                        #
#                                Núria Mercadé                                 #
#                                 28-10-2022                                   #
# ============================================================================ #

## Connect to database + packages
source("./AURUM_CDM_connection.R")

## Packages, variables, tables and functions needed
source("./s01_forFeatureExtraction.R")

# =========================== TABLES OF INTEREST ============================= # ----
# For feature extraction we are interested in the following tables from cdm:
# - drug_era
# - condition_occurrence
# - visit_occurrence
# Also we want
# - Covid-19 testing (we get it from 2 cohorts)
# - Age and gender, we get it from person table from cdm
# - Location (region + care_site_id) we get it from care_site, person and location tables from cdm

print("STEP 1: get tables from database")
## DRUG ERA ---------------------------
drug_patients <- cdm$drug_era %>% 
  select(person_id,drug_concept_id,drug_era_start_date,drug_era_end_date) %>%
  inner_join(list_id) %>%
  compute()

# Concept ID as covariate ID + concept name as covariate name
drug_era_id <- drug_patients %>%
  select(drug_concept_id) %>%
  rename("concept_id" = "drug_concept_id") %>%
  distinct() %>%
  left_join(concept_db %>% select(concept_id,concept_name)) %>%
  collect %>%
  rename("covariateId" = "concept_id") %>%
  rename("covariateName" = "concept_name") %>%
  mutate(FeatureExtractionId = as.numeric(paste0(covariateId,"001"))) %>%
  mutate(AnalysisId = 1)

# Add to patients table the featureExtractionID (covariate ID + window1: 001), Era_start_date, Era_end_date
drug_patients <- drug_patients %>%
  collect() %>%
  rename("covariateId" = "drug_concept_id") %>%
  left_join(drug_era_id) %>%
  select(person_id,FeatureExtractionId,drug_era_start_date,drug_era_end_date) %>%
  rename("Era_start_date" = "drug_era_start_date") %>%
  rename("Era_end_date" = "drug_era_end_date") %>%
  mutate(Era_start_date = as.Date(Era_start_date)) %>%
  mutate(Era_end_date = as.Date(Era_end_date))


print("Drug era done")

## CONDITION OCCURRENCE ---------------
CO_patients <- cdm$condition_occurrence %>%
  select(person_id,condition_concept_id,condition_start_date) %>%
  inner_join(list_id) %>%
  compute()

CO_id <- CO_patients %>%
  select(condition_concept_id) %>%
  rename("concept_id"="condition_concept_id") %>%
  distinct() %>%
  left_join(concept_db %>% select(concept_id,concept_name)) %>%
  collect() %>%
  rename("covariateId" = "concept_id") %>%
  rename("covariateName" = "concept_name") %>%
  mutate(FeatureExtractionId = as.numeric(paste0(covariateId,"002"))) %>%
  mutate(AnalysisId = 2)

CO_patients <- CO_patients %>%
  collect() %>%
  rename("covariateId"="condition_concept_id") %>%
  left_join(CO_id) %>%
  select(person_id,FeatureExtractionId,condition_start_date) %>%
  rename("Event_date"="condition_start_date") %>%
  mutate(Event_date = as.Date(Event_date))

save(list = c("CO_patients", "CO_id"), file = "./RData/s01_CO_forNCO.RData")

print("Condition occurrence done")

## VISITS -----------------------------
VI_patients <- cdm$visit_occurrence %>%
  select(person_id,visit_start_date) %>%
  inner_join(list_id) %>%
  distinct() %>%
  collect() %>%
  rename("Event_date"="visit_start_date") %>%
  mutate(FeatureExtractionId = 581477003)

VI_id <- tibble(FeatureExtractionId = 581477003,covariateId = 581477, covariateName = "Visit to GP", AnalysisId = 3)

print("Visits done")


## COVID19 diagnosis and testing -------
test_patients1 <- list_id %>%
  inner_join(cohorts_db %>% filter(cohort_definition_id %in% c(test_covid_any_id ,test_covid_pcr_id)) %>% rename("person_id"="subject_id")) %>%
  collect() %>%
  mutate(FeatureExtractionId = ifelse(cohort_definition_id==test_covid_any_id,1006,2006)) %>%
  select(-cohort_definition_id) %>%
  rename("Event_date" = "cohort_start_date")
test_patients2 <- list_id %>%
  inner_join(cohorts_db %>% filter(cohort_definition_id == covid_diagonsis_id) %>% rename("person_id"="subject_id")) %>%
  collect() %>%
  mutate(FeatureExtractionId = 3006) %>%
  select(-cohort_definition_id) %>%
  rename("Event_date" = "cohort_start_date")
test_patients <- rbind(test_patients1,test_patients2)
test_id <- tibble(covariateId = c(NA,NA,NA), covariateName = c("Any covid19 test","PCR covid19 test","Diag+ covid19"), FeatureExtractionId = c(1006,2006,3006), AnalysisId = c(6,6,6))


print("Covid-19 diagnosis and testing done")


## AGE --------------------------------
age_group <- list_id %>%
  left_join(person_db) %>%
  select(person_id,year_of_birth) %>%
  collect() %>%
  mutate(month_of_birth = 1) %>%
  mutate(day_of_birth   = 1) %>%
  mutate(dob = as.Date(dmy(paste(day_of_birth,month_of_birth,year_of_birth,sep="-")))) %>%
  mutate(age = floor(as.numeric(difftime(date1,dob,unit="days"))/365.25)) %>%
  mutate(age_group = cut(age, c(75,80,85,90,95,100,105,110,115,120),labels = c("75 to 79 years","80 to 84 years","85 to 89 years","90 to 94 years","95 to 99 years","100 to 104 years","105 to 109 years","110 to 114 years","115 to 119 years"),include.lowest = TRUE, right = FALSE)) %>%
  mutate(agegid = as.numeric(age_group)) %>%
  mutate(FeatureExtractionId = 1008)
age_group_patients <- age_group %>% select(person_id,FeatureExtractionId,agegid) %>% rename("value"="agegid")
age_group_id       <- age_group %>% select(FeatureExtractionId,agegid,age_group) %>% distinct() %>% arrange(agegid) %>% mutate(AnalysisId = 8) %>% rename("value"="agegid") %>% rename("valueName"="age_group")

print("Age done")


## GENDER -----------------------------
gender_patients <- list_id %>%
  left_join(person_db) %>%
  select(person_id,gender_concept_id) %>%
  collect() %>%
  mutate(value = if_else(gender_concept_id==8532,1,2)) %>%
  mutate(FeatureExtractionId = 1009) %>%
  select(-gender_concept_id)
gender_id <- tibble(FeatureExtractionId = c(1009,1009), value = c(1,2), valueName = c("Female","Male"), AnalysisId = c(9,9))
AnalysisRef  <- rbind(AnalysisRef,c(9,"Gender","Factor"))

print("Gender done")

## LOCATION ID ------------------------ 
location_ids <- list_id %>%
  left_join(person_db    %>% select(person_id,care_site_id))   %>%
  left_join(care_site_db %>% select(care_site_id,location_id)) %>%
  mutate(location_id = if_else(location_id==10,8,location_id)) %>%
  left_join(location_db  %>% select(location_id,location_source_value)) %>%
  collect() %>%
  mutate(location_source_value = if_else(location_id==8,"South East-Central",location_source_value))

csi                <- location_ids %>% select(care_site_id) %>% distinct() %>% arrange(care_site_id) %>% mutate(value = 1:length(care_site_id)) %>% mutate( valueName = paste("care site id =",care_site_id)) %>% mutate(FeatureExtractionId = 1010)
location_patients1 <- location_ids %>% left_join(csi) %>% select(person_id,value,FeatureExtractionId) 
loc                <- location_ids %>% select(location_id,location_source_value) %>% distinct() %>% arrange(location_id) %>% mutate(value = location_id) %>% rename("valueName" = "location_source_value") %>% mutate(FeatureExtractionId = 2010)
location_patients2 <- location_ids %>% left_join(loc) %>% select(person_id,value,FeatureExtractionId) 
location_patients  <- rbind(location_patients1,location_patients2)
location_id        <- rbind(loc %>% select(-location_id) %>% mutate(AnalysisId = 10), csi %>% select(-care_site_id) %>% mutate(AnalysisId = 10))

print("Location done")

print("STEP 2: get individual tables")
## INDIVIDUAL TABLES -------------------
# Here the function 'getIndividualTabs' from the script forFeatureExtraction is used to get tables with the
# the format: person_id, covariate, value

drug_table      <- getIndividualTabs(drug_era_id, drug_patients, individuals_id, "drug_era", 2, TRUE)
co_table        <- getIndividualTabs(CO_id, CO_patients, individuals_id, "condition_ocurrence", 3, TRUE)
vi_table        <- getIndividualTabs(VI_id, VI_patients, individuals_id, "visit_ocurrence", 3, FALSE)
test_table      <- getIndividualTabs(test_id, test_patients, individuals_id, "covid19_and_testing", 3, FALSE)
age_group_table <- getIndividualTabs(age_group_id, age_group_patients, individuals_id, "age_group", 0, TRUE)
gend_table      <- getIndividualTabs(gender_id, gender_patients, individuals_id, "gender", 0, TRUE)
loc_table       <- getIndividualTabs(location_id, location_patients, individuals_id, "location", 0, TRUE)
age_table       <- age_group %>% mutate(covariate = "age") %>% select(person_id, covariate, value = age) %>% distinct()
# Table with the same format containing the years in observation:
obs_years_table <- observation_period_db %>% inner_join(list_id) %>% select(person_id,observation_period_start_date) %>% collect() %>%
  inner_join(individuals_id) %>% mutate(value = as.integer(index_date-observation_period_start_date)/365.25) %>% 
  filter(value > 0) %>% mutate(covariate = "observation_years") %>% select(person_id,covariate,value)

save(list = c("drug_table", "co_table", "vi_table", "test_table",
              "age_group_table", "gend_table", "loc_table", "age_table", 
              "obs_years_table"),
     file = "RData/s01_FE_individualTabs_2.Rdata")

print("STEP 3: compute discrete, continuos, and factor tables, + eliminate non frequent covariates")
discrete_table   <- drug_table %>% union_all(co_table) 
continuous_table <- vi_table %>% union_all(test_table) %>% union_all(age_table) %>% union_all(obs_years_table) %>% ungroup()
factor_table     <- age_group_table %>% union_all(gend_table) %>% union_all(loc_table)

# Only covariates from the discrete table with more than 0.5% frequency
minimum.freq.covariates <- 0.005
min_cov <- floor(nrow(individuals_id)*minimum.freq.covariates)
exclude <- discrete_table %>% group_by(covariate) %>% tally() %>% filter(n<=min_cov) %>% select(covariate)

print(paste0("DISCRETE: Exclude ", exclude %>% tally() %>% as.numeric(), " covariates out of ", length(unique(discrete_table$covariate))))

discrete_table         <- discrete_table %>% anti_join(exclude) 

save(list = c("discrete_table", "continuous_table", "factor_table"), 
     file = "RData/s01_FE_prep.Rdata")


# ============================= PROPENSITY SCORE ============================= # ----
## LASSO REGRESSION ---------------------
# Data for LASSO regression --> ALL except those compulsory
# discrete table: ALL
# continuous table: none (all are compulsory)
# factor table: none (all are complusory)

print("STEP 4: lasso regression")
lasso_data <- individuals_id %>% select(person_id) %>% inner(discrete_table) %>%
  pivot_wider(names_from = covariate, values_from = value, values_fill = 0) %>%
  left_join(setGroups(cohort_id_groups, groupNames)) %>% mutate(group = as.factor(group)) %>%
  mutate(group = relevel(group, ref = groupNames[2]))

X <- as.matrix(lasso_data %>% select(-person_id,-group))
Y <- lasso_data$group

print("Data for lasso obtained, starting regression")

lambdas   <- 10^seq(2, -3, by = -.1)
lasso_reg <- cv.glmnet(X, Y, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5, family="binomial")
coef.lasso_reg <- coef(lasso_reg, s = lasso_reg$lambda.1se)
var.lasso_reg  <- names(coef.lasso_reg[(coef.lasso_reg[,1]!=0),1])
var.lasso_reg  <- gsub("`", "", var.lasso_reg)

save(list = c("lasso_data","var.lasso_reg"), file = "RData/s01_FE_lasso_data.Rdata")
save(list = c("lasso_reg"), file = "RData/s01_FE_lasso_regression.Rdata")

print("STEP 5: get PROPENSITY SCORES")
## PROPENSITY SCORES ---------------------
nam             <- colnames(lasso_data)
id_col          <- c(1,length(nam),which(nam %in% var.lasso_reg)) # Take person_id, group and lasso resulting features

# Get one of the columns from linear depended covariates so ensure that predictor variables are linear independent:
linearIndependent <- tibble(covariate = c(unique(factor_table$covariate[grep("care_site",factor_table$covariate)])[1],
                                          unique(factor_table$covariate[grep("region",factor_table$covariate)])[1], 
                                          unique(factor_table$covariate[grep("age_g",factor_table$covariate)])[1],
                                          unique(factor_table$covariate[grep("gender",factor_table$covariate)])[1]))

# Get the table with the compulsory variables
compulsory_features <- individuals_id %>% select(person_id) %>% 
  left_join(continuous_table %>% union_all(factor_table %>% anti_join(linearIndependent))) %>% 
  pivot_wider(names_from = covariate, values_from = value, values_fill = 0) 

# Join data from lasso with the compulsory features + add column age^2
data_base <- lasso_data[,id_col] %>% 
  left_join(compulsory_features) %>% left_join(individuals_id %>% mutate(index_date = as.numeric(index_date)))
data_base$age2  <- data_base$age^2

# PS regression model created with a subset:
data_base_subset <- data_base[sample(nrow(data_base), 200000, replace = FALSE),]

print("Data for PS obtained, starting regression + prediction")
model_base    <- glm(formula = group ~ . - person_id, family = 'binomial', data = data_base_subset)
PS     <- predict(model_base,data_base,type = "response")

propensityScores <- data_base %>% select(person_id,group) %>% mutate(ps = PS)

save(list=c("model_base"), file = "RData/s01_FE_PS_model_base.Rdata")
save(list=c("compulsory_features","data_base","PS", "propensityScores"), file = "RData/s01_FE_PS_data.Rdata")
write.csv(propensityScores, file = "./Data/s01_ps.csv")

# Save individuals with their weights:
individuals_weight <- individuals_cru %>% inner_join(propensityScores) %>%
  mutate(weight = ifelse(group == "vaccinated", 1-ps, ps)) %>% select(-ps)
save("individuals_weight", file = "./RData/s01_individuals_weight.Rdata")

## Plot PS ----
ggplot(propensityScores, aes(x=ps, fill=group, color = group)) +
  geom_density(alpha=0.3) + theme_minimal() +
  labs(x = "Propensity scores", y = "Density", fill = "Cohort", color = "Cohort") +
  scale_fill_discrete(labels = c("Unvaccinated", "Vaccinated")) +
  scale_color_discrete(labels = c("Unvaccinated", "Vaccinated")) 

ggsave(width = 5, height = 4, dpi = 300,
       filename = "./BalancePlots/s01_PS.png")

print("STEP 6: get SMD to evaluate balance")

# =============================== BALANCE (SMD) ============================== # ----
## SMD before PS ----
print("SMD before PS")
# Get SMD for discrete covariates:
discrete_dat1 <- lasso_data %>% mutate(weight = 1)
discrete_dat2 <- individuals_id %>% select(person_id) %>% left_join(factor_table) %>% 
  pivot_wider(names_from = covariate, values_from = value, values_fill = 0) %>% mutate(weight = 1)
discrete_dat <- discrete_dat1 %>% inner_join(discrete_dat2)

discrete_smd_before <- compute_discrete_smd(discrete_dat,groupNames)

# Get SMD for continuous covariates:
continuous_dat <- individuals_id %>% select(person_id) %>%
  left_join(continuous_table) %>% pivot_wider(names_from = covariate, values_from = value, values_fill = 0) %>%
  inner_join(setGroups(cohort_id_groups, groupNames)) %>% mutate(weight = 1)

continuous_smd_before <- compute_continuous_smd(continuous_dat, groupNames)

# Get the overall SMD for  age, care site, and region
age_group_factor <- age_group_table %>% mutate(value = gsub("age_group_","",covariate), covariate = "age_group")
gp_factor        <- loc_table[grepl("care_site",loc_table$covariate), ] %>% mutate(value = gsub("care_site_","",covariate), covariate = "care_site")
region_factor    <- loc_table[grepl("region",loc_table$covariate), ] %>% mutate(value = gsub("region_","",covariate), covariate = "region")

overall_factor_table <-  age_group_factor %>% union_all(gp_factor) %>% union_all(region_factor)
overall_factor_dat <- individuals_id %>% select(person_id) %>% left_join(overall_factor_table) %>%
  inner_join(setGroups(cohort_id_groups, groupNames)) %>%
  pivot_wider(names_from = covariate, values_from = value, values_fill = NA) %>% mutate(weight = 1)

vars <- colnames(overall_factor_dat)[3:5]
w.data <- svydesign(ids = ~1, data = overall_factor_dat, weights = ~ weight) 
tabWeighted <- svyCreateTableOne(vars = vars, strata = "group", data = w.data, test = FALSE)
overall_factor_smd_before <- tibble(covariate = vars, smd = as.numeric(ExtractSmd(tabWeighted)))

# Table with all factors and their smd
balance_tab_before <- discrete_smd_before %>% select(covariate, smd) %>%
  union_all(continuous_smd_before %>% select(covariate, smd)) %>% union_all(overall_factor_smd_before)

print(paste0("Proportion of balanced (all) ", nrow(balance_tab_before %>% filter(smd <= 0.1 & smd >= -0.1))/nrow(balance_tab_before)))

save(list = c("balance_tab_before","continuous_smd_before","discrete_smd_before", "overall_factor_smd_before"), file = "RData/s01_FE_balance_before.Rdata")
write.csv(balance_tab_before,"./Data/s01_balance_table_before.csv")

## SMD after PS ----
print("SMD after PS")
# Do the same but with the weights (vax: 1-ps, unvax: ps)
discrete_dat       <- discrete_dat %>% inner_join(propensityScores) %>%
  mutate(weight = ifelse(group == "vaccinated", 1-ps, ps)) %>% select(-ps)
discrete_smd_after <- compute_discrete_smd(discrete_dat,groupNames)

continuous_dat       <- continuous_dat %>% inner_join(propensityScores) %>%
  mutate(weight = ifelse(group == "vaccinated", 1-ps, ps)) %>% select(-ps)
continuous_smd_after <- compute_continuous_smd(continuous_dat, groupNames)

overall_factor_dat   <- overall_factor_dat %>% inner_join(propensityScores) %>%
  mutate(weight = ifelse(group == "vaccinated", 1-ps, ps)) %>% select(-ps)
w.data <- svydesign(ids = ~1, data = overall_factor_dat, weights = ~ weight) 
tabWeighted <- svyCreateTableOne(vars = vars, strata = "group", data = w.data, test = FALSE)
overall_factor_smd_after <- tibble(covariate = vars, smd = as.numeric(ExtractSmd(tabWeighted)))

balance_tab_after <- discrete_smd_after %>% select(covariate, smd) %>%
  union_all(continuous_smd_after %>% select(covariate, smd)) %>% union_all(overall_factor_smd_after)

print(paste0("Proportion of balanced (all) after PS: ", nrow(balance_tab_after %>% filter(smd <= 0.1 & smd >= -0.1))/nrow(balance_tab_after)))

save(list = c("balance_tab_after","continuous_smd_after","discrete_smd_after", "overall_factor_smd_after"), file = "RData/s01_FE_balance_after.Rdata")

## Plot SMD ----
plotBalance <- balance_tab_before %>% select(covariate, smd_before = smd) %>% 
  inner_join(balance_tab_after %>% select(covariate, smd_after = smd) )

ggplot(plotBalance, aes(x = smd_before, y = smd_after)) + ylab("SMD after weighting") +
  xlab("SMD before weighting") +  geom_point(size = 1, col = rgb(0, 0, 1, alpha = 0.5)) +
  ylim(c(0,0.35)) + xlim(c(0,0.35)) +
  geom_vline(xintercept = 0.1, col = "grey") +
  geom_hline(yintercept = 0.1, col = "grey")

ggsave(width = 7, height = 4, dpi = 300,
       filename = "./BalancePlots/s01_smd.png")


## NCO ----
print("STEP 7: get NCO")
NCO <- read_csv(file=here("Data","NCO.csv"),show_col_types = FALSE) # File with the conditions chosen to assess nco

# Use the table with all extracted conditions to get the NCOs
NCO_dates <- NCO %>% rename(covariateId = ConceptId) %>% 
  inner_join(CO_patients %>% inner_join(CO_id) %>% select(person_id, covariateId,Event_date))

# Filter conditions after index date and get the first occurrence
data_NCO_aux <- individuals_cru %>% 
  left_join(NCO_dates) %>%
  filter(Event_date > index_date | is.na(Event_date)) %>%
  group_by(person_id,covariateId) %>%
  filter(Event_date == min(Event_date)) %>%
  ungroup() %>%
  distinct() %>%
  select(person_id, OutcomeName, Event_date) %>% mutate(OutcomeName = gsub(" ", "_", OutcomeName))

# data_NCO: person_id, index_date, group, next_vaccine, leave_date, death_date, gender, age, weight, date_NCO1, date_NCO2...
data_NCO <- individuals_cru %>% 
  left_join(data_NCO_aux %>% pivot_wider(names_from = OutcomeName, values_from = Event_date, values_fill = NA))

# Get names of the NCO conditions for the for loop
names <- colnames(data_NCO)[10:length(colnames(data_NCO))]
groupNames <- c("vaccinated", "unvaccinated")
# Save survival analysis results in these tables:
toPlot_NCO_before <- tibble(Outcome = names, RelativeRisk = NA, SD = NA)
toPlot_NCO_after  <- tibble(Outcome = names, RelativeRisk = NA, SD = NA)


for (nam in names) {
  # For each NCO condition we do fine and grey:
  data_nco_regression <- data_NCO[,c(colnames(data_NCO)[1:9], nam)]

  data_1 <- data_nco_regression[data_nco_regression$group == "vaccinated",]
  data_2 <- data_nco_regression[data_nco_regression$group == "unvaccinated",]
  minimTime_1 <- pmin(data_1$death_date, data_1$leave_date,
                      data_1[[nam]], na.rm = TRUE)
  minimTime_2 <- pmin(data_2$death_date, data_2$leave_date,
                      data_2[[nam]], data_2$next_vaccine,na.rm = TRUE)
  data_1 <- data_1 %>% mutate(time = as.integer(minimTime_1 - index_date))
  data_1$event <- 0
  data_1$event[minimTime_1 == data_1[[nam]]]     <- 1 # Event of interest
  data_1$event[minimTime_1 == data_1$death_date] <- 2 # Competing event
  
  data_2 <- data_2 %>% mutate(time = as.integer(minimTime_2 - index_date))
  data_2$event <- 0
  data_2$event[minimTime_2 == data_2[[nam]]]     <- 1 # Event of interest
  data_2$event[minimTime_2 == data_2$death_date] <- 2 # Competing event
  
  data_nco_regression <- data_1 %>% union_all(data_2) %>% 
    select(person_id, group, time, event, weight)
  
  data_nco_regression$weight <- 1
  data_nco_regression$event  <- factor(data_nco_regression$event, levels = c(0,1,2), labels = c("Censored","PACS","Death"))
  data_nco_regression$group  <- as.factor(data_nco_regression$group) %>% relevel(group, ref = groupNames[2])
  
  # Results from Fine and Gray without weights (weight = 1)
  aux_before <- fgNCO(data_nco_regression)
  
  # Save results into the dataframe
  toPlot_NCO_before[toPlot_NCO_before$Outcome == nam, 2] <- aux_before[1]
  toPlot_NCO_before[toPlot_NCO_before$Outcome == nam, 3] <- aux_before[2]
  
  # Results from Fine and Gray with weights 
  data_nco_regression_after <- data_nco_regression %>% inner_join(propensityScores) %>%
    mutate(weight = ifelse(group == "vaccinated", 1-ps, ps)) %>% select(-ps)
  
  aux_after <- fgNCO(data_nco_regression_after)
  
  toPlot_NCO_after[toPlot_NCO_after$Outcome == nam, 2] <- aux_after[1]
  toPlot_NCO_after[toPlot_NCO_after$Outcome == nam, 3] <- aux_after[2]
}

save(list = c("toPlot_NCO_before","data_NCO","toPlot_NCO_after"), file = "./RData/s01_FE_NCO.Rdata")
write.csv(toPlot_NCO_before,"./Data/s01_NCO_before.csv")
write.csv(toPlot_NCO_after,"./Data/s01_NCO_after.csv")


# Plot NCO ----
plotCiCalibrationEffect_NMB(log(toPlot_NCO_before$RelativeRisk), toPlot_NCO_before$SD, rep(0,nrow(toPlot_NCO_before)), 
                            fileName = "./BalancePlots/s01_NCO_before.png")

plotCiCalibrationEffect_NMB(log(toPlot_NCO_after$RelativeRisk), toPlot_NCO_after$SD, rep(0,nrow(toPlot_NCO_after)), 
                            fileName = "./BalancePlots/s01_NCO_after.png")

