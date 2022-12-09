# ============================================================================ #
#              Feature Extraction - SMD - PS  (Vaccine strata: AZ)             #
#                                Núria Mercadé                                 #
#                                 28-10-2022                                   #
# ============================================================================ #

## Connect to database
source("./AURUM_CDM_connection.R")

## Packages, variables, tables and functions needed
source("./s01_forFeatureExtraction_PF.R")

# =========================== TABLES OF INTEREST ============================= # ----
print("Get NCO for AZ")
load(file = "./RData/s01_CO_forNCO.RData")
CO_patients<- CO_patients %>% inner_join(individuals_id %>% select(person_id))

load(file = "RData/s01_FE_individualTabs_2.Rdata")
c("drug_table", "co_table", "vi_table", "test_table",
  "age_group_table", "gend_table", "loc_table", "age_table", 
  "obs_years_table")
drug_table      <- drug_table %>% inner_join(individuals_id %>% select(person_id))
co_table        <- co_table %>% inner_join(individuals_id %>% select(person_id))
vi_table        <- vi_table %>% inner_join(individuals_id %>% select(person_id))
test_table      <- test_table %>% inner_join(individuals_id %>% select(person_id))
age_group_table <- age_group_table %>% inner_join(individuals_id %>% select(person_id))
gend_table      <- gend_table %>% inner_join(individuals_id %>% select(person_id))
loc_table       <- loc_table %>% inner_join(individuals_id %>% select(person_id))
age_table       <- age_table %>% inner_join(individuals_id %>% select(person_id))
obs_years_table <- obs_years_table %>% inner_join(individuals_id %>% select(person_id))

# Get most frequent covariates
discrete_table   <- drug_table %>% union_all(co_table) 
continuous_table <- vi_table %>% union_all(test_table) %>% union_all(age_table) %>% union_all(obs_years_table) %>% ungroup()
factor_table     <- age_group_table %>% union_all(gend_table) %>% union_all(loc_table)

minimum.freq.covariates <- 0.005
min_cov <- floor(nrow(individuals_id)*minimum.freq.covariates)
exclude <- discrete_table %>% group_by(covariate) %>% tally() %>% filter(n<=min_cov) %>% select(covariate)

print(paste0("DISCRETE: Exclude ", exclude %>% tally() %>% as.numeric(), " covariates out of ", length(unique(discrete_table$covariate))))

discrete_table <- discrete_table %>% anti_join(exclude) 

save(list = c("discrete_table", "continuous_table", "factor_table"), 
     file = "RData/s01_FE_prep_PF.Rdata")

print('Start lasso regression for PS')
# ============================= PROPENSITY SCORE ============================= # ----
## LASSO REGRESSION ---------------------
# Data for LASSO regression --> ALL except those compulsory
# discrete table: ALL
# continuous table: none (all are compulsory)
# factor table: none (all are complusory)

# lasso_data <- individuals_id %>% select(person_id) %>% inner(discrete_table) %>%
#   pivot_wider(names_from = covariate, values_from = value, values_fill = 0) %>%
#   left_join(setGroups(cohort_id_groups, groupNames)) %>% mutate(group = as.factor(group)) %>%
#   mutate(group = relevel(group, ref = groupNames[2]))
lasso_data_aux <- discrete_table %>%
  pivot_wider(names_from = covariate, values_from = value, values_fill = 0) 
lasso_data <- individuals_id %>% select(person_id) %>% left_join(lasso_data_aux) %>%
  left_join(setGroups(cohort_id_groups, groupNames)) %>% mutate(group = as.factor(group)) %>% 
  mutate(group = relevel(group, ref = groupNames[2]))
lasso_data[is.na(lasso_data)] <- 0

X <- as.matrix(lasso_data %>% select(-person_id,-group))
Y <- lasso_data$group

print("Data for lasso obtained, starting regression")

lambdas   <- 10^seq(2, -3, by = -.1)
lasso_reg <- cv.glmnet(X, Y, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5, family="binomial")
coef.lasso_reg <- coef(lasso_reg, s = lasso_reg$lambda.1se)
var.lasso_reg  <- names(coef.lasso_reg[(coef.lasso_reg[,1]!=0),1])
var.lasso_reg  <- gsub("`", "", var.lasso_reg)

save(list = c("lasso_data","var.lasso_reg"), file = "RData/s01_FE_lasso_data_PF.Rdata")
save(list = c("lasso_reg"), file = "RData/s01_FE_lasso_regression_PF.Rdata")

print("Get PROPENSITY SCORES for PF")
## PROPENSITY SCORES ---------------------
nam             <- colnames(lasso_data)
id_col          <- c(1,length(nam),which(nam %in% var.lasso_reg)) # Take person_id, group and lasso resulting features

linearIndependent <- tibble(covariate = c(unique(factor_table$covariate[grep("care_site",factor_table$covariate)])[1],
                                          unique(factor_table$covariate[grep("region",factor_table$covariate)])[1], 
                                          unique(factor_table$covariate[grep("age_g",factor_table$covariate)])[1],
                                          unique(factor_table$covariate[grep("gender",factor_table$covariate)])[1]))

compulsory_features <- individuals_id %>% select(person_id) %>% 
  left_join(continuous_table %>% union_all(factor_table %>% anti_join(linearIndependent))) %>% 
  pivot_wider(names_from = covariate, values_from = value, values_fill = 0) 

data_base <- lasso_data[,id_col] %>% 
  left_join(compulsory_features) %>% left_join(individuals_id %>% mutate(index_date = as.numeric(index_date)))
data_base$age2  <- data_base$age^2

data_base_subset <- data_base[sample(nrow(data_base), 200000, replace = FALSE),]

print("Data for PS obtained, starting regression + prediction")
model_base    <- glm(formula = group ~ . - person_id, family = 'binomial', data = data_base_subset)
PS     <- predict(model_base,data_base,type = "response")

propensityScores <- data_base %>% select(person_id,group) %>% mutate(ps = PS)

save(list=c("model_base"), file = "RData/s01_FE_PS_model_base_PF.Rdata")
save(list=c("compulsory_features","data_base","PS", "propensityScores"), file = "RData/s01_FE_PS_data_PF.Rdata")
write.csv(propensityScores, file = "./Data/s01_ps_PF.csv")

## Plot PS ----
ggplot(propensityScores, aes(x=ps, fill=group, color = group)) +
  geom_density(alpha=0.3) + theme_minimal() +
  labs(x = "Propensity scores", y = "Density", fill = "Cohort", color = "Cohort") +
  scale_fill_discrete(labels = c("Unvaccinated", "Vaccinated")) +
  scale_color_discrete(labels = c("Unvaccinated", "Vaccinated")) 

ggsave(width = 7, height = 4, dpi = 300,
       filename = "./BalancePlots/s01_PS_PF.png")

print("STEP 6: get SMD to evaluate balance")

# =============================== BALANCE (SMD) ============================== # ----
## SMD before PS ----
print("SMD before PS")
discrete_dat1 <- lasso_data %>% mutate(weight = 1)
discrete_dat2 <- individuals_id %>% select(person_id) %>% left_join(factor_table) %>% 
  pivot_wider(names_from = covariate, values_from = value, values_fill = 0) %>% mutate(weight = 1)
discrete_dat <- discrete_dat1 %>% inner_join(discrete_dat2)

discrete_smd_before <- compute_discrete_smd(discrete_dat,groupNames)

continuous_dat <- individuals_id %>% select(person_id) %>%
  left_join(continuous_table) %>% pivot_wider(names_from = covariate, values_from = value, values_fill = 0) %>%
  inner_join(setGroups(cohort_id_groups, groupNames)) %>% mutate(weight = 1)

continuous_smd_before <- compute_continuous_smd(continuous_dat, groupNames)

# Age, gp and region overall
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

balance_tab_before <- discrete_smd_before %>% select(covariate, smd) %>%
  union_all(continuous_smd_before %>% select(covariate, smd)) %>% union_all(overall_factor_smd_before)

print(paste0("Proportion of balanced (all) ", nrow(balance_tab_before %>% filter(smd <= 0.1 & smd >= -0.1))/nrow(balance_tab_before)))

save(list = c("balance_tab_before","continuous_smd_before","discrete_smd_before", "overall_factor_smd_before"), file = "RData/s01_FE_balance_before_PF.Rdata")
write.csv(balance_tab_before,"./Data/s01_balance_table_before_PF.csv")

## SMD after PS ----
print("SMD after PS")
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

save(list = c("balance_tab_after","continuous_smd_after","discrete_smd_after", "overall_factor_smd_after"), file = "RData/s01_FE_balance_after_PF.Rdata")

## Plot SMD ----
plotBalance <- balance_tab_before %>% select(covariate, smd_before = smd) %>% 
  inner_join(balance_tab_after %>% select(covariate, smd_after = smd) )

ggplot(plotBalance, aes(x = smd_before, y = smd_after)) + ylab("SMD after weighting") +
  xlab("SMD before weighting") +  geom_point(size = 1, col = rgb(0, 0, 1, alpha = 0.5)) +
  ylim(c(0,0.35)) + xlim(c(0,0.35)) +
  geom_vline(xintercept = 0.1, col = "grey") +
  geom_hline(yintercept = 0.1, col = "grey")

ggsave(width = 7, height = 4, dpi = 300,
       filename = "./BalancePlots/s01_smd_PF.png")


print("get NCO")
## NCO ----
NCO <- read_csv(file=here("Data","NCO.csv"),show_col_types = FALSE)
load("./RData/s01_individuals_weight.Rdata")
individuals_cru <- individuals_cru %>% select(-index_date) %>% inner_join(individuals_id, by = "person_id")

NCO_dates <- NCO %>% rename(covariateId = ConceptId) %>% 
  inner_join(CO_patients %>% inner_join(CO_id) %>% select(person_id, covariateId,Event_date))

data_NCO_aux <- individuals_cru %>% 
  left_join(NCO_dates) %>%
  filter(Event_date > index_date | is.na(Event_date)) %>%
  group_by(person_id,covariateId) %>%
  filter(Event_date == min(Event_date)) %>%
  ungroup() %>%
  distinct() %>%
  select(person_id, OutcomeName, Event_date) %>% mutate(OutcomeName = gsub(" ", "_", OutcomeName))
data_NCO <- individuals_cru %>% 
  left_join(setGroups(cohort_id_groups, groupNames)) %>% mutate(group = as.factor(group)) %>% 
  mutate(group = relevel(group, ref = groupNames[2])) %>%
  left_join(data_NCO_aux %>% pivot_wider(names_from = OutcomeName, values_from = Event_date, values_fill = NA))


names <- colnames(data_NCO)[10:length(colnames(data_NCO))]
groupNames <- c("vaccinated", "unvaccinated")
toPlot_NCO_before <- tibble(Outcome = names, RelativeRisk = NA, SD = NA)
toPlot_NCO_after <- tibble(Outcome = names, RelativeRisk = NA, SD = NA)


for (nam in names) {
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
  
  aux_before <- fgNCO(data_nco_regression)
  
  toPlot_NCO_before[toPlot_NCO_before$Outcome == nam, 2] <- aux_before[1]
  toPlot_NCO_before[toPlot_NCO_before$Outcome == nam, 3] <- aux_before[2]
  
  data_nco_regression_after <- data_nco_regression %>% inner_join(propensityScores) %>%
    mutate(weight = ifelse(group == "vaccinated", 1-ps, ps)) %>% select(-ps)
  
  aux_after <- fgNCO(data_nco_regression_after)
  
  toPlot_NCO_after[toPlot_NCO_after$Outcome == nam, 2] <- aux_after[1]
  toPlot_NCO_after[toPlot_NCO_after$Outcome == nam, 3] <- aux_after[2]
}

save(list = c("toPlot_NCO_before","data_NCO","toPlot_NCO_after"), file = "./RData/s01_FE_NCO_PF.Rdata")
write.csv(toPlot_NCO_before,"./Data/s01_NCO_before_PF.csv")
write.csv(toPlot_NCO_after,"./Data/s01_NCO_after_PF.csv")

individuals_PF <- individuals_cru %>% inner_join(propensityScores) %>%
  mutate(weight = ifelse(group == "vaccinated", 1-ps, ps))
write.csv(individuals_PF %>% select(group, ps), file = "./Data/s01_individuals_PF.csv")
individuals_PF <- individuals_PF %>% select(-ps)
save("individuals_PF", file = "./RData/s01_individuals_PF.Rdata")

# Calibrate NCO ----
model <- fitSystematicErrorModel(log(toPlot_NCO_after$RelativeRisk), toPlot_NCO_after$SD, rep(0,nrow(toPlot_NCO_after)))
result <- calibrateConfidenceInterval(log(toPlot_NCO_after$RelativeRisk), toPlot_NCO_after$SD, model,ciWidth = 0.95)
toPlot_NCO_calibrated <- toPlot_NCO_after[,1] %>% cbind(exp(result))
colnames(toPlot_NCO_calibrated)[2:5] <- c("RelativeRisk","low95", "upper95","SD")
write.csv(toPlot_NCO_calibrated, file = "./Data/s01_NCO_calibrated_PF.csv")


# Plot NCO ----
logRrtoSE <- function(logRr, alpha, mu, sigma) {
  phi <- (mu - logRr)^2 / qnorm(alpha / 2)^2 - sigma^2
  phi[phi < 0] <- 0
  se <- sqrt(phi)
  return(se)
}
plotCiCalibrationEffect_NMB <- function(logRr,
                                        seLogRr,
                                        trueLogRr,
                                        legacy = FALSE,
                                        model = NULL,
                                        xLabel = "Relative risk",
                                        title,
                                        fileName = NULL) {
  alpha <- 0.05
  if (is.null(model)) {
    model <- fitSystematicErrorModel(
      logRr = logRr,
      seLogRr = seLogRr,
      trueLogRr = trueLogRr,
      estimateCovarianceMatrix = FALSE,
      legacy = legacy
    )
  } else {
    legacy <- (names(model)[3] == "logSdIntercept")
  }
  d <- data.frame(
    logRr = logRr,
    seLogRr = seLogRr,
    trueLogRr = trueLogRr,
    trueRr = exp(trueLogRr),
    logCi95lb = logRr + qnorm(0.025) * seLogRr,
    logCi95ub = logRr + qnorm(0.975) * seLogRr
  )
  d <- d[!is.na(d$logRr), ]
  d <- d[!is.na(d$seLogRr), ]
  if (nrow(d) == 0) {
    return(NULL)
  }
  d$Group <- as.factor(d$trueRr)
  d$Significant <- d$logCi95lb > d$trueLogRr | d$logCi95ub < d$trueLogRr
  
  temp1 <- aggregate(Significant ~ trueRr, data = d, length)
  temp2 <- aggregate(Significant ~ trueRr, data = d, mean)
  temp1$nLabel <- paste0(formatC(temp1$Significant, big.mark = ","), " estimates")
  temp1$Significant <- NULL
  temp2$meanLabel <- paste0(
    formatC(100 * (1 - temp2$Significant), digits = 1, format = "f"),
    "% of CIs includes ",
    temp2$trueRr
  )
  temp2$Significant <- NULL
  dd <- merge(temp1, temp2)
  
  breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
  theme <- ggplot2::element_text(colour = "#000000", size = 10)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 10, hjust = 1)
  
  d$Group <- paste("True", tolower(xLabel), "=", d$trueRr)
  dd$Group <- paste("True", tolower(xLabel), "=", dd$trueRr)
  
  x <- seq(log(0.1), log(10), by = 0.01)
  calBounds <- data.frame()
  for (i in 1:nrow(dd)) {
    mu <- model[1] + model[2] * log(dd$trueRr[i])
    if (legacy) {
      sigma <- exp(model[3] + model[4] * log(dd$trueRr[i]))
    } else {
      sigma <- model[3] + model[4] * abs(log(dd$trueRr[i]))
    }
    calBounds <- rbind(
      calBounds,
      data.frame(
        logRr = x,
        seLogRr = logRrtoSE(x, alpha, mu, sigma),
        Group = dd$Group[i]
      )
    )
  }
  plot <- ggplot2::ggplot(d, ggplot2::aes(x = .data$logRr, y = .data$seLogRr)) +
    ggplot2::geom_vline(xintercept = log(breaks), colour = "#AAAAAA", lty = 1, size = 0.5) +
    ggplot2::geom_area(
      fill = rgb(1, 0.5, 0, alpha = 0.5),
      color = rgb(1, 0.5, 0),
      size = 1,
      alpha = 0.5, data = calBounds
    ) +
    ggplot2::geom_abline(ggplot2::aes(intercept = (-log(.data$trueRr)) / qnorm(0.025), slope = 1 / qnorm(0.025)), colour = rgb(0, 0, 0), linetype = "dashed", size = 1, alpha = 0.5, data = dd) +
    ggplot2::geom_abline(ggplot2::aes(intercept = (-log(.data$trueRr)) / qnorm(0.975), slope = 1 / qnorm(0.975)), colour = rgb(0, 0, 0), linetype = "dashed", size = 1, alpha = 0.5, data = dd) +
    ggplot2::geom_point(
      shape = 16,
      size = 2,
      alpha = 0.5,
      color = rgb(0, 0, 0.8)
    ) +
    ggplot2::geom_hline(yintercept = 0) +
    # ggplot2::geom_label(x = log(0.15), y = 0.95, alpha = 1, hjust = "left", ggplot2::aes(label = .data$nLabel), size = 3.5, data = dd) +
    # ggplot2::geom_label(x = log(0.15), y = 0.8, alpha = 1, hjust = "left", ggplot2::aes(label = .data$meanLabel), size = 3.5, data = dd) +
    ggplot2::scale_x_continuous(xLabel, limits = log(c(0.1, 10)), breaks = log(breaks), labels = breaks) +
    ggplot2::scale_y_continuous("Standard Error") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    # ggplot2::facet_grid(. ~ Group) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = themeRA,
      axis.text.x = theme,
      axis.title = theme,
      legend.key = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      strip.text.x = theme,
      strip.text.y = theme,
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    )
  
  if (!is.null(fileName)) {
    ggsave(width = 5, height = 4, dpi = 600,
           filename = fileName)
  }
  return(plot)
}

plotCiCalibrationEffect_NMB(log(toPlot_NCO_before$RelativeRisk), toPlot_NCO_before$SD, rep(0,nrow(toPlot_NCO_before)), 
                            fileName = "./BalancePlots/s01_NCO_before_PF.png")

plotCiCalibrationEffect_NMB(log(toPlot_NCO_after$RelativeRisk), toPlot_NCO_after$SD, rep(0,nrow(toPlot_NCO_after)), 
                            fileName = "./BalancePlots/s01_NCO_after_PF.png")

