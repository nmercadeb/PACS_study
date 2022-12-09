## ----------------------------- LOAD PACKAGES ------------------------------ ##
library("lubridate")
library("tidyr")
library("glmnet")
library("readr")
library("survey")
library("tableone")
library("ggplot2")
library("random")

## ------------------------------- FUNCTIONS -------------------------------- ##
## getIndividualTables
getIndividualTabs      <- function(tab_id, tab_patients, population_tab, tabName, numOfWindows, value1) {
  if (numOfWindows == 2) {
    aux_table <- tab_patients %>%
      inner_join(population_tab, by = "person_id") %>% 
      filter(index_date-Era_start_date <= days(180) & index_date-Era_start_date > days(0)) %>%
      mutate(window = if_else(index_date-Era_start_date >= days(1) & index_date-Era_start_date <= days(30), "1to30", "31to180")) %>%
      inner_join(tab_id, by = "FeatureExtractionId") %>%
      mutate(covariate = gsub(" ","_", paste(tabName,covariateName,window, sep = "_"))) %>%
      select(person_id, covariate) %>%
      distinct() %>%
      mutate(value = 1)
  } else if (numOfWindows == 3) {
    aux_table <- tab_patients %>%
      inner_join(population_tab, by = "person_id") %>% 
      filter(index_date-Event_date > days(0)) %>%
      mutate(window = NA)
    aux_table$window[aux_table$index_date-aux_table$Event_date >= days(1) & aux_table$index_date-aux_table$Event_date <= days(30)]    <- "1to30"
    aux_table$window[aux_table$index_date-aux_table$Event_date >= days(31) & aux_table$index_date-aux_table$Event_date <= days(180)]  <- "31to180"
    aux_table$window[aux_table$index_date-aux_table$Event_date >= days(181)] <- "+181"
    
    if (value1)
      aux_table <- aux_table %>%
      inner_join(tab_id, by = "FeatureExtractionId") %>%
      mutate(covariate = gsub(" ","_", paste(tabName,covariateName,window, sep = "_"))) %>%
      select(person_id, covariate) %>%
      distinct() %>%
      mutate(value = 1) 
    else {
      aux_table <- aux_table %>%
        inner_join(tab_id, by = "FeatureExtractionId") %>%
        mutate(covariate = gsub(" ","_", paste(tabName,covariateName,window, sep = "_"))) %>%
        select(person_id, covariate) %>%
        group_by(person_id,covariate) %>%
        count() %>%
        rename(value = n)
    }
  } else if (numOfWindows == 0) {
    if (value1) {
      aux_table <- tab_patients %>%
        inner_join(population_tab, by = "person_id") %>% 
        inner_join(tab_id, by = c("FeatureExtractionId", "value")) %>%
        select(-value) %>%
        mutate(covariate = gsub(" ","_", paste(tabName,valueName, sep = "_")), value = 1) %>%
        select(person_id, covariate, value) %>%
        distinct()
      
      if(tabName == "location") {
        aux_table <- aux_table %>%
          mutate(covariate = ifelse(grepl("care_site",covariate),gsub("location_","",covariate),gsub("location_","region_",covariate)))
        aux_table <- aux_table[!grepl("region_NA", aux_table$covariate),]
        aux_table <- aux_table[!grepl("region_Northern_Ireland", aux_table$covariate),]
      }
    } else {
      aux_table <- tab_patients %>% inner_join(population_tab, by = "person_id") %>% filter(index_date-Event_date > days(0)) %>%
        inner_join(tab_id, by = "FeatureExtractionId") %>%
        mutate(covariate = gsub(" ","_", paste(tabName,covariateName,"allTime", sep = "_"))) %>% select(person_id, covariate) %>%
        group_by(person_id,covariate) %>% count() %>% rename(value = n)
    }
  }
  return(aux_table)
}

## setGroups 
setGroups              <- function(cohort_ids, groupNames) {
  id1 <- cohort_ids[1]
  id2 <- cohort_ids[2]
  aux_tab_1 <- cohorts_db %>%
    filter(cohort_definition_id == id1) %>%
    select(person_id = subject_id) %>%
    collect()
  aux_tab_1$group <- groupNames[1]
  aux_tab_2 <- cohorts_db %>%
    filter(cohort_definition_id == id2) %>%
    select(person_id = subject_id) %>%
    collect()
  aux_tab_2$group <- groupNames[2]
  return(aux_tab_1 %>% union_all(aux_tab_2))
}

## SMD for discrete variables
compute_discrete_smd   <- function(data_bin,groupNames){
  
  namt <- names(data_bin)
  namt <- namt[!(namt %in% c("person_id","group","weight"))]
  
  data_g1 <- data_bin %>% filter(group == groupNames[1])
  mean1   <- data_g1 %>% summarise_at(namt, Hmisc::wtd.mean, weights = data_g1$weight, normwt = FALSE)
  var1    <- data_g1 %>% summarise_at(namt, Hmisc::wtd.var,  weights = data_g1$weight)
  
  data_g2 <- data_bin %>% filter(group == groupNames[2])
  mean2   <- data_g2 %>% summarise_at(namt, Hmisc::wtd.mean, weights = data_g2$weight)
  var2    <- data_g2 %>% summarise_at(namt, Hmisc::wtd.var,  weights = data_g2$weight)
  
  table1  <- rbind(mean1,var1,mean2,var2)
  table1  <- tibble(covariate = names(table1), mean1 = t(table1[1,])[,1], var1 = t(table1[2,])[,1], mean2 = t(table1[3,])[,1], var2 = t(table1[4,])[,1]) %>%
    mutate(smd = abs(mean1-mean2)/sqrt((var1+var2)/2))
  
  return(table1)
}

## SMD for continuous variables
compute_continuous_smd <- function(data_cont, groupNames){
  
  namt <- names(data_cont)
  namt <- namt[!(namt %in% c("person_id","group","weight"))]
  
  data_g1 <- data_cont %>% filter(group == groupNames[1])
  mean1   <- data_g1 %>% summarise_at(namt, Hmisc::wtd.mean, weights = data_g1$weight)
  var1    <- data_g1 %>% summarise_at(namt, Hmisc::wtd.var,  weights = data_g1$weight)
  
  data_g2 <- data_cont %>% filter(group == groupNames[2])
  mean2   <- data_g2 %>% summarise_at(namt, Hmisc::wtd.mean, weights = data_g2$weight)
  var2    <- data_g2 %>% summarise_at(namt, Hmisc::wtd.var,  weights = data_g2$weight)
  
  table1  <- rbind(mean1,var1,mean2,var2)
  table1  <- tibble(covariate = names(table1), mean1 = t(table1[1,])[,1], var1 = t(table1[2,])[,1], mean2 = t(table1[3,])[,1], var2 = t(table1[4,])[,1]) %>%
    mutate(smd = abs(mean1-mean2)/sqrt(var1+var2))
  
  return(table1)
}

## Fine and Gray for NCO
fgNCO                  <- function(table) {
  fg_data <- finegray(Surv(time, event) ~ ., data=table,  weights = weight) 
  fg_regression <- coxph(Surv(fgstart, fgstop, fgstatus) ~ group, weight=fgwt, data=fg_data)
  coef <- summary(fg_regression)$coefficients[c(2,3)]
  return(coef)
} 

# logRrtoSE + plotCiCalibrationEffect_NMB,  used to plot NCO
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
    ggplot2::facet_grid(. ~ Group) +
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
    ggsave(width = 7, height = 4, dpi = 600,
           filename = fileName)
  }
  return(plot)
}

## ------------------------------- VARIABLES -------------------------------- ##
# Cohorts_ID
AZvaccine_id       <- 2
PFvaccine_id       <- 3
covid_diagonsis_id <- 4
test_covid_any_id  <- 175
test_covid_pcr_id  <- 176
unexposed_id       <- 177
exposed_id         <- 186


cohort_id_groups   <- c(exposed_id,unexposed_id)

# List of individuals 
individuals_id <- cohorts_db %>%
  filter(cohort_definition_id %in% c(exposed_id, unexposed_id)) %>%
  select(subject_id, cohort_start_date) %>%
  rename("person_id" = "subject_id", "index_date" = "cohort_start_date" ) %>%
  compute()
list_id <- individuals_id %>%
  select(person_id) %>%
  compute()

IDvac_collected <- individuals_id %>% 
  inner_join(cohorts_db %>%
               filter(cohort_definition_id == exposed_id) %>%
               select(subject_id, cohort_start_date) %>%
               rename("person_id" = "subject_id", "index_date" = "cohort_start_date" )) %>% collect()
UNvac_collected <- individuals_id %>% 
  inner_join(cohorts_db %>%
               filter(cohort_definition_id == unexposed_id) %>%
               select(subject_id, cohort_start_date) %>%
               rename("person_id" = "subject_id", "index_date" = "cohort_start_date" )) %>% collect()

# Re-do sorting:
nbins <- table(IDvac_collected$index_date)
dates <- names(nbins)
names(nbins) <- NULL
prob  <- nbins/sum(nbins)
set.seed(randomNumbers(n = 1, col = 1))
UNvac_collected$index_date <- as.Date(sample(dates, nrow(UNvac_collected), replace = TRUE, prob=prob))

individuals_id <- IDvac_collected %>% union_all(UNvac_collected)
write.csv(individuals_id, file = "./Data/s01_individuals_PF.csv")

# Variables of interest
date1      <- as.Date("2022-01-01") # For age
groupNames <- c("vaccinated", "unvaccinated")
numOfPeople <- nrow(individuals_id) # Total number of individuals
personsPerGroup <- individuals_id %>% 
  inner_join(setGroups(cohort_id_groups,groupNames)) %>%
  mutate(present = 1) %>%
  group_by(group) %>%
  summarize(population = sum(present))
numExposed   <- personsPerGroup$population[personsPerGroup$group == groupNames[1]]
numUnexposed <- personsPerGroup$population[personsPerGroup$group == groupNames[2]]










