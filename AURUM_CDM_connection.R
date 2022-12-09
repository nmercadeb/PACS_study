# ============================================================================ #
#                                CDM CONNECTION                                #
#                                Núria Mercadé                                 #
#                                 10-10-2022                                   #
# ============================================================================ #


# -------------------------- Connect to the database ------------------------- #
library(CDMConnector)
library(DatabaseConnector)
library("here") 
library("DBI") 
library("RPostgres") # Set connection details
library("dplyr")
library("dbplyr")
library("EmpiricalCalibration")

server     <- Sys.getenv("DB_SERVER_name_database") 
user       <- Sys.getenv("DB_USER") 
password   <- Sys.getenv("DB_PASSWORD") 
port       <- Sys.getenv("DB_PORT") 
host       <- Sys.getenv("DB_HOST") 
server_dbi <- Sys.getenv("DB_SERVER_DBI_name_database")

# Connect to the database 
con <- DBI::dbConnect(RPostgres::Postgres(), dbname = server_dbi, port = port, 
                      host = host, user = user, password = password)

# Some consideration about our database:
targetDialect              <- "postgresql"
cdm_database_schema        <- "public"
vocabulary_database_schema <- "public"
results_database_schema    <- "results"
cohort_tab_name            <- "estudi_pacs_03"

# Connect to cdm
cdm <- cdm_from_con(con, cdm_schema = cdm_database_schema,
                    # cdm_tables = c("person"),
                    write_schema    = results_database_schema,
                    cohort_tables   = cohort_tab_name)

# Databases that will be used       
person_db             <- cdm$person
concept_db            <- cdm$concept
care_site_db          <- cdm$care_site
location_db           <- cdm$location
death_db              <- cdm$death
observation_period_db <- cdm$observation_period
cohorts_db            <- cdm[[cohort_tab_name]]


# marti <- "lcve_mc_index"
# cdm <- cdm_from_con(con, cdm_schema = cdm_database_schema, 
#                     write_schema    = results_database_schema,
#                     cohort_tables   = marti)
# cohorts_study <- cdm$lcve_mc_index %>% inner_join(tibble(cohort_definition_id = 1:8), by = "cohort_definition_id", copy = TRUE) %>%
#   mutate(cohort_definition_id = ifelse(cohort_definition_id == 1,177, 
#                                        ifelse(cohort_definition_id == 2,178, 
#                                               ifelse(cohort_definition_id == 3,179, 
#                                                      ifelse(cohort_definition_id == 4,180, 
#                                                             ifelse(cohort_definition_id == 5,181, 
#                                                                    ifelse(cohort_definition_id == 6,182, 
#                                                                           ifelse(cohort_definition_id == 7,183, 
#                                                                                  ifelse(cohort_definition_id == 8,184, NA)))))))))
 
# marti <- "lcve_mc_index"
# cdm <- cdm_from_con(con, cdm_schema = cdm_database_schema,
#                     write_schema    = results_database_schema,
#                     cohort_tables   = marti)
# cohorts_study <- cdm$lcve_mc_index %>% inner_join(tibble(cohort_definition_id = c(9,10,13,14,17,18,21,22)), by = "cohort_definition_id", copy = TRUE) %>%
#   mutate(cohort_definition_id = ifelse(cohort_definition_id == 9,185,
#                                        ifelse(cohort_definition_id == 10,186,
#                                               ifelse(cohort_definition_id == 13,187,
#                                                      ifelse(cohort_definition_id == 14,188,
#                                                             ifelse(cohort_definition_id == 17,189,
#                                                                    ifelse(cohort_definition_id == 18,190,
#                                                                           ifelse(cohort_definition_id == 21,191,
#                                                                                  ifelse(cohort_definition_id == 22,192, NA)))))))))
# 
# library("SqlUtilities")
# appendPermanent(x = cohorts_study, name = 'estudi_pacs_03', schema = "results")


# connectionDetails <-DatabaseConnector::downloadJdbcDrivers("postgresql",here::here())
# connectionDetails <-DatabaseConnector::createConnectionDetails(dbms = "postgresql",
#                                                                server = server,
#                                                                user = user,
#                                                                password = password,
#                                                                port = port,
#                                                                pathToDriver = here::here())






