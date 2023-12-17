# COVID-19 vaccine effectiveness against cardiac and thromboembolic complications following SARS-CoV-2 infection
Codes listed here are for the first staggered cohort study to assess the vaccine effectivness against Post Acute Covid-19 Sequelae. Notice that they assume that some cohorts are already instanciated cohorts. .JSON files are given as well as .csv with the cohort_definition_id for each.

## Briew description
Scripts:
  - AURUM_CDM_connection.R : These code is executed within the other scripts to stablish connection to the database and load some libraries.
  - outcome_cohorts.R: in this script the outcome cohorts are created and instanciated to the results schema in the database.
  - s01_dataAnalysis.R: this code does all the analysis once the data is obtained.
  - s01_featureExtraction.R: this code extract covariates for the propensity scores and evaluated the balance with ASMD and NCO before and afetr the        
                             OW weighting.
  - s01_forFeatureExtraction.R: this code is executed within the 's01_featureExtraction.R' script and contains libraries, functions and some 
                                variables to run the main script.
  - s01_featureExtraction_PF.R / s01_forFeatureExtraction_PF.R: same but the vaccination cohort just includes AstraZeneca vaccine  
  - s01_featureExtraction_PF.R / s01_forFeatureExtraction_PF.R: same but the vaccination cohort just includes Pfizer vaccine    
  
 Folders: 
  - Data: contains 'NCO.csv' which is a list of condition occurrences for evaluating the Negative Control Outcomes, and 'estudi_pacs_03.csv' which 
    contains the names of the table for each cohort_definition_id.
  - BalancePlots: empty, it will be filled with plots created in the featureExtraction scripts.
  - s01_RiskEstimates: contains multiple empty folders that will be filled with '.csv' files created in the script 's01_dataAnalysi' and contain the
                       sHR, 95%CI, and SE.
  - cohorts: cohorts as .JSON files.
 
  ### Considerations
  - The script 's01_forFeatureExtraction' saves the population weighted that will be used later for 's01_forFeatureExtraction_AZ' and 
    's01_forFeatureExtraction_PF'.
  - I would recomend to instanciate the cohorts following the legend in 'estudi_pacs_03.csv'.
  - Once the .JSON cohorts are instanciated, you can instanciate both population and outcome cohorts with the scripts for it. 
