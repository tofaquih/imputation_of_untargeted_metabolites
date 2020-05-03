[![DOI](https://zenodo.org/badge/235815834.svg)](https://zenodo.org/badge/latestdoi/235815834)
Tariq Faquih. tofaquih/imputation_of_untargeted_metabolites v1.0. April 2020. doi:10.5281/zenodo.3778920.

# Introduction
A script designed to impute missing values in Metabolon metabolomics datasets. Two imputation methods can be used: MICE and KNN
This script is designed to impute missing values in Metabolon HD4 datasets using KNN and MICE imputations. It is technically possible to use it with any metabolomic dataset.

# Workflow: 
- The user provides two lists of metabolites: *group1* to be imputed using MICE-pmm or kNN-obs-sel; *group2* to be impute with zero.
- The user must also provide the variables t be used in the analysis after the imputation including the outcome.
- A correlation matrix is created for all the *group1* metabolites.
- The *group1* metabolites are split to complete cases metabolites(ccm) and incomplete cases metabolites(icm).
- For each icm 10 ccm with the highest absolute R correlation are selected. These will be used to impute the missing values in icm.
  - if the icm has more than 90% missing values OR if the number of non-missing values in less than the number of predictor variables + 20. 
  - This is done to because of two reasons:
  - to eliminate possibly mis-annotated metabolites or unannotated metabolites that are xenobiotic in nature, 
  - To ensure the availability of enough cases to perform the imputation. 
  - This issue prevents the MICE package from performing the imputation all together.
  - The invalid cases will be imputed to zero.
- The imputed results are returned with 3 objects; The imputed data, the summary of the imputation, the mean R of the ccm used for each icm.

# How to use:
1. Import the script <code>source(UnMetImp.R)</code>
2. Create a vector with the names of the *group1* (endogenous and/or unannotated metabolites) and *group2* (xenobiotics) metabolites.
3. Use the <code>UnMetImp</code> function.
  - **Usage**
    - <code>UnMetImp(DataFrame , 
                    imp_type = 'mice' , 
                    number_m = 5 , 
                    group1 , 
                    group2 = NULL , 
                    outcome=NULL,
                    covars=NULL, 
                    fileoutname = NULL , 
                    use_covars = FALSE , 
                    logScale = TRUE )</code>
  - **Arguments**
    - **DataFrame**: The full dataframe to be used with all the metabolites, covariables and the outcome. Must numeric. Must be a dataframe.
    - **imp_type**: String. Type of imputation to be used: <code>mice</code> or <code>knn</code>. Default is mice.
    - **number_m**: Numeric. For __imp_type == "mice"__ only. Number of imputations to be used. Default = 5.    
    - **group1**: Vector. Required. Vector with the names of metabolite columns. Will be imputed using the provided __imp_type__.
    - **group2**: Vector. Optional. Vector with the names of metabolite columns. Will be imputed to zero.
    - **outcome**: String. Required. The outcome variable to be used in the future analysis.
    - **covars**: Vector. Recommended. variables used in the future analysis. Will be returned with the imputed data.
    - **fileoutname**: String value. Optional. Saves the imputed output to a file.
    - **use_covars**: Logical. Optional. Whether the __covars__ will be used to impute the missing values in the metabolites. Default = FALSE.
    - **logScale**: Logical. Optional. Whether the values need to be log and scaled for the imputation. if TRUE, the values will be log and scaled then un-log and unscaled before returning the imputed output. If FALSE, script will assume you have log the values. Default = TRUE.
# Examples of running the scripts:
> - mydata: your data table in dataframe class format.
> - endoids: a user created vector containing the column names of the endogenous metabolites 
> - unknowns: a user created vector containing the column names of the unannotated metabolites 
> - xeono: a user created vector containing the column names of the xenobiotic metabolites 

## Running default knn:
  <code>source('Master_Script.r')
  
  knnimp <- UnMetImp(DataFrame = mydata,  group1 = c( endoids , unknowns )  , group2= xeon , covars = c('age', 'sex'), imp_type = 'knn' ,outcome = c('BMI') , logScale = TRUE )</code>

## Running default MICE:

 <code>source('Master_Script.r')
  
 miceimp <- UnMetImp(DataFrame = mydata,  group1 = c( endoids , unknowns )  , group2= xeon , covars = c('age', 'sex'), imp_type = 'mice', number_m = 5, outcome = c('BMI') ,  use_covars = FALSE , logScale = TRUE)</code>
 
> The mice package stores the output from the imputation step into the object class mids by default. This stores information about the imputation process used and the imputation datasets created. The with() and pool() need the object class mids as input to run the analysis on the datasets, calculate the estimate for each dataset then pools the estimates and standard errors using Rubin’s Rules.
To convert the object class mids to a “long” format: 

<code>require(‘mice’)
IMP <- miceimp$mids
       
Longformat =<-complete(IMP ,  action = 'long' , include = TRUE)</code>
  
To convert the Longformat back to mids class:

<code>IMP <- as.mids(Longformat)</code>

To run the analysis on the mids class and pool the estimates by Rubin's Rules:

<code>Model_Formula = as.formula('BMI~age+sex+...')
  
 Mysummary <- summary(pool( with(data = IMP,  expr = lm( formula = Model_Formula ) ))  ,conf.int = TRUE) </code>
  
