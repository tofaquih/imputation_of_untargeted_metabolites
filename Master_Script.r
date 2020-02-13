#########################
#This script is designed to impute missing values in Metabolon HD4 datasets. However, it is technically possible to use it with any dataset.
#The script works as follows: 1- The variables in group are supplied by the user. These variables will be imputed using MICE or KNN
#2- Variables in group 2 will be imputed with zero values
#3- A pearon's correlation matrix is created using all the group 1 variables
#4a- group 1 variables will be split into completed cases (no missing values), incomplete cases, and invalid cases where the number of missing is too high. 
#4b- The cutoff for this is set at 90% missingness OR if the number of non missing values in less than the number of predictor variables + 20. This is done to because of two reasons: 1- to elimnate possibly mis-annotated metabolites or unannotated metabolites that are xenobiotic in nature, 2- to ensure the avaliability of 
#5- 

#Function that selects the 10 strongest correlated metabolites
get_aux_mets <- function(data_cor ,  Met , ccm) {
    #print(Met)
    cor_x = data_cor[which(colnames(data_cor) == Met) ,]
    cor_x <- abs(cor_x[names(cor_x) %in% ( ccm )])
    #print(cor_x)
    Preds = c(na.omit(sort(cor_x , decreasing = T)[1:10]))
    return(Preds)
}

#a function that scales the data (z-transform) and stores the values for reversing the scaling
getscaleval <- function(x, dat) {
    y <- scale(dat[x], center = T, scale = T) 
    return(data.frame('scale' = attr(y , 'scaled:scale') , 'center' = attr(y , 'scaled:center') ))
}
#function that unsclaes the values using the orginal mean
unscale <- function (x, d) {
    #d = dataframe, x = variable name 
    d[x] * scale_info[x,'scale'] + scale_info[x,'center']
}
#function for knn imputatiion, requires the metabolite with missing valuesm the main data, the correlation matrix, the names of metabolites with no missing values
#the 10 strongest metabolites are selected and used to make a subdata with X
#the knn function from VIM is used to impute
useknn2 <- function(X , dataframe , data_cor, ccm  ) {
    require(VIM)
    Preds = get_aux_mets(data_cor , X ,ccm  )
    mean_qual = mean(Preds , na.rm = T)
    QualSummary <<- rbind(QualSummary , data.frame('Metabolite' = X , 'Mean_Quality' = mean_qual))
    cluster = c(names(Preds) , X )
    minidf <- dataframe[cluster]

    output <- as.data.frame(kNN(data = minidf,
                               variable = X ,
                               dist_var = cluster[1:length(cluster)-1] ,
                               numFun = mean ,
                               k=10))[X]
    #print(head(output[X]))
    return(output[X])
    }
#function that runs the mice imputation, requires X, the main data, the correlation matrix, the metabolites with no missing values
#long to adjust the number of variables in the output of this function (TRUE only if the loop iteration is n>1)
#optional argumnets are the O for the outcome, co_vars for the covariants that will be needed for the analysis model(s),
#m is the number of imputations to be done, use_co_vars is true only if you want to use these variables as predcitors for the imputation
usemice2 <- function(X , dataframe , data_cor, O, co_vars, ccm, long , m , use_co_vars = FALSE ) {
    #step 1 creating a subset for X and the predictor variables to be used to calculate the imputation values
    
    Preds = get_aux_mets(data_cor , X ,ccm  )
    mean_qual = mean(Preds , na.rm = T)
    if (use_co_vars == TRUE){
        cluster = c(O, co_vars , names(Preds) , X )}
    else {cluster = c(O , names(Preds) , X )
         co_vars = NULL}
    minidf <- dataframe[cluster]
    
    #the mean of the correlations is added to QualSummary
    QualSummary <<- rbind(QualSummary , data.frame('Metabolite' = X , 'Mean_Quality' = mean_qual))

    #Where matrix specifies which missing values to impute, we state that only missing values in X are to be imputed (set to TRUE) 
    #other values are set to FALSE
    where_matrix <- make.where(minidf)
    Dont_Impute <- cluster[!cluster %in% c(X)]
    where_matrix[, Dont_Impute] = FALSE
    
    #The actual imputation step using the mice package
    IMP <- mice(minidf ,printFlag = TRUE , m= m , where = where_matrix )
    
    #Long is set to FALSE if this the first variable in the loop, in this case only we select the .imp (imputation number), .id (the original row number), the outcome, the covars (only if use_co_vars is TRUE), X1.
    #For the remaining iterations in the loop we do not need these columns again so we only extract Xn

    if (long) {output <- complete(IMP ,  action = 'long' , include = TRUE)[X]}
    
    else {output <- complete(IMP ,  action = 'long' , include = TRUE)[c('.imp','.id', O , co_vars , X )]
         output[X] <- exp(as.data.frame(do.call(cbind , lapply(X , FUN =  unscale, d = output))))
         #output[X] <- exp(output[X])
         }
    return(output)
}



# the main function to be used by the user
#required args are the DataFrame and imp_type (mice or knn), group1:  a vector of the variable NAMES to be imputed with mice or knn
#group2 are variables to imputed with zero (xenobiotics)
#optional args: number_m: number of multiple imputations for mice
#outcome variable to used with mice imputation
#covars: variables that will be used in the later analysis after the mice imputation
#use_covars: set to TRUE if the co variables should be used for the imputation (not required for mice, IMPORTANT: will NOT work if the co vars have any missing values)
#fileoutname: name of the output name for the knn or mice imputation as an csv sheet

Main <- function(DataFrame , imp_type = 'mice' , number_m = 5 , group1 , group2 = NULL , outcome=NULL,
                 covars=NULL, fileoutname = NULL , use_covars = FALSE) {
    
    require(mice)
    require(dplyr)
    #ptm <- proc.time()
    #Step 1 log: normal dist, outliers effect reduced
    assign(x = 'scale_info' ,value = do.call(rbind , lapply(group1 , FUN =getscaleval, dat = log( DataFrame[group1] ) )),
                        envir =.GlobalEnv)
    #QualSummary will be used to store the mean correlation of each metabolite with missing values
    assign(x = 'QualSummary' ,value = data.frame() ,
                        envir =.GlobalEnv)
    
    
    #Step2: calculate number of missing values in each metabolites, split them into complete ccm and incomplete icm
    l =lapply(DataFrame[group1] , function (x) {sum(is.na(x))})
    ccm <- names(l)[l == 0]
    #If the missing is too large in the group1 metabolites then imputation becomes difficult or impossible.
    #The metabolite is most likely a xenobiotic and truly missing. These will be imputed to zero
    miss90 = nrow(DataFrame) * 0.9
    maxpreds = (20 + length(covars))
    if (miss90 > maxpreds) {cutoff = miss90} else {cutoff = maxpreds}
    icm <- l[l > 0 & l < cutoff ]
    icm_names <- names(icm)
    invalids <- names(l)[l >= cutoff ]
    
    
    #icm_names <- icm_names[!(invalids %in% icm_names) ]

    #Step3: create a correlation matrix
    data_cor = cor(DataFrame[c(ccm, icm_names)] , use="p")
    
    #Step4: check if there is actually missing values
    if (length(icm) == 0) {return('There are no missing data in the group 1 metabolites')}

    #Step5: set group2 missing values to zero and store their names
    group2_summary = 0
    if (is.null(group2) == FALSE ) {
        group2_summary = colnames(DataFrame[group2][!complete.cases(t(DataFrame[group2]))])
        DataFrame[group2][is.na(DataFrame[group2])] <- 0
    }   
    #Step 6: set the "invalid" (>90% missing etc) variables to zero
    if (length(invalids) > 0 ) {
        DataFrame[invalids][is.na(DataFrame[invalids])] <- 0
    }
    
    #Step 7: scale and log the ccm and icm metabolites ONLY
    DataFrame[c(ccm , icm_names)] <- scale(log( DataFrame[c(ccm , icm_names)] ))
    #Summary of the imputation used for which variables 
    msummary <- list('Imputed' = icm_names, 'Bad Case' = invalids, 'Zero' = group2_summary )
    
    #knn imputation option
    if (imp_type == 'knn') {
        #This calls the useknn2 function above. output is un-scaled and exponentialized and returned as object. 
        #can also be saved to a csv file if fileoutname is provided
        output = do.call(cbind, lapply(icm_names, FUN = useknn2, dataframe = DataFrame , data_cor = data_cor, ccm = ccm))
        output <- exp(as.data.frame(do.call(cbind,lapply(colnames(output),FUN =  unscale, d = output)))) 
        #the imputed variables replace their orginals in the main dataframe 
        DataFrame[colnames(output)] <- output
        
        if (is.null(fileoutname)) {
            #write.csv(DataFrame , file = paste('Imputed_Data_KNN_',Sys.Date(),'.csv', sep = '') , row.names = FALSE)
        }
        else {write.csv(DataFrame , file = paste(fileoutname , '_Imputed_Data_KNN','.csv', sep = '') , row.names = FALSE)}
        return(list(mids = DataFrame , Msummary = msummary, QS = QualSummary))
    }
    
    #mice imputation section
    else if ( imp_type == 'mice' ){
        #Only the first variable with missing values is imputed here (see usemice2 fucntion for details)
        firstmids <- usemice2(X = icm_names[1] , dataframe = DataFrame ,
                              data_cor = data_cor,
                              ccm = ccm ,
                              O= outcome,
                              co_vars = covars,
                              long = FALSE ,
                              use_co_vars = use_covars,
                              m= number_m
                             )

        #if , for whatever reason, there is only one metabolite with missingness, the next step will not be run
        #Otherwise the remaining variables will be imputed then merged with the togther and with the first variable
        if (length(icm_names) > 1) {
            #print('step 2')
            allmids <- as.data.frame(do.call(cbind , lapply(
                                            icm_names[2:length(icm_names)] , FUN = usemice2 ,
                                            dataframe = DataFrame ,
                                            data_cor = data_cor,
                                            ccm = ccm ,
                                            O=outcome,
                                            co_vars = covars,
                                            long = TRUE ,
                                            use_co_vars = use_covars,
                                            m=number_m)))
            
            allmids <- cbind(firstmids , exp(as.data.frame(do.call(cbind,lapply(colnames(allmids),
                                                                                FUN =  unscale, d = allmids)))) )
            #else {allmids  <- cbind(firstmids , exp(allmids))}
            }
            #}
        #only used if there is one variable with missing values
        else{allmids <- firstmids}

        #output is un-scaled and exponentialized and returned as object
        DataFrame[ccm] <-  exp(as.data.frame(do.call(cbind , lapply(ccm ,FUN =  unscale, d = DataFrame[ccm]))) )
        
        #the other variables from ccm and group2 must be merged with the multiple imputation data generated 
        #duplicate the rows of other metabolies number_m times then bind them to the allmids data frame

        othersDF <- DataFrame[c(ccm, invalids, group2)]
        #allmids  <- cbind(allmids ,othersDF[rep(1:nrow(DataFrame) , number_m+1) , colnames(othersDF) ])
        allmids[colnames(othersDF)] = othersDF[rep(1:nrow(DataFrame) , number_m+1) , colnames(othersDF) ]
        
        #we also need to add the variables that will be used in the analysis model(s). if use_covars was set to TRUE, then
        #these variables would have been added in the usemice2 function. otherwise, it is done here
        if (use_covars == FALSE) {
            covarsDF <- DataFrame[covars]
            allmids[covars] = covarsDF[rep(1:nrow(DataFrame) , number_m+1) , colnames(covarsDF) ]
            }


        
        #can also be saved to a csv file if fileoutname is provided
        if (is.null(fileoutname)) {
            #write.csv(allmids , file = paste('Imputed_Data_MICE_',Sys.Date(),'.csv', sep = '') , row.names = FALSE) 
        }
        else {write.csv(allmids , file = paste(fileoutname , '_Imputed_Data_MICE','.csv', sep = '') , row.names = FALSE)}
        
        #convert allmids to a "mids" object, the object format required by the mice package to run the analysis
        allmids <- as.mids(allmids)

        return(list(mids = allmids , Msummary = msummary , QS = QualSummary))
        }

}