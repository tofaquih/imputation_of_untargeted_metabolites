get_metabo_values <-function (P, Mname , DF, Mode = 'mar', seed = 1) {
    set.seed(seed)
    totalnum_miss = P * nrow(DF)
    
    Q <- quantile(DF[[Mname]]  , probs = c(0.05 , .25,.75))

    
    if( Mode == 'mar') {
        
        miss_per_section <-  c(totalnum_miss * 0.40 ,  totalnum_miss * 0.50 , totalnum_miss * 0.05)
        section1miss = sample(which(DF[Mname] <= Q[[2]]) , size =  miss_per_section[1])
        section2miss = sample(which(DF[Mname] > Q[[2]] & DF[Mname] <= Q[[3]]) , size =  miss_per_section[2])
        section3miss = sample(which(DF[Mname] > Q[[3]]) , size =  miss_per_section[3])
        
        DF[ c(section1miss, section2miss , section3miss ) , Mname] <- NA
        }

    
    else if (Mode == 'mcar') {
        allmiss = sample(1:nrow(DF) , size = nrow(DF) * P)
        DF[ allmiss , Mname] <- NA
        }
    
    else if (Mode == 'mnar') {
        LODmiss = which(DF[Mname] <= Q[[1]])
        DF[ LODmiss , Mname] <- NA
    }

    return(DF[Mname])
    }


