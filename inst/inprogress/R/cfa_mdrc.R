# author: James McCarthy
# date: 5/11/22
# purpose: function to run confirmatory factor analysis

# factor_names: vector of factor names - should directly correlate with formulas
# factor_formulas: vector of formulas for each factor name
# data: df where factors come from
# amys_output:logical with default=FALSE. If TRUE, simplified output of cfa summary. 
#   comparative fit index (CFI) and
#   root mean squared error of applicaiton (RMSEA) are only output fit statistics. 
#   Standardized factor loadings for each item are output as well (std.all) 
# note: names in factor_names cannot be the same names as existing variables 
cfa_mdrc <- function(factor_names, factor_formulas, data, amys_output = F) {
  
  #issue warning if lenghts of factor_names and factor_formulas are not the same
  if (length(factor_names) != length(factor_formulas)) {
    print("Warning: Input vectors are not the same length")
  }
  
  #create factor model 
  #loop through input vectors to construct the model
  for (i in seq_along(factor_names)) {
    #grab the i'th factor's name and formula
    name <- factor_names[i]
    formula <- factor_formulas[i]
    
    #create overall lavaan formula pasting together name and formula 
    factor <- paste(name, formula, sep = " =~ ")
    
    if (i == 1) { #if first iteration, assign model to first factor equation
      model <- factor
    } else { #if more factors, paste new formulas to first model, separated by \n
      model <- paste(model, factor, sep = "\n")
    }
    
  }
  
  # fit the cfa model to the input data
  fit <- lavaan::cfa(model, data=data)
  
  if (amys_output == F) {
    #summarize all the output
    #this will be printed when the function is called
    return(summary(fit, fit.measures = T, standardized = T))
    
    
  } else {
    #assign the fit summary to an object
    fit_summary <- summary(fit, fit.measures = T, standardized = T)
    
    #grab the fit statistics
    fit_stats <- as_tibble(fit_summary$FIT) %>%
      rename("Val" = "value")
    #grab the fit statistics names
    fit_names <- as_tibble(names(fit_summary$FIT)) %>%
      rename("Statistic" = "value")
    
    #bind fits statistics to names, round, keep only cfi and rmsea stats
    fit_overall <- bind_cols(fit_names, fit_stats) %>%
      mutate(across(where(is.numeric), ~round(., 3))) %>%
      filter(Statistic %in% c("cfi", "rmsea"))
    
    #now interested in the itemized factor loadings
    #get tibble of these from summary
    factor_loadings <- as_tibble(fit_summary$PE) %>%
      select(lhs, rhs, std.all) %>%
      filter(lhs %in% factor_names & !(rhs %in% factor_names))
    
    output <- list(fit_overall, factor_loadings)
    
    return(output)
  }
  
}