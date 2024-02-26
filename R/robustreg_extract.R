#' Robust Regression Extract
#'
#' Performs robust regression analysis.
#'
#' @param .dataset Input dataframe containing all needed variables
#' @param .TmntVar Character vector of treatment variables
#' @param .DepVars Character vector of dependent variables
#' @param .TreatmentNA Value of treatment variable to coerce to NA
#' @param .CoVars Character vector of covariates
#' @param .asfactorcovs character vector of .covariates that should be treated as
#' factors/categorical variables
#' @param .RobustSE Boolean determining if should be done with robust SE
#' @param .ClusterSE Character vector of variables to cluster by if
#' clustered SE desired
#' @param .StrataVar Not yet enabled
#' @param .ByX Variable with values for analysis running by
#' @param .WtVar Weights for linear model
#' @param .WhereX Character vector to subset data by
#' @param .Total Not yet enabled
#' @param .Rate Not yet enabled
#' @param .Intercept Boolean determining if the implied intercept term should be included
#' @param .OtherModelOptions Not yet enabled
#' @param .RptProc Not yet enabled
#' @param .RptOut Not yet enabled
#' @param .output_path Excel file path
#' @param .output_file Excel file name
#' @param .Replace Not yet enabled
#' @param .Round Not yet enabled
#' @param .RoundProb Not yet enabled
#' @param .DepLabels Labels for dependent variables, not enabled
#'
#' @returns Tibble containing results of linear regression analysis on sample
#' @export
robustreg_extract <- function(.dataset,
                              .TmntVar,
                              .DepVars,
                              .TreatmentNA = NA_character_,
                              # Opt
                              .DepLabels = NA_character_,
                              .CoVars = NA_character_,
                              .asfactorcovs = NA_character_,
                              .RobustSE = TRUE,
                              .ClusterSE = NA_character_,
                              .StrataVar = NA_character_,
                              .ByX = NA_character_,
                              .WtVar = NA_character_,
                              .WhereX = NA_character_,
                              # Processing
                              .Total = NA_character_,
                              .Rate = NA_character_,
                              .Intercept = T,
                              .OtherModelOptions = NA_character_,
                              .RptProc = NA_character_,
                              .RptOut = NA_character_,
                              .output_path = NA_character_,
                              .output_file = NA_character_,
                              .Replace = NA_character_,
                              .Round = NA_character_,
                              .RoundProb = NA_character_){

  # Checking required arguments
  if(missing(.dataset)){
    stop("Input Dataset Missing")
  }

  if(missing(.TmntVar)){
    stop("Treatment Variable Missing")
  }

  if(missing(.DepVars)){
    stop("Dependent Variable Missing")
  }

  # Pulling out the levels
  if (!missing(.ByX)) {
    .ByValues <- unique(.dataset[[.ByX]])
  } else{
    .ByValues <- NULL
  }

  # Check for if the variable is a factor
  # Print out message for if they have been converted to a factor.
  # If not, we will convert to a factor.

  if(!is.na(.asfactorcovs)){

    # check which ones are factors
    are_factors <- .dataset %>%
      dplyr::select(all_of(.asfactorcovs)) %>%
      dplyr::summarise(across(everything(), is.factor))

    if(sum(are_factors) < length(.asfactorcovs)){

      # pulling out non factors
      nonFactorVars <- .asfactorcovs[which(are_factors == 0)]

      # print the non factors
      print(paste0("The following variables are not factorized. We have factorized them for you. ", nonFactorVars))

      # convert to factors
      .dataset <- .dataset %>%
        dplyr::mutate(across(all_of(nonFactorVars), factor))

    } # factorizing the non factor variables
  }

  # creating a list of data frames,
  # dependent variables and subgroup values to robust regression by
  variables_list <- tidyr::expand_grid(".dataset" = list(.dataset),
                                       ".DepVar" = .DepVars,
                                       ".ByX" = .ByX,
                                       ".ByValue" = .ByValues)

  # Adding non varying additional parameters to the expanded list
  # to be supplied to the function
  variables_list$.TmntVar <- .TmntVar
  variables_list$.TreatmentNA <- as.character(.TreatmentNA)
  variables_list$.CoVars <- list(.CoVars)
  variables_list$.asfactorcovs <- .asfactorcovs # parameter setting for predictors which might be factors
  variables_list$.ClusterSE <- .ClusterSE
  variables_list$.RobustSE <- .RobustSE
  variables_list$.WhereX <-.WhereX
  variables_list$.WtVar <- .WtVar
  variables_list$.StrataVar <- .StrataVar

  models_df <- tidyr::as_tibble(variables_list)

  # Running all the models in one step with pmap and list_rbind
  # with the use of explicit functions
  OutputDS <- models_df %>%
    purrr::pmap(
      filter_helper
    ) %>%
    purrr::list_rbind()

  # Putting the unique identifier upfront
  # and sorting by unique identifier
  OutputDS <-
    OutputDS %>%
      dplyr::relocate(UniqueId) %>%
      dplyr::arrange(UniqueId)

  # if the subgroup is not missing, change the variable name to SubGrpVal
  if(!is.na(.ByX)){

    OutputDS <-
      OutputDS %>%
        dplyr::rename(SubGrpVal = {{.ByX}})
  }

  # Outputting to an excel spreadsheet
  if (!is.na(.output_path) & !is.na(.output_file)) {

    create_excel(
      .output_path = .output_path,
      .output_file = .output_file,
      .x = OutputDS
    )
  }

  return(OutputDS)
}

filter_helper <- function(
    .dataset,
    .DepVar,
    .ByX,
    .ByValue,
    .TmntVar,
    .TreatmentNA,
    .CoVars,
    .asfactorcovs,
    .ClusterSE,
    .WhereX,
    .WtVar,
    .StrataVar,
    .RobustSE
  ){

  # Placeholder values that may not be used but called upon
  # for reporting later on.

  .StrataVar <- .StrataVar

  # filtering out where dependent variable is missing
  # by default, the lm model would take care of missing values
  # but since we may want to determine the sample size,
  # there is some utility to this set up
  .dataset01 <- .dataset %>%
    dplyr::filter(!is.na({{.DepVar}}))

  # if we are not missing levels of a subgroup
  if (!missing(.ByValue)) {

    # if the level is a missing level, filter the
    # data set to the missing level
    if (is.na(.ByValue)) {
      .dataset02 <- .dataset01 %>%
        dplyr::filter(is.na(!!sym(.ByX)))
    } else {

    # otherwise, filter the data set to the level
      .dataset02 <- .dataset01 %>%
        dplyr::filter(!!sym(.ByX) == {{.ByValue}})
    }
  } else {

      .dataset02 <- .dataset01

  }

  # If value passed to .TreatmentNA, remove those values as well
  # Removing 1 arm from the data set of the ones you do not wish to test
  # ZH: Look into the literature for the specification of arms and how
  # arms in general are specified
    if(!is.na(.TreatmentNA)){
      .dataset02 <- .dataset02 %>%
      dplyr::filter(!!rlang::ensym(.TmntVar) != .TreatmentNA)
    }

  # Converting treatment variable into a numeric for ease of model generation
  .dataset03 <- .dataset02 %>%
    dplyr::mutate(
      # Note, {{}} which is equivalent to !!enquo is not sufficient here.
      # Actually creates a univariate factor variable of the string .TmnVar
      # Unquoting an ensym specifies that a column variable is desired
      Treatment = factor(!!rlang::ensym(.TmntVar)), # 1 and 2 factorization. You are converting to 0 and 1.
      Treatment = as.numeric(Treatment),
      Treatment = Treatment - 1
    )

  # creating formula for linear model
  modFormula.obj <- create_formula(
    .TmntVar = .TmntVar,
    .CoVars = unlist(.CoVars),
    .DepVar = .DepVar
  )

  # pulling out the formula object
  modFormula <- modFormula.obj$modFormula

  # fitting the linear model on the filtered data set
  lmfit <- stats::lm(modFormula,
                     data = .dataset03)

  # Apply emmeans to linear model
  # Change reference grid to treat only categorical variables as factors
  if (length(.asfactorcovs) > 0) {

    covkeep <- c(.TmntVar, .asfactorcovs)

  } else {

    covkeep <- .TmntVar

  }

  # Input emmeans here
  lmfit_rg <- emmeans::ref_grid(lmfit, cov.keep = covkeep)

  ## NOTE: In addition to providing the specs parameter (as with a single call), the data frame must also be specified
  ##       Otherwise, the emmeans functions do not know which data to use for the emmeans
  emm <- emmeans::emmeans(lmfit_rg, data = .dataset03, specs = .TmntVar)

  # Apply tidy to emmeans result
  ## Drop unnecessary columns
  emmtidy <-
    broom::tidy(emm) %>%
    dplyr::select(all_of(.TmntVar), "estimate", "std.error")


  # Explicit standard errors
  # ZH Note: Reading on Clustered SE (the intuitive motivation is good)
  # https://en.wikipedia.org/wiki/Clustered_standard_errors

  # ZH Note: Robust SE: How do we know whether the results are going to be heteroskedastic or not?
  # https://economictheoryblog.com/2016/08/07/robust-standard-errors/#:~:text=%E2%80%9CRobust%E2%80%9D%20standard%20errors%20is%20a,linear%20unbiased%20estimator%20(BLUE).


  if(!is.na(.ClusterSE)){

    # Calculating the degree of freedom for coeftest
    .degFreedomInput <-
      .dataset03 %>%
      dplyr::pull({{.ClusterSE}})

    .degFreedom <- length(unique(.degFreedomInput)) - 1

    # vcovCL{sandwich} -: pulling the column by which we want to cluster by
    CL <- .dataset03[[.ClusterSE]]

    # estimating the clustered standard error
    lmSE <- lmtest::coeftest(
      x = lmfit, # model object needs to be passed in
      vcov = sandwich::vcovCL, # covariance matrix specification
      type = "HC1",
      df = .degFreedom,
      cluster = CL) # supply a column by which we want to cluster

  } else if(.RobustSE){
    ## Take robust standard error if TRUE and not taking clustered
    lmSE <- lmtest::coeftest(lmfit, vcov = sandwich::vcovHC, type = "HC1")
  } else{
    ## Otherwise take unadjusted standard errors
    lmSE <- lmtest::coeftest(lmfit)
  }

  ###############################
  # Estimating marginal means the old manual way but we are going to keep it for std errors adjustment
  ###############################

  # collecting results in tidy format
  Estimates <- lmSE %>%
    broom::tidy()

  # transfering estimatse into variables
  Treatment.Estimates <- Estimates %>%
    dplyr::filter(term == "treatment")
  Intercept.Estimates <- Estimates %>%
    dplyr::filter(term == "(Intercept)")

  lmsummary <- stats::summary.lm(lmfit)

  # generating Naive Model
  ## naive Formula
  ## essentially no control variables
  naive.formula.obj <- create_formula(
    .TmntVar = .TmntVar,
    .DepVar = .DepVar
  )

  # pulling out the formula object
  naive.formula <- naive.formula.obj$modFormula

  LmNaive <- stats::lm(naive.formula, data = .dataset03)

  ## Naive Standard Error, same logic as non-naive
  if(!is.na(.ClusterSE)){

    # Calculating the degree of freedom for coeftest
    .degFreedomInput <-
      .dataset03 %>%
      dplyr::pull({{.ClusterSE}})

    .degFreedom <- length(unique(.degFreedomInput)) - 1

    # vcovCL{sandwich} -: pulling the column by which we want to cluster by
    CL <- .dataset03[[.ClusterSE]]

    # estimating the clustered standard error
    LmNaiveSE <-
        lmtest::coeftest(
          LmNaive,
          sandwich::vcovCL,
          type = "HC1",
          df = .degFreedom,
          cluster = CL)

  } else if(.RobustSE){

    ## Take robust standard error if TRUE and not taking clustered
    LmNaiveSE <- lmtest::coeftest(LmNaive,
                                  vcov = sandwich::vcovCL,
                                  type = "HC1")
  } else {

      ## Otherwise take unadjusted standard errors
      LmNaiveSE <- lmtest::coeftest(LmNaive)
  }

  # Collecting naive results in tidy format
  Naive.Estimates <- LmNaiveSE %>%
    broom::tidy()
  Naive.Treatment.Estimates <- Naive.Estimates %>%
    dplyr::filter(term == "treatment")
  Naive.Intercept.Estimates <- Naive.Estimates %>%
    dplyr::filter(term == "(Intercept)")
  Naive.lmsummary <- stats::summary.lm(LmNaive)

  # Collecting F Statistic
  Fstat <- lmsummary$fstatistic
  Naive.Fstat <- Naive.lmsummary$fstatistic

  # Collecting number of observations
  if(all(is.na(unlist(.CoVars)))){

    Obs <- .dataset03 %>%
      tidyr::drop_na({{.DepVar}}, Treatment) %>%
      dplyr::count(Treatment)

  } else {

    Obs <- .dataset03 %>%
      tidyr::drop_na({{.DepVar}}, Treatment, tidyselect::any_of(.CoVars)) %>%
      dplyr::count(Treatment)

  }

  # Collecting standard deviations
  VarStdDev <- .dataset03 %>%
    tidyr::drop_na({{.DepVar}}, Treatment) %>%
    dplyr::group_by(Treatment) %>%
    dplyr::summarise(StdDev = sd(!!rlang::ensym(.DepVar)))

  # Putting everything in a tibble/table format
  # Generating Output row here.
  # Generating as blank any columns that don't do anything yet
  # Potentially issue I see here is that this necessitates that the treatment
  # variable be binary  to be interpretable - BB

  # Doing this by initializing tibble and then using mutate so we can use columns
  # immediately after creation
  Output <- tibble::tibble(
    "Dependent" = .DepVar
  )

  Output <- Output %>%
    dplyr::mutate(
      "AdjMean_1" = emmtidy$estimate[emmtidy$treatment == 1],
      "AdjMean_0" = emmtidy$estimate[emmtidy$treatment == 0],
      "Impact_1_0" = Treatment.Estimates$estimate,
      "Stars_1_0"  = create_stars(Treatment.Estimates$p.value),
      "ProbT_1_0" = Treatment.Estimates$p.value,
      "StdErr_1_0" = Treatment.Estimates$std.error,
      "ProbF" = stats::pf(Fstat[1], Fstat[2], Fstat[3], lower.tail = F),
      "Stars_ProbF" = create_stars(stats::pf(Fstat[1], Fstat[2], Fstat[3], lower.tail = F)),
      "RSquare" = lmsummary$r.squared,
      "AdjRSquare" = lmsummary$adj.r.squared,
      "NObs_1" = Obs[[which(Obs$Treatment == 1), "n"]],
      "NObs_0" = Obs[[which(Obs$Treatment == 0), "n"]],
      "NObs" = sum(Obs$n),
      "NClusters" = ifelse(!is.na(.ClusterSE),
                           dplyr::n_distinct(dplyr::pull(.dataset03, {{.ClusterSE}})),
                           NA_integer_),
      "UnAdj_Mean_1" = Naive.Treatment.Estimates$estimate + Naive.Intercept.Estimates$estimate,
      "UnAdj_Mean_0" = Naive.Intercept.Estimates$estimate,
      "Unadj_Impact_1_0" = Naive.Treatment.Estimates$estimate,
      "UnadjProbT_1_0" =  Naive.Treatment.Estimates$p.value,
      "Unadj_ProbF" = stats::pf(Naive.Fstat[1], Naive.Fstat[2], Naive.Fstat[3], lower.tail = F),
      "DepVarMean" = mean(dplyr::pull(tidyr::drop_na(.dataset03, {{.DepVar}}, Treatment), {{.DepVar}})),
      "DepVarStdDev" = stats::sd(dplyr::pull(tidyr::drop_na(.dataset03, {{.DepVar}}, Treatment), {{.DepVar}})),
      "StdDev_1" = VarStdDev[[which(VarStdDev$Treatment == 1), "StdDev"]],
      "StdDev_0" = VarStdDev[[which(VarStdDev$Treatment == 0), "StdDev"]],
      "EfSzAllSD_1_0" = Impact_1_0/StdDev_0,
      "EfSzOverall" = Impact_1_0/DepVarStdDev, # ZH: Added Effect Size Overall
      "EfSzPooled" = Impact_1_0/sqrt((((NObs_0 - 1) * StdDev_0) + ((NObs_1 - 1) * StdDev_1))/(NObs_0 + NObs_1 - 2)), # ZH: Added Effect Size Pooled
      "TmntVar" = .TmntVar,
      "Wherex" = base::ifelse(is.na(.WhereX), NA_character_, .WhereX),
      "Wtvar" = base::ifelse(is.na(.WtVar), NA_character_, .Wtvar),
      "Equation" = modFormula.obj$modFormulaString,
      "Covariates" = base::ifelse(all(is.na(unlist(unique(.CoVars)))), NA_character_,
                                  paste(unlist(unique(.CoVars)), collapse = ", ")),
      "ClusterVar" = .ClusterSE,
      "StrataVar" = .StrataVar # Not enabled
  ) %>%
    dplyr::bind_cols(trail())

  # Attaching by column if specified
  if (!missing(.ByValue)) {

    Output <-
      Output %>%
        dplyr::mutate(
          !!.ByX := {{ .ByValue }},
          UniqueId = paste(.ByValue,
                           Dependent, sep = "_")
        ) %>%
        dplyr::relocate(UniqueId, {{.ByX}}) %>%
        dplyr::arrange(UniqueId)
  } else {

    Output <-
      Output %>%
        dplyr::mutate(
          UniqueId = paste(Dependent)
        ) %>%
        dplyr::arrange(UniqueId)
  }

    Output <-
      Output %>%
        dplyr::relocate(UniqueId) %>%
        dplyr::arrange(UniqueId)

  return(Output)
}

create_formula <- function(
    .TmntVar = NA_character_,
    .CoVars = NA_character_,
    .DepVar = NA_character_
){

  # Creating Formula for linear model
  ## String of covariates
  CV <- paste_na(.CoVars, collapse = " + ")

  ## RHS of base formula
  X <- paste_na(.TmntVar, CV, sep = " + ")

  ## Full Formula as string
  modFormulaString <- paste_na(.DepVar, X, sep = " ~ ")

  ## string to formula
  modFormula <- stats::as.formula(modFormulaString)

return(list(modFormula = modFormula,
            modFormulaString = modFormulaString))

}



