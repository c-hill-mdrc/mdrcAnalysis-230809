#' Surveyreg Extract
#'
#' Performs regression analysis for sample survey data.
#'
#' @param .dataset Input data
#' @param .TmntVar Character vector of treatment variables
#' @param .DepVars Character vector of dependent variables
#' @param .TreatmentNA Value of treatment variable to coerce to NA
#' @param .CoVars Character vector of covariates
#' @param .RobustSE Boolean determining if should be done with robust SE
#' @param .ClusterSE Character vector of variables to cluster by if
#' clusterd SE desired
#' @param .StrataVar Not yet enabled
#' @param .ByX Variable to separate extract by
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
#' @returns Tibble containing results of linear regression analysis on sample
#' @export
surveyreg_extract <- function(.dataset,
                             .TmntVar,
                             .DepVars,
                             .TreatmentNA = NA_character_,
                             # Opt
                             .DepLabels = NA_character_,
                             .CoVars = NULL,
                             .RobustSE = TRUE,
                             .ClusterSE = NULL,
                             .StrataVar = NULL,
                             .ByX = NULL,
                             .WtVar = NULL,
                             .WhereX = NULL,
                             # Processing
                             .Total = NULL,
                             .Rate = NULL,
                             .Intercept = T,
                             .OtherModelOptions = NULL,
                             .RptProc = NULL,
                             .RptOut = NULL,
                             .output_path,
                             .output_file,
                             .Replace = NULL,
                             .Round = NULL,
                             .RoundProb = NULL){

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

  # Filtering for non-missing cluster var if set
  # .dataset is not going to have .ClusterSE

  if (!is.null(.ClusterSE)) {
    .dataset <- .dataset %>%
      tidyr::drop_na({{ .ClusterSE }})
  }
  # Creating a row for all .ByX value and dependent variable combinations
  # Generating list of all combinations

  # Pulling out the levels
  if (!missing(.ByX)) {
    .ByValues <- unique(.dataset[[.ByX]])
  } else{
    .ByValues <- NULL
  }

  variables_list <- tidyr::expand_grid(
                                       ".dataset" = list(.dataset),
                                       ".DepVar" = .DepVars,
                                       ".ByX" = .ByX,
                                       ".ByValue" = .ByValues)

  # So what is being piped in.
  # You have fixed data set here.
  # In general, it's harder to debug if structure this way.

  variables_list %>%
    purrr:pmap_dfr(
      nameofthefunction
    )

  OutputDS <- variables_list %>%
    purrr::pmap_dfr(function(
      .dataset,
      .DepVar,
      .ByValue){

      # Filtering out where dependent variable is missing
      # the lm model could take care of it
      # but want to determine what our sample size is
      # so there is so utility to these series of functions
      .dataset01 <- .dataset %>%
        dplyr::filter(!is.na({{.DepVar}}))

      # If value passed to .ByValues, filter for it here
      if (!missing(.ByValue)) {
        if (is.na(.ByValue)) {
          .dataset02 <- .dataset01 %>%
            dplyr::filter(is.na(!!sym(.ByX)))
        } else {
          .dataset02 <- .dataset01 %>%
            dplyr::filter(!!sym(.ByX) == {{.ByValue}})
        }
      }

      # If value passed to .TreatmentNA, remove those values as well
      # Moving 1 arm from the data set of the ones you do not wish to test
      # Maybe come back to it later for conceptually.
      if(!is.na(.TreatmentNA)){
        .dataset <- .dataset %>%
          dplyr::filter(!!rlang::ensym(.TmntVar) != .TreatmentNA)
      }

      # Converting treatment variable into a numeric for ease of model generation
      .dataset <- .dataset %>%
        dplyr::mutate(
          # Note, {{}} which is equivalent to !!enquo is not sufficient here.
          # Actually creates a univariate factor variable of the string .TmnVar
          # Unquoting an ensym specifies that a column variable is desired
          Treatment = factor(!!rlang::ensym(.TmntVar)), # 1 and 2 facotrization. You are converting to 0 and 1.
          Treatment = as.numeric(Treatment),
          Treatment = Treatment - 1
        )
      # Creating Formula for linear model
      ## String of covariates
      CV <- paste_na(.CoVars, collapse = " + ")
      ## RHS of base formula
      X <- paste_na("Treatment", CV, sep = " + ")
      ## Full Formula as string
      formula <- paste_na(.DepVar, X, sep = " ~ ")
      ## string to formula
      formula <- stats::as.formula(formula)

      # # Used to keep name to better read extract
      # ZH: This is to paste it into an excel extract
      X1 <- paste_na(.TmntVar, CV, sep = " + ")
      formula1 <- paste_na(.DepVar, X1, sep = " ~ ")

      # Generate Linear model, taking subset if .WhereX is passed
      # ZH: Needs to check here.
      # Why Weighted Least Square?
      # https://stats.stackexchange.com/questions/359205/what-are-the-correct-ways-of-weighting-linear-regression-model

      # ZH: Unsure of this eval(parse(text = .WhereX)) data subsetting.
      # Can look for another implementation method
      # It's grabbing the character variable then evaluating it which is essentially making it an object?
      # Without a column vector specifying where we are subsetting, not sure this option is working?
      # It has not been run this will have to be a separate test.
      # Usually I would subset the data directly data = subset(.dataset, .WhereX == SOMECONDITION)
      # Can dig deeper into the documentation.

      if(is.null(.WhereX)){
        lmfit <- stats::lm(formula, data = .dataset, weights = .WtVar)
      } else{
        lmfit <- stats::lm(formula, data = .dataset, weights = .WtVar,
                           # David implemented these optional parameters with
                           # eval parse. I know there are some risks to this so
                           # have some other implementation methods possible
                           # ZH: discard this subset because we already filter above.
                           # MvG: project staff can filter the data on their own.
                           subset = eval(parse(text = .WhereX)))
      }

      # Explicit standard errors
      # ZH Note: Reading on Clustered SE (the intuitive motivation is good)
      # https://en.wikipedia.org/wiki/Clustered_standard_errors

      # ZH Note: Robust SE: How do we know whether the results are going to be heteroskedastic or not?
      # https://economictheoryblog.com/2016/08/07/robust-standard-errors/#:~:text=%E2%80%9CRobust%E2%80%9D%20standard%20errors%20is%20a,linear%20unbiased%20estimator%20(BLUE).

      if(!is.null(.ClusterSE)){
        ## Take clustered standard error if parameter passed
        CL <- paste0("~", .ClusterSE) %>%
          stats::as.formula()
        lmSE <- lmtest::coeftest(
          lmfit, vcov = sandwich::vcovCL,
          type = "HC1",
          df = dplyr::n_distinct(dplyr::pull(.dataset, {{.ClusterSE}})) - 1,
          cluster = CL
        )
      } else if(.RobustSE){
        ## Take robust standard error if TRUE and not taking clustered
        lmSE <- lmtest::coeftest(lmfit, vcov = sandwich::vcovHC, type = "HC1")
      } else{
        ## Otherwise take unadjusted standard errors
        lmSE <- lmtest::coeftest(lmfit)
      }

      # collecting results in tidy format
      Estimates <- lmSE %>%
        broom::tidy()
      Treatment.Estimates <- Estimates %>%
        dplyr::filter(term == "Treatment")
      Intercept.Estimates <- Estimates %>%
        dplyr::filter(term == "(Intercept)")
      lmsummary <- stats::summary.lm(lmfit)

      # Generating Naive Model
      ## Naive Formula
      ## Essentially no Control variables
      naive.formula <- as.formula(paste(.DepVar, "~", "Treatment"))
      ## Generating model, subsetting if parameter passed
      if(is.null(.WhereX)){
        LmNaive <- stats::lm(naive.formula, data = .dataset, weights = .WtVar)
      } else{
        LmNaive <- stats::lm(naive.formula, data = .dataset, weights = .WtVar,
                             subset = eval(parse(text = .WhereX)))
      }

      ## Naive Standard Error, same logic as non-naive
      if(!is.null(.ClusterSE)){
        CL <- paste("~", .ClusterSE) %>%
          stats::as.formula()
        LmNaiveSE <- lmtest::coeftest(
          LmNaive, sandwich::vcovCL,
          type = "HC1",
          df = dplyr::n_distinct(dplyr::pull(.dataset, {{.ClusterSE}})) - 1,
          cluster = CL
        )
      } else if(.RobustSE){
        LmNaiveSE <- lmtest::coeftest(LmNaive, vcov = sandwich::vcovCL, type = "HC1")
      } else{
        LmNaiveSE <- lmtest::coeftest(LmNaive)
      }

      # Collecting naive results in tidy format
      Naive.Estimates <- LmNaiveSE %>%
        broom::tidy()
      Naive.Treatment.Estimates <- Naive.Estimates %>%
        dplyr::filter(term == "Treatment")
      Naive.Intercept.Estimates <- Naive.Estimates %>%
        dplyr::filter(term == "(Intercept)")
      Naive.lmsummary <- stats::summary.lm(LmNaive)

      # Collecting F Statistic
      Fstat <- lmsummary$fstatistic
      Naive.Fstat <- Naive.lmsummary$fstatistic

      # Collecting number of observations
      Obs <- .dataset %>%
        tidyr::drop_na({{.DepVar}}, Treatment, tidyselect::any_of(.CoVars)) %>%
        dplyr::count(Treatment)
      # Collecting standard deviations
      VarStdDev <- .dataset %>%
        tidyr::drop_na({{.DepVar}}, Treatment) %>%
        dplyr::group_by(Treatment) %>%
        dplyr::summarise(StdDev = sd(!!rlang::ensym(.DepVar)))

      # Generating Output row here.
      # Generating as blank any columns that don't do anything yet
      # Potentially issue I see here is that this necessitates that the treatment
      # variable be binary  to be interpretable - BB

      # Doing this by initializing tibble and then using mutate so we can use columns
      # immediately after creation
      Output <- tibble::tibble(
        "Dependent" = .DepVar
      )

      Output <- Output %>% dplyr::mutate(
        "AdjMean_P" = Treatment.Estimates$estimate + Intercept.Estimates$estimate,
        "AdjMean_C" = Intercept.Estimates$estimate,
        "Impact_P_C" = Treatment.Estimates$estimate,
        "Stars_P_C"  = create_stars(Treatment.Estimates$p.value),
        "ProbT_P_C" = Treatment.Estimates$p.value,
        "StdErr_P_C" = Treatment.Estimates$std.error,
        "ProbF" = stats::pf(Fstat[1], Fstat[2], Fstat[3], lower.tail = F),
        "Stars_ProbF" = create_stars(stats::pf(Fstat[1], Fstat[2], Fstat[3], lower.tail = F)),
        "RSquare" = lmsummary$r.squared,
        "AdjRSquare" = lmsummary$adj.r.squared,
        "NObs_P" = Obs[[which(Obs$Treatment == 1), "n"]],
        "NObs_C" = Obs[[which(Obs$Treatment == 0), "n"]],
        "NObs" = sum(Obs$n),
        "NClusters" = ifelse(!is.null(.ClusterSE),
                                     dplyr::n_distinct(dplyr::pull(.dataset, {{.ClusterSE}})),
                                     NA_integer_),
        "UnAdj_Mean_P" = Naive.Treatment.Estimates$estimate + Naive.Intercept.Estimates$estimate,
        "UnAdj_Mean_C" = Naive.Intercept.Estimates$estimate,
        "Unadj_Impact_P_C" = Naive.Treatment.Estimates$estimate,
        "UnadjProbT_P_C" =  Naive.Treatment.Estimates$p.value,
        "Unadj_ProbF" = stats::pf(Naive.Fstat[1], Naive.Fstat[2], Naive.Fstat[3], lower.tail = F),
        "DepVarMean" = mean(dplyr::pull(tidyr::drop_na(.dataset, {{.DepVar}}, Treatment), {{.DepVar}})),
        "DepVarStdDev" = stats::sd(dplyr::pull(tidyr::drop_na(.dataset, {{.DepVar}}, Treatment), {{.DepVar}})),
        "StdDev_P" = VarStdDev[[which(VarStdDev$Treatment == 1), "StdDev"]],
        "StdDev_C" = VarStdDev[[which(VarStdDev$Treatment == 0), "StdDev"]],
        "EfSzAllSD_P_C" = Impact_P_C/StdDev_C,
        "TmntVar" = .TmntVar,
        "Wherex" = .WhereX,
        "Wtvar" = .WtVar,
        "Equation" = formula1,
        "Covariates" = paste(.CoVars, collapse = ", "),
        "ClusterVar" = .ClusterSE,
        "StrataVar" = .StrataVar # Not enabled
      ) %>%
        dplyr::bind_cols(trail())
      # Attaching by column if specified
      if (!missing(.ByValue)) {
        Output <- Output %>%
          dplyr::mutate(
            !!.ByX := {{ .ByValue }},
            UniqueID = paste(Dependent, .ByValue, sep = "_")
          )
      }
      return(Output)
    })
  if (!missing(.DepLabels)) {
    # If labels provided as a character vector
    if (is.character(.DepLabels)) {
      labeldb <- tibble::tibble(
        COLUMN_NAME = .DepVars,
        COLUMN_LABEL = .DepLabels)
    } else {
      labeldb <- .DepLabels
    }
    OutputDS <- OutputDS %>%
      left_join(labeldb, by = c("Dependent" = "COLUMN_NAME")) %>%
      rename(Dependent_Label = "COLUMN_LABEL") %>%
      relocate(Dependent, Dependent_Label)
  }

  if (!missing(.output_file)) {
    create_excel(rlang::maybe_missing(.output_path),
                 .output_file,
                 OutputDS)
  }
  return(OutputDS)
}
