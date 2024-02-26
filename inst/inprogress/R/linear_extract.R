#' Creates linear extracts
#' @description Modeled after the SURVEYREG procedure from SAS, computes linear
#' regression model with user provided specifications and generates an extract
#' in Excel if a file path is provided
#' @param .df Input data
#' @param .tmnt_vars Stirng specifying treatment variable
#' @param .dep_vars Character vector of dependent variables
#' @param .treatment_na Value of treatment variable to coerce to NA
#' @param .covars Optional, Character vector of covariates
#' @param .robust_se Boolean determining if should be done with robust SE
#' @param .cluster_se Optional, Character vector of variables to cluster by
#' @param .strata_var Optional, Not yet enabled
#' @param .by_x Optional, Not yet enabled
#' @param .where_x Optional, Used to subset linear model
#' @param .wt_var Optional, Weights for linear model
#' @param .total Optional, Not yet enabled
#' @param .rate Optional, Not yet enabled
#' @param .intercept Optional, Boolean determining if the implied intercept
#'  term should be included
#' @param .othermodeloptions Not yet enabled
#' @param .rpt_proc Not yet enabled
#' @param .rpt_out Not yet enabled
#' @param .file_path output File Path
#' @param .file_name output File Name
#' @param .replace Not yet enabled
#' @param .round Not yet enabled
#' @param .round_prob Not yet enabled
linear_extract <- function(.df,
                           .tmnt_vars,
                           .dep_vars,
                           .treatment_na = NA_character_,
                           # Opt
                           .covars = NULL,
                           .robust_se = T,
                           .cluster_se = NULL,
                           .strata_var = NULL,
                           .by_x = NULL,
                           .wt_var = NULL,
                           # Processing
                           .total = NULL,
                           .rate = NULL,
                           .intercept = T,
                           .othermodeloptions = NULL,
                           .rpt_proc = NULL,
                           .rpt_out = NULL,
                           .file_path = NULL,
                           .file_name = NULL,
                           .replace = NULL,
                           .round = NULL,
                           .round_prob = NULL) {
  # Checking required arguments
  if (is.missing(.df)) {
    stop("Input Dataset Missing")
  }

  if (!(is.character(.tmnt_vars) & length(.tmnt_vars) == 1)) {
    stop("Treatment Variable Missing")
  }

  if (!is.character(.dep_vars)) {
    stop("Dependent Variable(s) Missing")
  }

  # Creating a row for all treatment and dependent variable combinations
  # Generating list of all combinations
  variables_list <- tidyr::expand_grid(".tmnt_var" = .tmnt_vars,
                                       ".dep_var" = .dep_vars)

  output_ds <- variables_list %>%
    purrr::pmap_dfr(function(.tmnt_var, .dep_var) {
       # Filtering out where dependent variable is missing
      .df <- .df %>%
        dplyr::filter(!is.na({{.dep_var}}))

      # If value passed to .treatment_na, remove those values as well
      if (!is.na(.treatment_na)) {
        .df <- .df %>%
          dplyr::filter(!!rlang::ensym(.tmnt_var) != .treatment_na)
      }

      # Converting treatment var into a numeric for ease of model generation
      .df <- .df %>%
        dplyr::mutate(
          # Note, {{}} which is equivalent to !!enquo is not sufficient here.
          # Actually creates a univariate factor variable of the string .TmnVar
          # Unquoting an ensym specifies that a column variable is desired
          Treatment = factor(!!rlang::ensym(.tmnt_var)),
          Treatment = as.numeric(Treatment),
          Treatment = Treatment - 1
        )
      # Creating Formula for linear model
      ## String of covariates
      cv <- paste(.covars, collapse = " + ")
      ## RHS of base formula
      x <- paste("Treatment", cv, sep = " + ")
      ## Full Formula as string
      formula <- paste(.dep_var, x, sep = " ~ ")
      ## string to formula
      formula <- stats::as.formula(formula)
      # # Used to keep name to better read extract
      x1 <- paste(.tmnt_var, cv, sep = " + ")
      formula1 <- paste(.dep_var, x1, sep = " ~ ")

      # Generate Linear model, taking subset if .where_x is passed
      if (is.null(.where_x)) {
        lmfit <- stats::lm(formula, data = .df, weights = .wt_var)
      } else{
        lmfit <- stats::lm(formula, data = .df, weights = .wt_var,
                           # David implemented these optional parameters with
                           # eval parse. I know there are some risks to this so
                           # have some other implementation methods possible
                           subset = eval(parse(text = .where_x)))
      }

      # Explicit standard errors
      if (!is.null(.cluster_se)) {
        ## Take clustered standard error if parameter passed
        cl <- paste0("~", .cluster_se) %>%
          stats::as.formula()
        lm_se <- lmtest::coeftest(
          lmfit, vcov = sandwich::vcovCL,
          type = "HC1",
          df = dplyr::n_distinct(dplyr::pull(.df, {{.cluster_se}})) - 1,
          cluster = cl
        )
      } else if (.robust_se) {
        ## Take robust standard error if TRUE and not taking clustered
        lm_se <- lmtest::coeftest(lmfit, vcov = sandwich::vcovHC, type = "HC1")
      } else{
        ## Otherwise take unadjusted standard errors
        lm_se <- lmtest::coeftest(lmfit)
      }

      # Collecting results in tidy format
      estimates <- lm_se %>%
        broom::tidy()
      treatment_estimates <- estimates %>%
        dplyr::filter(term == "Treatment")
      intercept_estimates <- estimates %>%
        dplyr::filter(term == "(Intercept)")
      lmsummary <- stats::summary.lm(lmfit)

      # Generating Naive Model
      ## Naive Formula
      naive_formula <- as.formula(paste(.dep_var, "~", "Treatment"))
      ## Generating model, subsetting if parameter passed
      if (is.null(.where_x)) {
        lm_naive <- stats::lm(naive_formula, data = .df, weights = .wt_var)
      } else{
        lm_naive <- stats::lm(naive_formula, data = .df, weights = .wt_var,
                             subset = eval(parse(text = .where_x)))
      }

      ## Naive Standard Error, same logic as non-naive
      if (!is.null(.cluster_se)) {
        cl <- paste("~", .cluster_se) %>%
          stats::as.formula()
        lm_naive_se <- lmtest::coeftest(
          lm_naive, sandwich::vcovCL,
          type = "HC1",
          df = dplyr::n_distinct(dplyr::pull(.df, {{.cluster_se}})) - 1,
          cluster = cl
        )
      } else if (.robust_se) {
        lm_naive_se <- lmtest::coeftest(lm_naive,
                                      vcov = sandwich::vcovCL,
                                      type = "HC1")
      } else{
        lm_naive_se <- lmtest::coeftest(lm_naive)
      }

      # Collecting naive results in tidy format
      naive_estimates <- lm_naive_se %>%
        broom::tidy()
      naive_treatment_estimates <- naive_estimates %>%
        dplyr::filter(term == "Treatment")
      naive_intercept_estimates <- naive_estimates %>%
        dplyr::filter(term == "(Intercept)")
      naive_lmsummary <- stats::summary.lm(lm_naive)

      # Collecting F Statistic
      fstat <- lmsummary$fstatistic
      naive_fstat <- naive_lmsummary$fstatistic

      # Collecting number of observations
      obs <- .df %>%
        tidyr::drop_na({{.dep_var}}, Treatment, tidyselect::any_of(.covars)) %>%
        dplyr::count(Treatment)
      # Collecting standard deviations
      var_std_dev <- .df %>%
        tidyr::drop_na({{.dep_var}}, Treatment) %>%
        group_by(Treatment) %>%
        summarise(StdDev = sd(!!rlang::ensym(.dep_var)))
      # Generating output row here.

      output <- tibble::tibble(
        "Dependent" = .dep_var
      )
      output <- output %>% dplyr::mutate(
        "DepVarLabel" = NA_character_,
        "AdjMean_P" = magrittr::sum(treatment_estimates$estimate,
                                    intercept_estimates$estimate),
        "AdjMean_C" = intercept_estimates$estimate,
        "Impact_P_C" = treatment_estimates$estimate,
        "Stars_P_C"  = create_stars(treatment_estimates$p.value),
        "ProbT_P_C" = treatment_estimates$p.value,
        "StdErr_P_C" = treatment_estimates$std.error,
        "ProbF" = stats::pf(fstat[1], fstat[2], fstat[3], lower.tail = F),
        "Stars_ProbF" = create_stars(stats::pf(fstat[1], fstat[2], fstat[3],
                                        lower.tail = F)),
        "RSquare" = lmsummary$r.squared,
        "AdjRSquare" = lmsummary$adj.r.squared,
        "NObs_P" = obs[[which(obs$Treatment == 1), "n"]],
        "NObs_C" = obs[[which(obs$Treatment == 0), "n"]],
        "NObs" = sum(obs$n),
        "NClusters" = dplyr::if_else(!is.null(.cluster_se),
                                     dplyr::n_distinct(
                                       dplyr::pull(.df, {{.cluster_se}})),
                                     NA_integer_),
        "UnAdj_Mean_P" = magrittr::sum(naive_treatment_estimates$estimate,
                                       naive_intercept_estimates$estimate),
        "UnAdj_Mean_C" = naive_intercept_estimates$estimate,
        "Unadj_Impact_P_C" = naive_treatment_estimates$estimate,
        "UnadjProbT_P_C" =  naive_treatment_estimates$p.value,
        "Unadj_ProbF" = stats::pf(naive_fstat[1], naive_fstat[2],
                                  naive_fstat[3], lower.tail = F),
        "DepVarMean" = mean(dplyr::pull(
          tidyr::drop_na(.df, {{.dep_var}}, Treatment),
          {{.dep_var}})),
        "dep_varstdDev" = stats::sd(dplyr::pull(
          tidyr::drop_na(.df, {{.dep_var}}, Treatment),
          {{.dep_var}})),
        "StdDev_P" = var_std_dev[[which(var_std_dev$Treatment == 1), "StdDev"]],
        "StdDev_C" = var_std_dev[[which(var_std_dev$Treatment == 0), "StdDev"]],
        "EfSzAllSD_P_C" = Impact_P_C / StdDev_C,
        "TmntVar" = .tmnt_var,
        "where_x" = .where_x,
        "wt_var" = .wt_var,
        "Equation" = formula1,
        "Covariates" = paste(.covars, collapse = ", "),
        "ClusterVar" = .cluster_se,
        "strata_var" = .strata_var # Not enabled
      ) %>%
        dplyr::bind_cols(trail())
    })
  if (!missing(.file_path) & !missing(.file_name)) {
    create_excel(.file_path, .file_name, output_ds)
  }

  return(output_ds)
}
