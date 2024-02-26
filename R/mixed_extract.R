#' Generate mixed extract given a dataframe object.
#'
#' @param .dataset Data to create model on. Must inherit class data.frame
#' @param .dependents Vector of dependent variables to create extracts on
#' @param .rhs Optional string to use as right hand side lme4::lmer formula.
#' Required for more complex models
#' @param .treatment character vector for .treatment/independent variable
#' @param .reference_level Optional string specifying which value of .treatment is
#' used as reference value. Defaults to first value.
#' @param .fixed_effects NOT IN USE Vector of fixed effects if .rhs not specified
#' @param .random_intercept NOT IN USE String specifying variable to use for random
#' intercept If vector of length > 1 passed, assumed to be mlm of > 2 levels.
#' Grouping structure of levels taken from order of vector with first character
#' being lowest level.
#' @param .alternate NOT IN USE Alternate method if mixed model fails to converge
#' @param .output_path File path to write excel to
#' @param .output_file File name to write excel to
#'
#' @return dataframe displaying relevant mixed model statistics
#'
#' @export mixed_extract
#'

mixed_extract <- function(
    .dataset,
    .dependents,
    .rhs = NA_character_,
    .treatment = NA_character_,
    .reference_level = NA_character_,
    .fixed_effects = NA_character_,
    .random_intercept = NA_character_,
    .alternate = c("cluster", "aggregate", "skip"),
    .output_path,
    .output_file) {

  ####################################################################################
  # Selecting alternate method - aggregate method is selected below.
  # ZH: How is only one alternative selected? Not getting how that is being passed on
  # ZH: What is match.arg doing? Need to discuss with Melvin
  ####################################################################################
  .alternate <- match.arg(.alternate)

  # Ensuring that .treatment is a factor
  if (length(.treatment) == 1 & !is.na(.treatment)) {

    if (!is.factor(.dataset[[.treatment]])) {

      .dataset <- .dataset %>%
        dplyr::mutate(!!.treatment := factor(.data[[.treatment]]))

    }

    # Saving reference value if doesn't exist
    if (!is.na(.reference_level)) {

      .dataset[[.treatment]] <- forcats::fct_relevel(.dataset[[.treatment]], .reference_level)

    }

    .reference_level <- levels(.dataset[[.treatment]])[[1]]

  } else {

    # Saving filler value to pass
    .reference_level <- character()

  }

  ############################################################
  # Code generating RHS if RHS is not provided
  ############################################################

  # Generates rhs of formula if not explicitly specified
  if (is.na(.rhs)) {

    # Print message notifying user that rhs contruction is not fully tested/implemented
    print("WARNING: Formula construction is not fully tested and implemented.")
    print("         For best results, please provide a properly specified formula.")

    # If random intercept passed as a vector length > 1, coerce to multi-level
    if (length(.random_intercept) > 1) {

        paste0(.random_intercept, sep = "/")

    }
    # Generating formula structure
    .random_intercept <- paste0("(1|",
                                paste(.random_intercept, collapse = "/"),
                                ")")
    # Ensuring that we don't repeat by_group
    if (length(.fixed_effects) > 0 & length(.treatment) > 0) {

      .fixed_effects <- .fixed_effects[.fixed_effects != .treatment]

    }
    # Generating RHS by pasting components
    .rhs <- stringr::str_c(
      .fixed_effects, .treatment,
      .random_intercept,
      sep = " + ", collapse = " + ")
  }

  #############################################################################
  # Creating a model without covariates to test if mixed model will converge
  # Creating the Null Model
  # But though the test is on convergence, the actual function used is for singularity check
  #############################################################################
  .test_rhs <- .rhs %>%
    # Extracting components inside of parentheses
    stringr::str_extract_all("\\(.*?\\)") %>%
    unlist() %>%
    paste(collapse = " + ")
  # Create a vector of method to call

  #############################################################################
  # if dependent variable is a categorical variable, go to aggregate
  # check aggregate below for further explanation of what it is doing
  #############################################################################

  # Not currently in use.

  #model_calls <-
  # purrr::map_chr(.dependents, function(.dependent) {
    # If dependent variable is categorical (a factor), use aggregate_dummied
    #if (is.factor(dplyr::pull(.dataset, {{ .dependent }}))) {
      #return("mixed_extract_aggregate")
  #  }

      #############################################################################
      # if dependent variable is a categorical variable, go to aggregate
      # check aggregate below for further explanation of what it is doing
      #############################################################################

      # test_formula <- paste(.dependent, "~", .test_rhs)

      # test_result <-
      #  #suppressMessages(lmerTest::lmer(
      #   stats::formula(test_formula), .dataset)) %>%
      #   lme4::isSingular() %>%
      #   # Returning alternate if is singular
      #   dplyr::if_else(.alternate, "mixed")
      #   # If aggregate, set to aggregate
      #   test_result <- dplyr::if_else(test_result == "aggregate",
      #                                 "aggregate", test_result)
      #
      # # Making the string callable method object!
      # test_result <- paste0("mixed_extract_", test_result)
      #
      #   return(test_result)
    # }) # End for anonymous function in purrr::map_chr (model_calls <- )

  #############################################################################
  # reordering the returned results
  # both model variance results and estimates have been returned
  #############################################################################

  extract_dataframe <-
    purrr::map2_dfr(
      "mixed_extract_mixed.data.frame", .dependents,
    ~do.call(.x, list(.dataset, .y, .rhs, .treatment, .reference_level))) %>%
    # Removing blank reference columns
    dplyr::select(-tidyselect::matches(
      paste0(c("P-Value", "Effect Size"), ".*", .reference_level)))

  # # Ensuring ordering if cluster_var introduced late
  # extract_dataframe <- extract_dataframe %>%
  #   dplyr::relocate(dplyr::all_of(dplyr::any_of(
  #     c("Dependent", "Method", "Formula", "Cluster_Var")))) %>%
  #   dplyr::relocate(dplyr::all_of(dplyr::any_of(
  #     c("sigma", "REMLcrit"))), .after = `BIC`)

  # Writing to specified excel
  if (!missing(.output_path) & !missing(.output_file)) {
    create_excel(.output_path, .output_file, extract_dataframe)
  }

  return(extract_dataframe)
}

## ZH: Not sure if the internal function here is serving any purpose!!

# Writing internal function as method
mixed_extract_mixed <- function(.object, ...) {
  UseMethod("mixed_extract_mixed")
}

#' Specific method for mixed_extract, default mixed
#' @param .dataset Input data
#' @param .dependent String specifying dependent variable
#' @param .rhs Right hand side of mixed model formula
#' @param .treatment String specifying .treatment
#' @param .reference_level String specifying reference value of \code{.treatment}
#' @return Tibble row summarizing mixed model results
mixed_extract_mixed.data.frame <- function(.dataset,
                                .dependent,
                                .rhs,
                                .treatment = NA_character_,
                                .reference_level = NA_character_) {


  # Setting browser for mixed_extract_mixed models

  # Storing call information
  output_row <- tibble::tibble(
    "Dependent" = .dependent,
    "Method" = "Mixed",
    "Formula" = .rhs
  )

  # Generating proper model
  model_formula <- stats::formula(paste(.dependent, "~", .rhs))
  model <- lmerTest::lmer(model_formula, .dataset)

  # Collecting base statistics
  output_row <- output_row %>%
    dplyr::bind_cols(broom.mixed::tidy(model)) %>%
    dplyr::bind_cols(broom.mixed::glance(model))

  ################################################################
  # Calculate adjusted means, keeping as emmeans object
  # when treatment does not exist for the null model, skip!
  ################################################################
  if(!is.na(.treatment)){

    model_means_raw <- try(emmeans::emmeans(model, .treatment,
                                            options = emmeans::get_emm_option("satterthwaite")))

  } else {

    model_means_raw <- ".treatment does not exist"

  }

  # Ensuring that singularity in emmeans does not cause extract to halt
  if (inherits(model_means_raw, "try-error")) {
    output_row <- output_row %>%
      dplyr::mutate(treatment = .treatment, Reference_Value = .reference_level,
                    .before = "nobs") %>%
      dplyr::mutate(Method = "Mixed (singularity in full model)")
    return(output_row)
  }

  ############################################################
  # if there is treatment and emmeans::emmeans exist
  ############################################################

  if(!is.na(.treatment)){

    .by_statistics <- mixed_by_compare(.dataset = .dataset,
                                       .model = model,
                                       .emmeans = model_means_raw,
                                       .treatment = .treatment,
                                       .reference_level = .reference_level,
                                       .dependent = .dependent)

    # Attaching .subgroup statistics to output_row
    output_row <- output_row %>%
      dplyr::filter(term == "treatment0") %>%
      dplyr::bind_cols(.by_statistics) %>%
      dplyr::select(-estimate, -std.error,
                    -statistic, -df,
                    -p.value, -nobs)

  } else {

    # return output row for null model
    ############################################################################
    # Rename column names to fall in standards with SAS and other column names
    ############################################################################
    output_row <- output_row %>%
      dplyr::mutate(Stars = create_stars(p.value))
  }

  return(output_row)
}

mixed_by_compare <- function(.object, ...) {
  UseMethod("mixed_by_compare")
}
#' Helper Function for extracts that use emmmeans::emmeans()
#' @param .model lm or lmerMod class object
#' @param .emmeans emmGrid class object
#' @param .subgroup String specifying by_group used in emmeans
#' @param .reference_level String specifying reference value to use in .subgroup
#' @param .dependent String specifying dependent variable
#' @keywords internal
mixed_by_compare.default <- function(
                               .dataset,
                               .model,
                               .emmeans,
                               .treatment,
                               .reference_level,
                               .dependent) {

  if(!is.character(.emmeans)){

    # Storing emmeans as tidy object
    model_means <- .emmeans %>%
      broom.mixed::tidy()
    # Calculating p-value of differences between groups
    model_pairs <- .emmeans %>%
      emmeans:::pairs.emmGrid() %>%
      broom.mixed::tidy() %>%
      # Creating Stars Columns
      dplyr::mutate(stars = create_stars(p.value))

    ## Mean Extract Insertion
    .unadjusted_statistics <-
      mean_extract(
        .dataset = .dataset,
        .subgroup = .treatment
      ) %>%
      dplyr::filter(Variable %in% .dependent) %>%
      dplyr::rename(
        Unadjs_Mean = Mean,
        Unadjs_StdDev = Stddev,
        N = N,
        N_miss = NMiss
      ) %>%
      dplyr::select(
        treatment,
        Unadjs_Mean,
        Unadjs_StdDev,
        N,
        N_miss
      )

    # Joining Relevant Statistics
    .by_statistics <- model_means %>%
      dplyr::select({{.treatment}}, estimate, std.error) %>%
      dplyr::full_join(.unadjusted_statistics, by = .treatment) %>%
      tidyr::pivot_wider(names_from = treatment, values_from = - treatment)


    # Add contrast calculation
    model_pairs <-
      model_pairs %>%
        dplyr::rename(
          !!paste0("Impact_", levels(.dataset[[.treatment]])[1], "_", levels(.dataset[[.treatment]])[2]) := estimate,
          !!paste0("StdErr_", levels(.dataset[[.treatment]])[1], "_", levels(.dataset[[.treatment]])[2]) := std.error,
          !!paste0("ProbT_", levels(.dataset[[.treatment]])[1], "_", levels(.dataset[[.treatment]])[2]) := p.value,
          !!paste0("Stars_", levels(.dataset[[.treatment]])[1], "_", levels(.dataset[[.treatment]])[2]) := stars
        ) %>%
        dplyr::select(
          !!paste0("Impact_", levels(.dataset[[.treatment]])[1], "_", levels(.dataset[[.treatment]])[2]),
          !!paste0("StdErr_", levels(.dataset[[.treatment]])[1], "_", levels(.dataset[[.treatment]])[2]),
          !!paste0("ProbT_", levels(.dataset[[.treatment]])[1], "_", levels(.dataset[[.treatment]])[2]),
          !!paste0("Stars_", levels(.dataset[[.treatment]])[1], "_", levels(.dataset[[.treatment]])[2])
        )

    # combine columns
    # relocate to rearrange the variables
    .by_statistics <-
      .by_statistics %>%
        dplyr::bind_cols(model_pairs)

    .by_statistics <- .by_statistics[, c(1,2,3,4,13,14,15,16,
                                         5,6,7,8,9,10,11,12)]


    return(.by_statistics)

  } # end of if for emmeans existing

}

# dplyr::rename(treatment = term,
#               AdjMean_1 = estimate_1,
#               AdjMean_0 = estimate_0,
#               Impact_1_0 = estimate,
#               StdErr_1_0 = std.error,
#               ProbF = p.value,
#               Stars_ProbF = Stars,
#               NObs = nobs,
#               UnAdj_Mean_1 = `Unadjs. Mean_1`,
#               UnAdj_Mean_0 = `Unadjs. Mean_0`,
#               UnAdj_StdDev_1 = `Unadjs. Std Dev_1`,
#               UnAdj_StdDev_0 = `Unadjs. Std Dev_0`,
#               Adj_Mean_1 = `Adj. Mean_1`,
#               Adj_Mean_0 = `Adj. Mean_0`,
#               Adj_Mean_Diff_1_0 = `Adj. Mean Diff_0`)



#' Alterantive mlm extract generated using clustered lm instead of mlm
#' @param .dataset data frame being passed by user
#' @param .dependent character value of desired outcome variable
#' @param .rhs model specification for predictors
#' @param .subgroup Required, String specifying .subgroup
#' @param .reference_level Required, String specifying reference value of .subgroup
#' @return Tibble row summarizing clustered linear model results
mixed_extract_cluster <- function(.dataset,
                                  .dependent,
                                  .rhs,
                                  .subgroup,
                                  .reference_level) {
  # Storing call information
  output_row <- tibble::tibble(
    "Dependent" = .dependent,
    "Method" = "Clustered Linear Model"
  )
  # Constructing rhs without random effects
  .linear_rhs <- .rhs %>%
    # Extracting separate terms
    stringr::str_split("\\s*\\+\\s*", simplify = TRUE) %>%
    # Removing random effects
    magrittr::extract(1, !stringr::str_detect(., "\\(|\\)")) %>%
    # Collapse remaining terms by " + "
    paste(collapse = " + ")
  output_row <- output_row %>%
    dplyr::mutate(Formula = .linear_rhs)

  # Assume we want to cluster "up" to highest grouping level, identifying var
  # Method depends on if 2+ levels specified as (1|level1:level2)
  # or (1|level1/level2)
  if (stringr::str_detect(
    # Checking for (1 +1 level1/level2)
    .rhs, "\\(1\\|[:alnum:][:graph:]*\\/[:alnum:][:graph:]*\\)")) {
    # Take string delimited by /
    .cluster_var <- stringr::str_extract(
      .rhs, "(?<=\\(1\\|).*(?=\\))") %>%
      # Separate by /
      stringr::str_split("\\/", simplify = TRUE) %>%
      # Take last word
      magrittr::extract(., length(.))
  } else{
    # Otherwise, assume it is of the form (1|level1:level2)
    .cluster_var <- .rhs %>%
      # Extract all random intercepts
      stringr::str_extract_all("(?<=\\(1\\|).*?(?=\\))",
                               simplify = TRUE) %>%
      # Split by :
      stringr::str_split(":", simplify = TRUE) %>%
      # Taking last word
      magrittr::extract(., length(.))
  }
  output_row <- output_row %>%
    dplyr::mutate(Cluster_Var = .cluster_var)

  # Creating base linear model
  # We could also use the surveyreg function to potentially do this
  # However, would need to ensure that they are interoperable
  model_formula <- stats::as.formula(paste(.dependent, "~", .linear_rhs))
  model <- stats::glm(model_formula, data = .dataset)

  output_row <- output_row %>%
    dplyr::bind_cols(broom::glance(model))

  # Creating by statistics if by_group specified
  if (is_scalar_character(.subgroup)) {
    # Creating emmeans object with clustered standard erros
    model_means_raw <- try(emmeans::emmeans(
      model, .subgroup,
      vcov = sandwich::vcovCL(
        model,
        cluster = stats::as.formula(paste("~", .cluster_var))
      )))
    # Ensuring that extract does not fail if emmeans is singular
    if (inherits(model_means_raw, "try-error")) {
      output_row <- output_row %>%
        dplyr::mutate(By_Group = .subgroup, Reference_Value = .reference_level,
                      .before = "nobs")
      return(output_row)
    }
    .by_statistics <- mixed_by_compare(model,
                                       model_means_raw,
                                       .subgroup,
                                       .reference_level,
                                       .dependent)
    # Attaching .subgroup statistics to output_row
    output_row <- output_row %>%
      dplyr::mutate(By_Group = .subgroup, Reference_Value = .reference_level,
                    .before = "nobs") %>%
      dplyr::bind_cols(.by_statistics)
  }
  return(output_row)
}

#' Alternative mlm extract generated using aggregated lm instead of mlm
#' @description Extract used as a multi-level model substitute for
#' categorical dependent variables. Does not follow work-flow process of other
#' functions because this process is so unique
#' @param .dataset data frame being passed by user
#' @param .dependent character value of desired outcome variable
#' @param .rhs model specification for predictors
#' @param .subgroup Required, String specifying .subgroup
#' @param .reference_level Required, String specifying reference value of .subgroup
#' @return Tibble row summarizing aggregated linear model results
mixed_extract_aggregate <- function(.dataset,
                                    .dependent,
                                    .rhs,
                                    .subgroup,
                                    .reference_level) {

  #browser()

  # Storing call information
  output_row <- tibble::tibble(
    "Dependent" = .dependent,
    "Method" = "Aggregated Linear Model"
  )
  # Assume we want to cluster "up" to highest grouping level, identifying var
  # Method depends on if 2+ levels specified as (1|level1:level2)
  # or (1|level1/level2)
  if (stringr::str_detect(
    # Checking for (1 +1 level1/level2)
    .rhs, "\\(1\\|[:alnum:][:graph:]*\\/[:alnum:][:graph:]*\\)")) {
    # Take string delimited by /
    .agg_var <- stringr::str_extract(
      .rhs, "(?<=\\(1\\|).*(?=\\))") %>%
      # Separate by /
      stringr::str_split("\\/", simplify = TRUE) %>%
      # Take last word
      magrittr::extract(., length(.))
  } else{
    # Otherwise, assume it is of the form (1|level1:level2)
    .agg_var <- .rhs %>%
      # Extract all random intercepts
      stringr::str_extract_all("(?<=\\(1\\|).*?(?=\\))",
                               simplify = TRUE) %>%
      # Split by :
      stringr::str_split(":", simplify = TRUE) %>%
      # Taking last word
      magrittr::extract(., length(.))
  }
  # Extracting covariates
  .covars <- .rhs %>%
    # Separating terms by +
    stringr::str_split("\\+", simplify = T) %>%
    stringr::str_remove_all("\\s") %>%
    # Taking only fixed effects
    stringr::str_extract("^[:alnum:].*") %>%
    purrr::discard(is.na)
  # Filtering where not NA
  .dataset <- .dataset %>%
    tidyr::drop_na(tidyselect::any_of(
      c(.covars, .agg_var, .dependent, .subgroup)
    ))
  # Creating dataframe aggregated to .agg_var
  # First making dependent variable into dummy variable
  .dummied_df <- fastDummies::dummy_cols(.dataset, .dependent) #%>%
  #janitor::clean_names()
  .agg_df <- .dummied_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(
      c(.dependent, .subgroup, .covars)
    ))) %>%
    # Aggregating up by Mean, assumes that .covars are consistent across
    dplyr::summarise(dplyr::across(
      tidyselect::starts_with(.dependent), ~mean(., na.rm = T)),
      .groups = "drop")
  # Making .subgroup into "dummy"
  # 1 if specific non-reference level, 0 if reference, NA otherwise
  # Extracting non_referrence group levels
  # Extracting group as vector
  .subgroup_vector <- .agg_df %>%
    dplyr::pull({{ .subgroup }})
  # Identifying non-reference levels
  .by_nonreference <- .subgroup_vector %>%
    levels() %>%
    purrr::discard(~. == .reference_level) %>%
    set_names(., .)
  .dummied_by_group <- purrr::map_dfc(
    .by_nonreference, ~dplyr::case_when(
      .subgroup_vector == . ~ 1L,
      .subgroup_vector == .reference_level ~ 0L,
      TRUE ~ NA_integer_
    )
  )
  # Attaching dummied .subgroup columns
  .agg_df <- dplyr::bind_cols(.agg_df, .dummied_by_group)

  # Identifying dependent levels as column names in aggregated model
  dependent_levels <- names(.agg_df) %>%
    purrr::keep(~str_detect(., stringr::fixed(.dependent))) %>%
    purrr::discard(~. == .dependent)

  # Formula for null model, contains only covars
  .null_rhs <- .covars %>%
    paste(collapse = " + ")
  # Formula for measure model, contains dummied measues
  .measure_rhs <- names(.agg_df) %>%
    purrr::keep(~stringr::str_detect(., stringr::fixed(.dependent))) %>%
    paste(collapse = " + ") %>%
    paste(.null_rhs, ., sep = " + ")
  .reference_level_pvals <- purrr::map_dfc(.by_nonreference, function(.treatment) {
    # Performing anova comparison between nuil and measure linear models
    c(.measure_rhs, .null_rhs) %>%
      paste(.treatment, "~", .) %>%
      # Generating linear model on null and measure formuilas
      purrr::map(~stats::lm(as.formula(.), data = .agg_df)) %>%
      # Performing anova
      purrr::reduce(stats::anova) %>%
      magrittr::use_series("Pr(>F)") %>%
      purrr::discard(is.na) %>%
      tibble::tibble("P" = .) %>%
      dplyr::mutate(Stars = create_stars(P)) %>%
      magrittr::set_names(paste0(c("P-Value_", "P-Value Stars_"), .treatment))
  })
  # Replicating p-vals until appropriate number of rows
  .reference_level_pvals <- purrr::map_dfr(dependent_levels, ~return(.y),
                                           .reference_level_pvals)
  # Adjusted mean and effect size calculated as if normal mixed model
  # Perform a call for each factor of dependent
  # # Identifying appropriate call for each dependent level

  dependent_level_calls <- dependent_levels %>%
    # Keeping only dummied form
    purrr::map_chr(function(.dependent_level) {
      # Generating null model
      test_formula <- stats::formula(paste(.dependent_level, "~", .rhs))
      test_result <- try(suppressMessages(
        lmerTest::lmer(test_formula, .dummied_df)
      ))
      if (inherits(test_result, "try-error")) {
        # Returning cluster if downdated VtV
        test_result <- "mixed_extract_cluster"
      } else{
        test_result <- test_result %>%
          lme4::isSingular() %>%
          # Returning cluster if singular
          dplyr::if_else("mixed_extract_cluster", "mixed_extract_mixed")
      }
      return(test_result)
    })
  # # Generating dependent levels statistics dataframe
  dependent_level_dataframe <- purrr::map2_dfr(
    dependent_level_calls, dependent_levels,
    ~do.call(.x, list(.dummied_df, .y, .rhs, .subgroup, .reference_level))
  ) %>%
    # Removing P-Values
    dplyr::select(-tidyselect::matches("P-Value"))
  output_rows <- dependent_level_dataframe %>%
    # Modifying Method for clarity
    dplyr::mutate(Method = paste("Aggregated Anova:", Method)) %>%
    # Attaching p-value columns, replicated
    dplyr::bind_cols(.reference_level_pvals) %>%
    dplyr::relocate(By_Group, Reference_Value, nobs,
                    .after = "Cluster_Var") %>%
    dplyr::relocate(tidyselect::contains("sigma"),
                    tidyselect::contains("REMLcrit"),
                    .after = "df.residual") %>%
    dplyr::relocate(tidyselect::matches("P-Value"),
                    tidyselect::matches("Effect Size"),
                    .after = dplyr::last_col()
    )
}
