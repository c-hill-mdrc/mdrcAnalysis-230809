#' Generate mixed extract given a dataframe object. Will fall back
#' to alternative analysis when mixed model fails to converge. Differing
#' behavior with mids and dataframe object / mixed and alternate analyses
#' currently handled with basic control-flow if state.ments. Could be
#' accomplished with S3 method-class if desirable. Priamry advantage would be
#' easier extendability moving forward.
#'
#' @param .df Data to create model on. Must inherit class data.frame
#' @param .dependents Vector of dependent variables to create extracts on
#' @param .rhs Optional string to use as right hand side lme4::lmer formula.
#' Required for more complex models
#' @param .fixed_effect Vector of fixed effects if .rhs not specified
#' @param .random_intercept String specifying variable to use for random
#' intercept If vector of length > 1 passed, assumed to be mlm of > 2 levels.
#' Grouping structure of levels taken from order of vector with first character
#' being lowest level.
#' @param .by_group String specifying variable performing comparisons across.
#' Must be specified in .rhs passing .rhs explicitly.
#' @param .by_reference Optional string specifying which value of .by_group is
#' used as reference value. Defaults to first value.
#' @param .alternate Alternate method if mixed model fails to converge
#' @param .file_path File path to write excel to
#' @param .file_name File name to write excel to
#' @return Dataframe displaying relevant mixed model statistics
#' @export
mixed_extract.data.frame <- function(
    .df,
    .dependents,
    .rhs,
    .fixed_effects = character(),
    .random_intercept = character(),
    .by_group = character(),
    .by_reference,
    .alternate = c("cluster", "aggregate", "skip"),
    .file_path,
    .file_name) {
  # Selecting alternate method
  .alternate <- match.arg(.alternate)

  # Ensuring that .by_group is a factor
  if (length(.by_group) == 1 & is.character(.by_group)) {
    if (!is.factor(.df[[.by_group]])) {
      .df <- .df %>%
        dplyr::mutate(!!.by_group := factor(.data[[.by_group]]))
    }
    # Saving reference value if doesn't exist
    if (!missing(.by_reference)) {
      .df[[.by_group]] <- forcats::fct_relevel(.df[[.by_group]], .by_reference)
    }
    .by_reference <- levels(.df[[.by_group]])[[1]]
  } else {
    # Saving filler value to pass
    .by_reference <- character()
  }
  # Generates rhs of formula if not explicitly specified
  if (missing(.rhs)) {
    # If random intercept passed as a vector length > 1, coerce to multi-level
    if (length(.random_intercept) > 1) {
      paste0(.random_intercept, sep = "/")
    }
    # Generating formula structure
    .random_intercept <- paste0("(1|",
                                paste(.random_intercept, collapse = "/"),
                                ")")
    # Ensuring that we don't repeat by_group
    if (length(.fixed_effects) > 0 & length(.by_group) > 0) {
      .fixed_effects <- .fixed_effects[.fixed_effects != .by_group]
    }
    # Generating RHS by pasting components
    .rhs <- stringr::str_c(
      .fixed_effects, .by_group,
      .random_intercept,
      sep = " + ", collapse = " + ")
  }
  # Creating a model without covariates to test if mixed model will converge
  .test_rhs <- .rhs %>%
    # Extracting components inside of parantheses
    stringr::str_extract_all("\\(.*?\\)") %>%
    unlist() %>%
    paste(collapse = " + ")
  # Create a vector of method to call

  # ZH: Why are we running a map function here?
  model_calls <- purrr::map_chr(.dependents, function(.dependent) {
    # If dependent variable is categorical (a factor), use aggregate_dummied
    if (is.factor(dplyr::pull(.df, {{ .dependent }}))) {
      return("mixed_extract_aggregate")
    }
    # Generating test model
    test_formula <- paste(.dependent, "~", .test_rhs)
    test_result <- suppressMessages(lmerTest::lmer(
      stats::formula(test_formula), .df)) %>%
      lme4::isSingular() %>%
      # Returning alternate if is singular
      dplyr::if_else(.alternate, "mixed")
    # If aggregate, set to numeric
    test_result <- dplyr::if_else(test_result == "aggregate",
                                  "aggregate", test_result)
    # Making the string callable
    test_result <- paste0("mixed_extract_", test_result)

    return(test_result)
  })

  extract_dataframe <- purrr::map2_dfr(
    model_calls, .dependents,
    ~do.call(.x, list(.df, .y, .rhs, .by_group, .by_reference))) %>%
    # Removing blank reference columns
    dplyr::select(-tidyselect::matches(
      paste0(c("P-Value", "Effect Size"), ".*", .by_reference)))
  # Ensuring ordering if cluster_var introduced late
  extract_dataframe <- extract_dataframe %>%
    dplyr::relocate(dplyr::all_of(dplyr::any_of(
      c("Dependent", "Method", "Formula", "Cluster_Var")))) %>%
    dplyr::relocate(dplyr::all_of(dplyr::any_of(
      c("sigma", "REMLcrit"))), .after = `BIC`)
  # Writing to specified excel
  if (!missing(.file_path) & !missing(.file_name)) {
    create_excel(.file_path, .file_name, extract_dataframe)
  }

  return(extract_dataframe)
}

# Writing internal function as method
mixed_extract_mixed <- function(.object, ...) {
  UseMethod("mixed_extract_mixed")
}

#' Specific method for mixed_extract, default mixed
#' @param .df Input data
#' @param .dependent String specifying dependent variable
#' @param .rhs Right hand side of mixed model formula
#' @param .by_group String specifying .by_group
#' @param .by_reference String specifying reference value of \code{.by_group}
#' @return Tibble row summarizing mixed model results
mixed_extract_mixed.data.frame <- function(.df,
                                           .dependent,
                                           .rhs,
                                           .by_group,
                                           .by_reference) {
  # Storing call information
  output_row <- tibble::tibble(
    "Dependent" = .dependent,
    "Method" = "Mixed",
    "Formula" = .rhs
  )

  # Generating proper model
  model_formula <- stats::formula(paste(.dependent, "~", .rhs))
  model <- lmerTest::lmer(model_formula, .df)

  # Collecting base statistics
  output_row <- output_row %>%
    dplyr::bind_cols(broom.mixed::glance(model))

  # Performing .by_group comparisons if set
  if (length(.by_group) == 1) {
    # Calculate adjusted means, keeping as emmeans object
    model_means_raw <- try(emmeans::emmeans(model, .by_group))
    # Ensuring that singularity in emmmeans does not cause extract to halt
    if (inherits(model_means_raw, "try-error")) {
      output_row <- output_row %>%
        dplyr::mutate(By_Group = .by_group, Reference_Value = .by_reference,
                      .before = "nobs") %>%
        dplyr::mutate(Method = "Mixed (singularity in full model)")
      return(output_row)
    }
    .by_statistics <- mixed_by_compare(model,
                                       model_means_raw,
                                       .by_group,
                                       .by_reference,
                                       .dependent)

    # Attaching .by_group statistics to output_row
    output_row <- output_row %>%
      dplyr::mutate(By_Group = .by_group, Reference_Value = .by_reference,
                    .before = "nobs") %>%
      dplyr::bind_cols(.by_statistics)
  }
  return(output_row)
}

#' Alterantive mlm extract generated using clustered lm instead of mlm
#' @return Tibble row summarizing clustered linear model results
mixed_extract_cluster <- function(.df,
                                  .dependent,
                                  .rhs,
                                  .by_group,
                                  .by_reference) {
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
  model <- stats::glm(model_formula, data = .df)

  output_row <- output_row %>%
    dplyr::bind_cols(broom::glance(model))

  # Creating by statistics if by_group specified
  if (is_scalar_character(.by_group)) {
    # Creating emmeans object with clustered standard erros
    model_means_raw <- try(emmeans::emmeans(
      model, .by_group,
      vcov = sandwich::vcovCL(
        model,
        cluster = stats::as.formula(paste("~", .cluster_var))
      )))
    # Ensuring that extract does not fail if emmeans is singular
    if (inherits(model_means_raw, "try-error")) {
      output_row <- output_row %>%
        dplyr::mutate(By_Group = .by_group, Reference_Value = .by_reference,
                      .before = "nobs")
      return(output_row)
    }
    .by_statistics <- mixed_by_compare(model,
                                       model_means_raw,
                                       .by_group,
                                       .by_reference,
                                       .dependent)
    # Attaching .by_group statistics to output_row
    output_row <- output_row %>%
      dplyr::mutate(By_Group = .by_group, Reference_Value = .by_reference,
                    .before = "nobs") %>%
      dplyr::bind_cols(.by_statistics)
  }
  return(output_row)
}

#' Alternative mlm extract generated using aggregated lm instead of mlm
#' @description Extract used as a multi-level model substitute for
#' categorical dependent variables. Does not follow work-flow process of other
#' functions because this process is so unique
#' @param .by_group Required, String specifying .by_group
#' @param .by_reference Required, String specifying reference value of .by_group
#' @return Tibble row summarizing aggregated linear model results
mixed_extract_aggregate <- function(.df,
                                    .dependent,
                                    .rhs,
                                    .by_group,
                                    .by_reference) {
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
  .df <- .df %>%
    tidyr::drop_na(tidyselect::any_of(
      c(.covars, .agg_var, .dependent, .by_group)
    ))
  # Creating dataframe aggregated to .agg_var
  # First making dependent variable into dummy variable
  .dummied_df <- fastDummies::dummy_cols(.df, .dependent) %>%
    janitor::clean_names()
  .agg_df <- .dummied_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(
      c(.dependent, .by_group, .covars)
    ))) %>%
    # Aggregating up by Mean, assumes that .covars are consistent across
    dplyr::summarise(dplyr::across(
      tidyselect::starts_with(.dependent), ~mean(., na.rm = T)),
      .groups = "drop")
  # Making .by_group into "dummy"
  # 1 if specific non-reference level, 0 if reference, NA otherwise
  # Extracting non_referrence group levels
  # Extracting group as vector
  .by_group_vector <- .agg_df %>%
    dplyr::pull({{ .by_group }})
  # Identifying non-reference levels
  .by_nonreference <- .by_group_vector %>%
    levels() %>%
    purrr::discard(~. == .by_reference) %>%
    set_names(., .)
  .dummied_by_group <- purrr::map_dfc(
    .by_nonreference, ~dplyr::case_when(
      .by_group_vector == . ~ 1L,
      .by_group_vector == .by_reference ~ 0L,
      TRUE ~ NA_integer_
    )
  )
  # Attaching dummied .by_group columns
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
  .by_reference_pvals <- purrr::map_dfc(.by_nonreference, function(.treatment) {
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
  .by_reference_pvals <- purrr::map_dfr(dependent_levels, ~return(.y),
                                        .by_reference_pvals)
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
    ~do.call(.x, list(.dummied_df, .y, .rhs, .by_group, .by_reference))
  ) %>%
    # Removing P-Values
    dplyr::select(-tidyselect::matches("P-Value"))
  output_rows <- dependent_level_dataframe %>%
    # Modifying Method for clarity
    dplyr::mutate(Method = paste("Aggregated Anova:", Method)) %>%
    # Attaching p-value columns, replicated
    dplyr::bind_cols(.by_reference_pvals) %>%
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

mixed_by_compare <- function(.object, ...) {
  UseMethod("mixed_by_compare")
}
#' Helper Function for extracts that use emmmeans::emmeans()
#' @param .model lm or lmerMod class object
#' @param .emmeans emmGrid class object
#' @param .by_group String specifying by_group used in emmeans
#' @param .by_reference String specifying reference value to use in .by_group
#' @param .dependent String specifying dependent variable
#' @keywords internal
mixed_by_compare.default <- function(.model,
                                     .emmeans,
                                     .by_group,
                                     .by_reference,
                                     .dependent) {
  # Storing emmeans as tidy object
  model_means <- .emmeans %>%
    broom.mixed::tidy()
  # Calculating p-value of differences between groups
  model_pairs <- .emmeans %>%
    emmeans:::pairs.emmGrid() %>%
    broom.mixed::tidy() %>%
    # Separating contrasts
    tidyr::separate(contrast,
                    into = c("Reference", .by_group),
                    sep = " - ",
                    remove = TRUE) %>%
    # Removing the prefixing pairs applies
    dplyr::mutate(
      dplyr::across(c("Reference", {{ .by_group }}),
                    ~stringr::str_remove(., stringr::fixed(.by_group))
      )) %>%
    # Taking only comparison with reference
    dplyr::filter(Reference == .by_reference) %>%
    # Selecting p-value and contrast level
    dplyr::select({{ .by_group }}, tidyselect::contains("p.value")) %>%
    dplyr::rename_with(~return("adj.p.value"),
                       tidyselect::contains("p.value")) %>%
    # Creating Stars Columns
    dplyr::mutate(stars = create_stars(adj.p.value))
  # Effect size: calculated as such:
  # # (effect_mean  - reference_mean)/reference_sd
  # # Calculating reference sd, dependent on type of model run
  # # Also extracting unadjusted here
  if (inherits(.model, "lmerModLmerTest")) {
    .df <- .model@frame
  } else if (inherits(.model, "lm")) {
    .df <- .model$model
  }
  .unadjusted_statistics <- .df %>%
    dplyr::group_by(!!ensym(.by_group)) %>%
    dplyr::summarise(
      `Unadjs. Mean` = mean(.data[[.dependent]]),
      `Unadjs. Std Dev` = sd(.data[[.dependent]]),
      `N` = sum(!is.na(.data[[.dependent]])),
      `N_miss` = sum(is.na(.data[[.dependent]])),
      .groups = "drop"
    )
  .pooled_statistics <- .df %>%
    dplyr::summarise(
      `Full Mean` = mean(.data[[.dependent]]),
      `Full Std Dev` = sd(.data[[.dependent]])
    )
  # Reference column sd
  .by_reference_sd <- .unadjusted_statistics %>%
    dplyr::filter(!!ensym(.by_group) == .by_reference) %>%
    dplyr::pull(`Unadjs. Std Dev`)

  # # Extracting reference mean
  .by_reference_mean <- model_means %>%
    dplyr::filter(!!ensym(.by_group) == .by_reference) %>%
    dplyr::pull(estimate)
  # Joining Relevant Statistics
  .by_statistics <- model_means %>%
    dplyr::select({{ .by_group }}, estimate, std.error) %>%
    dplyr::full_join(model_pairs, by = .by_group) %>%
    dplyr::full_join(.unadjusted_statistics, by = .by_group) %>%
    # Creating Effect Size
    dplyr::mutate(
      "Adj. Mean Diff" = estimate - .by_reference_mean,
      "Effect Size" =  `Adj. Mean Diff` / .by_reference_sd
    ) %>%
    # Renaming
    dplyr::rename("Adj. Mean" = estimate,
                  "Std. Error" = std.error,
                  "P-Value" = adj.p.value,
                  "P-Value Stars" = stars) %>%
    # Pivoting to a single row
    tidyr::pivot_wider(
      names_from = .by_group,
      names_glue = paste0("{.value}_{", .by_group, "}"),
      values_from = c("Unadjs. Mean", "Unadjs. Std Dev",
                      "Std. Error", "Adj. Mean", "Adj. Mean Diff",
                      "P-Value", "P-Value Stars",
                      "Effect Size", "N", "N_miss"),
      names_vary = "slowest"
    ) %>%
    dplyr::bind_cols(.pooled_statistics) %>%
    dplyr::relocate(contains("Unadjs. Mean"),
                    contains("Unadjs. Std Dev"),
                    contains("Adj. Mean"),
                    contains("Adj. Mean Diff")) %>%
    dplyr::relocate(contains("Effect Size"), .after = tidyselect::last_col())
}
