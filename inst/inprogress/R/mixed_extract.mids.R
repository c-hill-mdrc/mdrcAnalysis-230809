#' Generate mixed extract given a dataframe or mids object. Will fall back
#' to alternative analysis when mixed model fails to converge. Differing
#' behavior with mids and dataframe object / mixed and alternate analyses
#' currently handled with basic control-flow if state.ments. Could be
#' accomplished with S3 method-class if desirable. Priamry advantage would be
#' easier extendability moving forward.
#'
#' @param .mids Data to create model on. Must inherit class mids
#' @param .dependents Vector of dependent variables to create extracts on
#' @param .rhs Optional string to use as right hand side lme4::lmer formula.
#' Required for more complex models
#' @param .fixed_effect Vector of fixed effects if .rhs not specified
#' @param .random_slopes Vector of random slopes to insert if .rhs unspecified
#' @param .random_intercept String specifying variable to use for random
#' intercept If vector of length > 1 passed, assumed to be mlm of > 2 levels.
#' Grouping structure of levels taken from order of vector
#' @param .by_group String specifying variable performing comparisons across.
#' Must be specified in .rhs passing .rhs explicitly. First level in .by_group
#' will be used as the reference value
#' @param .alternate Alternate method if mixed model fails to converge
#' @param .file_path File path to write excel to
#' @param .file_name File name to write excel to
#' @return Dataframe displaying relevant mixed model statistics
mixed_extract.mids <- function(.mids,
                               .dependents,
                               .rhs,
                               .fixed_effects = character(),
                               .random_slopes = character(),
                               .random_intercept = character(),
                               .by_group = character(),
                               .alternate = c("cluster", "aggregate", "skip"),
                               .file_path,
                               .file_name) {
  # Selecting alternate method
  .alternate <- match.arg(.alternate)
  # Creating extractable dataframe
  .df <- mice::complete(.mids, "long", TRUE)
  # Ensuring that .by_group is a factor if defined
  if (is_scalar_character(.by_group)) {
    # Ensuring that .by_group is a factor
    if (!is.factor(.df[[.by_group]])) {
      .df <- .df %>%
        dplyr::mutate(!!ensym(.by_group) := factor(.data[[.by_group]]))
    }
    # Saving "first" value as reference value
    .by_reference <- levels(.df[[.by_group]])[[1]]
    .mids <- suppressWarnings(mice::as.mids(.df))
  } else {.by_reference <- character()}
  # Generates rhs of formula if not explicitly specified
  if (missing(.rhs)) {
    # If random intercept passed as a vector length > 1, coerce to multi-level
    if (length(.random_intercept) > 1) {
      paste0(.random_intercept, sep = "/")
    }
    # Generating RHS by pasting components
    .rhs <- paste(.fixed_effects,
                  .by_group,
                  paste0("(1 + ", .random_slopes, ")"),
                  paste0("(1|", .random_intercept, ")"),
                  sep = " + ", collapse = " + ")
  }
  # Creating a model without covariates to test if mixed model will converge
  .test_rhs <- .rhs %>%
    # Extracting components inside of parantheses
    stringr::str_extract_all("\\(.*?\\)") %>%
    unlist() %>%
    paste(collapse = " + ")
  # Create a vector of method to call
  model_calls <- purrr::map_chr(.dependents, function(.dependent) {
    # If dependent variable is categorical (a factor), use aggregate_dummied
    if (is.factor(dplyr::pull(.df, {{ .dependent }}))) {
      return("mixed_extract_aggregate_dummied")
    }
    # Generating test formula, calling formula to capture environment
    test_formula <- stats::formula(paste(.dependent, "~", .test_rhs))
    # Checking for singularity in any model generation
    test_result <- .mids$m %>%
      seq_len() %>%
      purrr::map_lgl(
        ~lme4::isSingular(lme4::lmer(test_formula, mice::complete(.mids, .)))
      ) %>%
      any() %>%
      # Returning alternate if is singular
      dplyr::if_else(.alternate, "mixed") %>%
      suppressMessages()
    # If aggregate, set to numeric
    test_result <- dplyr::if_else(test_result == "aggregate",
                                  "aggregate_numeric", test_result)
    # Making the string callable
    test_result <- paste0("mixed_extract_", test_result)
    test_result <- "mixed_extract_mixed"
    return(test_result)
  })
  extract_dataframe <- purrr::map2_dfr(
    model_calls, .dependents,
    ~do.call(.x, list(.mids, .y, .rhs, .by_group, .by_reference))
  ) %>%
    # Removing blank reference columns
    dplyr::select(-tidyselect::matches(
      paste0(c("P-Value", "Effect Size"), ".*", .by_reference)))
  # Ensuring ordering if cluster_var introduced late
  extract_dataframe <- extract_dataframe %>%
    dplyr::relocate(dplyr::all_of(dplyr::any_of(
      c("Dependent", "Method", "Formula", "Cluster_Var")))) %>%
    dplyr::relocate(dplyr::all_of(dplyr::any_of(
      c("sigma", "REMLcrit"))), .after = "BIC")
  # Writing to specified excel
  if (!missing(.file_path) & !missing(.file_name)) {
    create_excel(.file_path, .file_name, extract_dataframe)
  }

  return(extract_dataframe)
}

#' Specific method for mixed_extract, default mixed
#' @param .mids Input data
#' @param .dependent String specifying dependent variable
#' @param .rhs Right hand side of mixed model formula
#' @param .by_group String specifying .by_group
#' @param .by_reference String specifying reference value of \code{.by_group}
#' @return Tibble row summarizing mixed model results
mixed_extract_mixed.mids <- function(.mids,
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
  model_formula <- paste(.dependent, "~", .rhs)
  model <- .mids %>%
    mice:::with.mids(lmerTest::lmer(stats::formula(model_formula)))
  # Not sure what base statistics to get here. Just nimp and nobs for now
  output_row <- output_row %>%
    dplyr::bind_cols(broom.mixed::glance(mice::pool(model)))
  # Performing .by_group comparisons if set
  if (is_scalar_character(.by_group)) {
    # Calculate adjusted means, keeping as emmeans object
    model_means_raw <- try(emmeans::emmeans(model, .by_group))
    # Ensuring that singularity in emmmeans does not cause extract to halt
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

mixed_by_compare.mira <- function(.model,
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
      dplyr::across(c("Reference", .by_group),
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
  if (inherits(.model$analyses[[1]], "lmerModLmerTest")) {
    # Pooling all datasets together
    .df <- purrr::map_dfr(.model$analyses, ~.@frame)
  } else if (inherits(.model$analyses[[1]], "lm")) {
    .df <- purrr::map_dfr(.model$analyses, ~.$model)
  }
  .unadjusted_statistics <- .df %>%
    dplyr::group_by(!!ensym(.by_group)) %>%
    dplyr::summarise(
      `Unadjs. Mean` = mean(.data[[.dependent]]),
      `Unadjs. Std Dev` = sd(.data[[.dependent]]),
      .groups = "drop"
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
      "Effect Size" =  (estimate - .by_reference_mean) / .by_reference_sd
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
                     "Std. Error", "Adj. Mean", "P-Value", "P-Value Stars",
                      "Effect Size"),
      names_vary = "slowest"
    ) %>%
    dplyr::relocate(contains("Unadjs. Mean"),
                    contains("Unadjs. Std Dev"),
                    contains("Adj. Mean")) %>%
    dplyr::relocate(contains("Effect Size"), .after = tidyselect::last_col())
}
