#' T-Test Extract
#'
#' Outputs frequencies, percents, and results of a t-test for two independent
#' samples to an excel sheet and dataframe
#' DOES NOT compare a mean to a constant value, perform a paired comparison,
#' or extract confidence levels
#'
#' @param .df Required, Input dataframe
#' @param .classvar Required, Categorical independent variable. Must have
#' exactly 2 different values
#' @param .depvars Optional, names of numeric dependent variables. Defaults to
#' all numeric variables
#' @param .depvar_labels Optional, character vector of labels to use for
#' .depvars
#' @param .byx Optional, variable to to split and group results by
#' @param .alternative Alternative hypothesis, defaults to two sided test
#' @param .mu Value indicating "true" difference in mean
#' @param .output_path Optional, path to output file location
#' @param .output_file Optional, output file name
#'
#' @returns Tibble describing difference within dependent variables between
#' binary categorical variable
#'
#' @export
ttest_extract <- function(.df,
                          .classvar,
                          .depvars,
                          .depvar_labels,
                          .byx = NULL,
                          .alternative = c("two.sided",
                                           "less",
                                           "greater"),
                          .mu = 0,
                          .output_path = NA_character_,
                          .output_file = NA_character_) {
  # Match alternative hypothesis parameter to valid options
  alternative <- match.arg(.alternative)

  # Filtering out missing classvar rows
  .df <- .df %>%
    dplyr::filter(!is.na(!!ensym(.classvar)))
  # Taking total observations
  nobs <- nrow(.df)
  # Taking .classvar as vector, ignoring NA
  .class_values <- .df %>%
    dplyr::pull({{ .classvar }}) %>%
    unique()

  # Confirming that .classvar has exactly 2 levels
  if (length(.class_values) != 2) {
    stop(".classvar must have exactly 2 distinct values")
  }
  .class_val_1 <- .class_values[which(order(.class_values) == 1)]
  .class_val_2 <- .class_values[which(order(.class_values) == 2)]

  # IF .depvars is missing, use all numeric that are unspecified by others
  if (missing(.depvars)) {
    .depvars <- dplyr::select(.df, tidyselect::where(is.numeric),
                              -{{ .classvar }}, -{{ .byx }}) %>%
      names()
  } else{
    # Otherwise ensure that .depvars are numeric
    .df <- .df %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(.depvars), as.numeric))
  }

  # If .byx defined, take all possible values ignoring NA
  if (!is.null(.byx)) {
    .by_vals <- .df %>%
      dplyr::pull({{ .byx }}) %>%
      unique() %>%
      purrr::discard(~is.na(.))
  } else {
    .by_vals <- NULL
  }

  # Generating labels matching if possible
  if (!missing(.depvar_labels)) {
    # If passed as a vector, coerce into a matchable dataframe
    if (rlang::is_bare_character(.depvar_labels)) {
      labeldf <- tibble::tibble("COLUMN_NAME" = .depvars,
                                "COLUMN_LABEL" = .depvar_labels)
    } else{
      labeldf <- .depvar_labels
    }
  } else{
    labeldf <- NULL
  }
  # Generating all combinations of .depvars and .byx
  parameter_grid <- tidyr::expand_grid(".dep_var" = .depvars,
                                       ".by_val" = .by_vals)
  # Writing function to loop through parameters
  ttest_results <- purrr::pmap_dfr(parameter_grid, function(.dep_var, .by_val) {
    if (!missing(.by_val)) {
      .df <- .df %>%
        dplyr::filter(!!sym(.byx) == {{ .by_val }})
    }
    .class_vector <- .df %>%
      dplyr::pull({{ .classvar }})
    .dep_vector <- .df %>%
      dplyr::pull({{ .dep_var }})
    .dep_vector_1 <- .dep_vector[.class_vector == .class_val_1]
    .dep_vector_2 <- .dep_vector[.class_vector == .class_val_2]
    ttest_equal_results <- stats::t.test(.dep_vector_1, .dep_vector_2,
                  alternative = .alternative,
                  mu = .mu,
                  var_equal = TRUE)
    ttest_unequal_results <- stats::t.test(.dep_vector_1, .dep_vector_2,
                  alternative = .alternative,
                  mu = .mu,
                  var_equal = FALSE)
    fprob <- stats::var.test(.dep_vector_1, .dep_vector_2)$p.value
    # Creating new names for ttest results tibble
    # Generate univarates on subsets and total vector
    univariate_stats <- list(.dep_vector_1, .dep_vector_2, .dep_vector) %>%
      purrr::map(function(.dep_subset, .index) {
        univariates <- tibble::tibble(
          "N" = sum(!is.na(.dep_subset), na.rm = TRUE),
          "Mean" = mean(.dep_subset, na.rm = TRUE),
          "NMISS" = sum(is.na(.dep_subset)),
          "STDDEV" = stats::sd(.dep_subset, na.rm = TRUE),
          "MEDIAN" = stats::median(.dep_subset, na.rm = TRUE),
          "MIN" = min(.dep_subset, na.rm = TRUE),
          "MAX" = max(.dep_subset, na.rm = TRUE))
      } )
    # Constructing final dataframe in order, skipping features to keep for all
    intermediate_output <- tibble::tibble("Variable" = .dep_var)
    if (!missing(.by_val)) {
      intermediate_output <- intermediate_output %>%
        dplyr::mutate(!!.byx := {{ .by_val }},
                      UniqueID = paste0({{ .by_val }}, "_", Variable))
    }
    # # Attaching labels if possible
    if (!is.null(labeldf)) {
      intermediate_output <- intermediate_output %>%
        dplyr::left_join(labeldf, by = c("Variable" = "COLUMN_NAME")) %>%
        dplyr::rename("VarLabel" = "COLUMN_LABEL")
    } else {
      intermediate_output <- intermediate_output %>%
        dplyr::mutate("VarLabel" = NA_character_)
    }
    final_output <- intermediate_output %>%
      dplyr::mutate(
        # Attaching N
        !!paste0("N_", .class_val_1) := univariate_stats[[1]]$N,
        !!paste0("N_", .class_val_2) := univariate_stats[[2]]$N,
        !!paste0("N_", "Total") := univariate_stats[[3]]$N,
        # Attaching Mean
        !!paste0("Mean_", .class_val_1) := univariate_stats[[1]]$Mean,
        !!paste0("Mean_", .class_val_2) := univariate_stats[[2]]$Mean,
        !!paste0("Diff") := magrittr::subtract(univariate_stats[[2]]$Mean,
                                               univariate_stats[[1]]$Mean),
        !!paste0("Mean_", "Total") := univariate_stats[[3]]$Mean,
        "EqualTValue" = ttest_equal_results$statistic,
        "UnequalTValue" = ttest_unequal_results$statistic,
        "EqualProbt" = ttest_equal_results$p.value,
        "UnequalProbt" = ttest_unequal_results$p.value,
        "ProbF" = fprob,
        "TValue" = dplyr::if_else(fprob < 0.05, UnequalTValue, EqualTValue),
        "ProbT" = dplyr::if_else(fprob < 0.05, UnequalProbt, EqualProbt),
        "EqualStars" = create_stars(EqualProbt),
        "UnequalStars" = create_stars(UnequalProbt),
        "Stars" = create_stars(ProbT),
        "ClassVar" = .classvar,
        "ClassVarLabel" = NA_character_, # Could add support for variable label
        !!paste0("NMISS_", .class_val_1) := univariate_stats[[1]]$NMISS,
        !!paste0("NMISS_", .class_val_2) := univariate_stats[[2]]$NMISS,
        !!paste0("STDDEV_", .class_val_1) := univariate_stats[[1]]$STDDEV,
        !!paste0("STDDEV_", .class_val_2) := univariate_stats[[2]]$STDDEV,
        !!paste0("MEDIAN_", .class_val_1) := univariate_stats[[1]]$MEDIAN,
        !!paste0("MEDIAN_", .class_val_2) := univariate_stats[[2]]$MEDIAN,
        !!paste0("MIN_", .class_val_1) := univariate_stats[[1]]$MIN,
        !!paste0("MIN_", .class_val_2) := univariate_stats[[2]]$MIN,
        !!paste0("MAX_", .class_val_1) := univariate_stats[[1]]$MAX,
        !!paste0("MAX_", .class_val_2) := univariate_stats[[2]]$MAX,
        "NMISS_Total" = univariate_stats[[3]]$NMISS,
        "NObs" = nobs,
        "STDDEV" = univariate_stats[[3]]$STDDEV,
        "MEDIAN" = univariate_stats[[3]]$MEDIAN,
        "MIN" = univariate_stats[[3]]$MIN,
        "MAX" = univariate_stats[[3]]$MAX
      )
  })

  # Outputting to an excel spreadsheet
  if (!is.na(.output_file) & !is.na(.output_path)){
    create_excel(.output_path = .output_path,
                 .output_file = .output_file,
                 .x = ttest_results)
  }
  return(ttest_results)
}
