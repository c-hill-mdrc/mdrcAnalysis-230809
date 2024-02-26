#' Correlation Extract
#'
#'
#' @description
#' Generates correlation extract dataframe using \link[psych]{corr.test}
#' which is itself a wrapper around \link[stats]{cor} . This function also
#' optionally attaches the p-value of a two-tailed probability of t and a Cronbach alpha
#' column if requested.
#'
#' @param .dataset Input tibble
#' @param .vars Optional, character vector of variables to correlate. Defaults to all
#' numeric variables not specified in other args (e.g. \code{.with_vars}, \code{.by_x})
#' @param .with_vars Optional, character vector of variables to correlate
#' with .vars
#' @param .by_x Optional, string of variable name to separate correlations for
#' @param .impute Optional, logical determining if missing values will be
#' imputed by mean replacement
#' @param .method Method to use for correlations, defaults to "pearson"
#' @param .use Determines if correlations will be calculated with "pairwise"
#' deletions or with "complete" data
#' @param .round Digits to round to, uses "exact" if NULL
#' @param .plot Boolean determining if correlation plot(s) included in Excel
#' @param .output_path Output path for excel file
#' @param .output_file Output name for excel file
#' @param .stats Character vector of other univariate methods to include in
#' output. Defaults to "mean"
#' @param .prob Boolean determining if t-test probability will be displayed. Calculated
#' using an alpha of 0.05.
#' @param .nobs Boolean determining if number of observations will be displayed
#' @param .alpha Boolean determining if Cronbach alphas will be calculated
#' @param .warnings Boolean determining if Warnings will be printed on run
#'
#' @returns Dataframe of correlations and other requested output
#'
#' @export
corr_extract <- function(.dataset,
                         .vars = NULL,
                         .with_vars = NULL,
                         .by_x = NULL,
                         .impute = FALSE,
                         .method = c("pearson",
                                     "spearman",
                                     "kendall"),
                         .use = c("pairwise.complete.obs",
                                  "all.obs",
                                  "complete.obs",
                                  "everything",
                                  "na.or.complete"),
                         .round = 4,
                         .stats = c("mean",
                                    "min",
                                    "max",
                                    "n",
                                    "sd",
                                    "sum"),
                         .prob = TRUE,
                         .nobs = TRUE,
                         .alpha = TRUE,
                         .plot = TRUE,
                         .output_path = NA_character_,
                         .output_file = NA_character_,
                         .warnings = FALSE) {
  # Parameter Setup and Checking
  # # Matching correlation missing handling
  .use <- match_arg_null(.use)
  # # Matching correlation coefficient method
  .method <- match_arg_null(.method)
  # # Ensuring that .stats contains acceptable values or is NULL
  # # If not explicitly passed, default to mean
  if (".stats" %in% names(match.call())) {
    .stats <- match_arg_null(.stats, .several_ok = TRUE)
  } else {
    .stats <- "mean"
  }

  # if .vars is not explicitly defined, set to be vector of all numeric except where
  # used in .with_vars and .byx
  if (is.null(.vars)) {
    .vars <- .dataset %>%
      # excluding any passed columns
      dplyr::select(-{{ .with_vars }}, -{{ .by_x }}) %>%
      # Keeping only numerics
      dplyr::select(tidyselect::where(is.numeric)) %>%
      names()
  }
  # create selection of data
  indata <- .dataset %>%
    dplyr::select({{ .vars }}, {{ .with_vars }}, {{ .by_x }})

  # Checking that all are numeric besides .by_x
  not_numeric <- indata %>%
    dplyr::select(-{{ .by_x }}) %>%
    purrr::map_lgl(~!is.numeric(.))
  if (any(not_numeric)) {
    stop(paste0("'", names(not_numeric)[which(not_numeric)], "'", collapse = ", "),
         " are not numeric")
  }

  # Displaying a warning if any constant and .warnings is TRUE
  if (.warnings) {
    constants <- indata %>%
      dplyr::select(-{{ .by_x }}) %>%
      purrr::map_lgl(~var(.) == 0)
    if (any(constants)) {
      warning(paste0("'", names(constants)[which(constants)], "'", collapse = ", "),
              " are constant")
    }
  }

  # If imputation requested, perform mean imputation on numeric
  # Also will force all integers into double
  if (.impute) {
    # Reconstruct indata on a column basis
    indata <- purrr::modify_if(
      indata,
      is.numeric,
      ~dplyr::if_else(is.na(.), mean(., na.rm = T), as.double(.)))
  }

  # If .by_x is specified, create iterable list of filtered dataframes .by_x
  # Otherwise, list of only unfiltered
  if (is.null(.by_x)) {
    indata_list <- list("No .by_x" = indata)
  } else {
    # Split into list .by_x, rows with missing .by_x implicitly dropped
    indata_list <- split(indata, dplyr::pull(indata, {{.by_x}})) %>%
      # Removing .by_x Variable from each dataframe
      purrr::modify(~dplyr::select(., -{{.by_x}}))
  }

  # Create Correlations for each element of indata_list
  corr_dfs <- indata_list %>%
    purrr::map(corr_extract_singleton,
               .vars = .vars,
               .with_vars = .with_vars,
               .warnings = .warnings,
               .method = .method,
               .use = .use,
               .prob = .prob,
               .nobs = .nobs,
               .alpha = .alpha,
               .stats = .stats,
               .round = .round)

  # Creating a excel output if requested
  if (!is.na(.output_path) & !is.na(.output_file)) {
    # Initializing workbook
    wb <- openxlsx::createWorkbook()
    # # Creating Style objects for workbooks
    correlation_style <- openxlsx::createStyle(
      border = c("top", "bottom", "left", "right"),
      borderColour = openxlsx::openxlsx_getOp("borderColour", "black"),
      borderStyle = openxlsx::openxlsx_getOp("borderStyle", "slantDashDot"),
      halign = "center",
      textRotation = -45
    )
    alpha_header_style <- openxlsx::createStyle(
      fontSize = 12,
      halign = "center",
      textDecoration = "bold",
      wrapText = FALSE)

    # setting up path to save the excel spreadsheet
    workbook_path <- paste0(.output_path, "/", .output_file, ".", "xlsx")
    # openxlsx works via side effects, creating a sheet for each .byx
    corr_dfs %>%
      purrr::iwalk(corr_extract_excel,
                   .vars = .vars,
                   .by_x = .by_x,
                   .alpha = .alpha,
                   .prob = .prob,
                   .wb = wb,
                   .plot = .plot,
                   .correlation_style = correlation_style,
                   .alpha_header_style = alpha_header_style)
    openxlsx::saveWorkbook(wb, workbook_path, overwrite = TRUE)
  }

  # Returning all correlations, binded by rows
  if (is.null(.by_x)) {
    corr_df <- corr_dfs[[1]]
  } else {
    corr_df <- corr_dfs %>%
      purrr::imap_dfr(~dplyr::mutate(.x, {{ .by_x }} := .y, .after = "Variable"))
  }

  return(corr_df)
}

#' Function to generate correlation extract on a particular value of \code{.byx} in
#' \link[mdrcAnalysis]{corr_extract}
#' @param .dataset Tibble to generate correlations on, only contains numeric data
#' @inheritParams corr_extract
#' @returns Dataframe containing correlation and, optionally, the univariate statistics
#' and Cronbach alphas based on passed parameters.
#' @keywords internal
corr_extract_singleton <- function(.dataset,
                                   .vars,
                                   .with_vars,
                                   .warnings,
                                   .method,
                                   .use,
                                   .prob,
                                   .nobs,
                                   .alpha,
                                   .stats,
                                   .round) {
  # Suppressing warnings in correlation calculation if desired
  if (.warnings) {
    corr_results <- psych::corr.test(.dataset, method = .method, use = .use, ci = FALSE)
  } else {
    corr_results <- suppressMessages(suppressWarnings(
      psych::corr.test(.dataset, method = .method, use = .use, ci = FALSE)))
  }
  # Coercing tibble from results based on requested output parameters
  corr_df <- tibble::as_tibble(corr_results$r, rownames = "Variable") %>%
    dplyr::mutate(Type = "Corr")
  # Fixing "bug" where correlation not always set to NA if sd is equal to 0.
  constant_cols <- .dataset %>%
    dplyr::select(tidyselect::where(~magrittr::or(var(., na.rm = TRUE) == 0,
                                                   is.na(var(., na.rm = TRUE))))) %>%
    names()
  corr_df %<>% dplyr::mutate(
    dplyr::across(dplyr::any_of(constant_cols), ~NA_real_),
    dplyr::across(!c(Variable, Type),
                  ~dplyr::if_else(Variable %in% constant_cols, NA_real_, .)))

  # Adding rows based on requested pairwise characteristics
  # # T-Score Probability
  if (.prob) {
    prob_df <- tibble::as_tibble(corr_results$p, rownames = "Variable") %>%
      dplyr::mutate(
        Type = "Prob",
        dplyr::across(!c(Variable, Type),
                      ~dplyr::if_else(Variable %in% constant_cols, NA_real_, .)))
    corr_df <- dplyr::bind_rows(corr_df, prob_df)
  }
  # # Number of Observations
  if (.nobs) {
    # Output from corr.test is a scalar numeric if all are identical
    if (is_scalar_double(corr_results$n)) {
      corr_results$n <- rep_len(corr_results$n,
                                magrittr::raise_to_power(nrow(corr_results$r), 2)) %>%
        matrix(nrow = nrow(corr_results$r)) %>%
        magrittr::set_colnames(colnames(corr_results$r)) %>%
        magrittr::set_rownames(rownames(corr_results$r))
    }
    n_df <- tibble::as_tibble(corr_results$n, rownames = "Variable") %>%
      dplyr::mutate(Type = "N")
    corr_df <- dplyr::bind_rows(corr_df, n_df)
  }
  # Performing .with_vars restructuring
  if (!is.null(.with_vars)) {
    corr_df <- corr_df %>%
      # Keeping all .vars
      dplyr::select(Variable, Type, {{ .vars }}) %>%
      # Removing all .vars rows
      dplyr::filter(!Variable %in% .vars)
  } else {
    # Otherwise reorder rows according to variable order
    corr_df <- corr_df %>%
      dplyr::arrange(factor(Variable, levels = .vars))
  }
  # Round values
  if (!is.null(.round)) {
    corr_df <- corr_df %>%
      dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~round(., .round)))
  }

  # Attaching Univariate Statistics as requested
  # # Selecting explicit rows present in case .with_vars removal occured
  if (!is.null(.stats)) {
    .dataset_long <- .dataset %>%
      dplyr::select(dplyr::all_of(unique(corr_df$Variable))) %>%
      tidyr::pivot_longer(cols = tidyselect::everything(),
                          names_to = "Variable",
                          values_to = "Value") %>%
      dplyr::group_by(Variable)
    univariates_df <- .stats %>%
      purrr::map(~dplyr::summarise(
        .dataset_long,
        {{ . }} := suppressWarnings(round(do.call(., list(Value, na.rm = T)), 3)),
        .groups = "drop")) %>%
      purrr::reduce(~dplyr::full_join(.x, .y, by = "Variable"))
    corr_df <- corr_df %>%
      dplyr::left_join(univariates_df, by = "Variable")
  }

  # Attaching Cronbach Alphas if desired
  if (.alpha) {

    # Skipping completely if only 1 or fewer non constant_cols
    if (length(setdiff(names(.dataset), constant_cols)) <= 1) {
      alpha_df <- tibble::tibble(Variable = names(.dataset)) %>%
        dplyr::mutate(Raw_Alpha = NA_real_, Std_Alpha = NA_real_)
      alpha_results <- list(total = c(NA_real_, NA_real_))
    } else {
      # # Suppressing warnings if desired
      if (.warnings) {
        alpha_results <- purrr::quietly(psych::alpha)(.dataset)
        warning(alpha_results$warnings)
        alpha_results <- alpha_results$result
      } else {
        alpha_results <- purrr::quietly(psych::alpha)(.dataset)$result
      }
      alpha_df <- alpha_results %>%
        magrittr::use_series("a") %>%
        tibble::as_tibble(rownames = "Variable") %>%
        dplyr::select(Variable, Raw_Alpha = raw_alpha, Std_Alpha = std.alpha)
    }

    corr_df <- corr_df %>%
      dplyr::left_join(alpha_df, by = "Variable") %>%
      dplyr::mutate(Total_Raw_Alpha = round(alpha_results$total[[1]], 4),
                    Total_Std_Alpha = round(alpha_results$total[[2]], 4))
  }
  # Reordering
  corr_df <- corr_df %>%
    dplyr::relocate(Variable, Type, dplyr::matches("Alpha")) %>%
    dplyr::relocate({{.vars}}, .after = dplyr::last_col())
  return(corr_df)
}

#' Creates excel output from \link[mdrcAnalysis]{corr_extract_singleton}
#' @inheritParams corr_extract
#' @param .datset Output from \link[mdrcAnalysis]{corr_extract_singleton}
#' @param .byval Specific value of .byx from \link[mdrcAnalysis]{corr_extract}
#' @param .wb Workbook to write to
#' @param .correlation_style Default styling for correlation table
#' @param .alpha_header_style Styling for alphas grouped header
#' @returns NULl
#' @keywords internal
corr_extract_excel <- function(.dataset,
                               .by_val,
                               .vars,
                               .by_x,
                               .alpha,
                               .prob,
                               .wb,
                               .plot,
                               .correlation_style,
                               .alpha_header_style) {
  # Clearing pairwise values where column and variable name is matching or right of diag
  corr_df <- .dataset %>%
    dplyr::mutate(dplyr::across(
      {{ .vars }},
      ~dplyr::if_else(dplyr::row_number() <= which.max(Variable == dplyr::cur_column()),
                      NA_real_, .)
    ))
  # Generating sheet name based on .by_val if .by_x was defined
  excelsn <- dplyr::if_else(is.null(.by_x),
                            "Correlations",
                            paste("Corr:", paste(.by_x, .by_val, sep = " = ")))

  # Truncating if name is too long
  excelsn <- stringr::str_trunc(excelsn, 80, "center")

  # Creating Sheet
  openxlsx::addWorksheet(.wb, sheetName = excelsn)
  openxlsx::mergeCells(.wb, excelsn, cols = 1:ncol(corr_df), rows = 1)
  openxlsx::writeData(.wb, excelsn, "Corrextract Output", startCol = 1, startRow = 1)
  # # If Alphas were calculated, extract total alphas and place as header
  if (.alpha) {
    total_raw_alpha <- corr_df %>%
      dplyr::pull(Total_Raw_Alpha) %>%
      unique()
    corr_df %<>% dplyr::select(-Total_Raw_Alpha)
    total_std_alpha <- corr_df %>%
      dplyr::pull(Total_Std_Alpha) %>%
      unique()
    corr_df %<>% dplyr::select(-Total_Std_Alpha)

    header2 <- paste0("Cronbach's Alpha   Raw = ", total_raw_alpha,
                      "    Standardized = ", total_std_alpha)
    openxlsx::mergeCells(.wb, excelsn, cols = 1:ncol(corr_df), rows = 2)
    openxlsx::writeData(.wb, excelsn, header2, startCol = 1, startRow = 2)
    openxlsx::addStyle(.wb, excelsn, style = .alpha_header_style, rows = 1:2, cols = 1)
  }
  # Writing correlations to sheet
  corr_rowstart <- 2 + .alpha
  openxlsx::writeData(.wb, sheet = excelsn, x = corr_df, startRow = corr_rowstart)
  start_col <- which.max(colnames(corr_df) %in% .vars)
  # Adding styling to correlations
  openxlsx::addStyle(.wb,
                     excelsn,
                     style = .correlation_style,
                     corr_rowstart,
                     start_col:ncol(corr_df))

  # Creating row number vector proxy for excel sheet to style correlations
  corr_rows <- corr_df %>%
    dplyr::mutate(rn = dplyr::row_number() + corr_rowstart) %>%
    dplyr::filter(Type == "Corr") %>%
    dplyr::pull(rn)
  # Per-Cell Styling
  for (i in corr_rows) {
    for (j in c(start_col:ncol(corr_df))) {
      row <- i - corr_rowstart
      col <- j
      value <- corr_df[row, col, drop = TRUE]
      pvalue <- ifelse(.prob,
                       corr_df[row + 1, col, drop = TRUE], 99)
      color <- dplyr::case_when(
        value < -0.75 ~ "#CC0004",
        value < -0.5 ~ "#F20000",
        value < -0.25 ~ "#FF7171",
        value < 0.0 ~ "#FFB9B9",
        value == 0.0 ~ "#ffffff",
        value < 0.25 ~ "#B4D79D",
        value < 0.5 ~ "#85BE5E",
        value < 0.75 ~ "#6BA743",
        TRUE ~ "#517D33")
      if (is.na(value)){
        next

      }
      if (pvalue <= 0.1) {
        newStyle <- openxlsx::createStyle(fgFill = color, textDecoration = "bold")
      } else {
        newStyle <- openxlsx::createStyle(fgFill = color)
      }

      openxlsx::addStyle(
        wb = .wb,
        sheet = excelsn,
        style = newStyle,
        cols = j,
        rows = i)
    } }

  # Styling .prob if created
  if (.prob) {
    prob_rows <- corr_df %>%
      dplyr::mutate(rn = dplyr::row_number() + corr_rowstart) %>%
      dplyr::filter(Type == "Prob") %>%
      dplyr::pull(rn)
    for (i in prob_rows) {
      # find first NA column
      stop_col <- which.max(is.na(corr_df[i - corr_rowstart, ])) - 1
      stop_col <- ifelse(is.na(stop_col), ncol(corr_df) - 1, stop_col)
      # as long as the first column isn't NA ...
      if (!is.na(corr_df[i - corr_rowstart, start_col]) & !is.na(stop_col)){
        openxlsx::conditionalFormatting(
          wb = .wb,
          sheet = excelsn,
          cols = start_col:stop_col,
          rows = i,
          rule = c(0, 0.01),
          type = "between",
          style = openxlsx::createStyle(bgFill = '#FFFA00')
        )}
      # Rules defined this way because between is inclusive
      openxlsx::conditionalFormatting(
        wb = .wb,
        sheet = excelsn,
        cols = start_col:ncol(corr_df),
        rows = i,
        rule = c(0.01 + 0.000000001, 0.05),
        type = "between",
        style = openxlsx::createStyle(bgFill = '#FFFF89')
      )

      openxlsx::conditionalFormatting(
        wb = .wb,
        sheet = excelsn,
        cols = start_col:ncol(corr_df),
        rows = i,
        rule = c(0.05 + 0.000000001,0.10),
        type = "between",
        style = openxlsx::createStyle(bgFill = '#FFFFC5')
      )
    }
  }

  # Adding filter for Typing
  openxlsx::addFilter(.wb, excelsn, rows = corr_rowstart, cols = 2)

  # Creating correlation plot if requested
  if (.plot) {
    plotsn <- dplyr::if_else(is.null(.by_x),
                             "Correlation Plot",
                             paste("Corr Plot:", paste(.by_x, .by_val, sep = " = ")))
    plotsn <- stringr::str_trunc(plotsn, 80, "center")
    openxlsx::addWorksheet(.wb, sheetName = plotsn)

    # Creating and storing "current plot"
    .dataset %>%
      dplyr::filter(Type == "Corr") %>%
      tibble::column_to_rownames("Variable") %>%
      dplyr::select(dplyr::all_of(.vars)) %>%
      as.matrix() %>%
      corrplot::corrplot()
    grDevices::recordPlot()

    openxlsx::insertPlot(.wb, sheet = plotsn, startRow = 6, startCol = 2)
  }
  return(NULL)
}

