#' Explanatory Factor Analysis Extract
#'
#' Generate exploratory factor analysis extract.
#' Removes collinear variables from analysis and displays Cronbach alphas, EFA
#' loadings, and EFA summary. EFA loadings are color-coded based on correlation
#' strength
#'
#' @param .data Tibble of data
#' @param .variable_names Character vector specifying column subset for EFA
#' @param .variable_labels Optional Character vector specifying verbose columns
#' for \code{.variable_names}
#' @param .cor_use String parameter passed to stats::cor.
#' Defaults to "pairwise.complete.obs"
#' @param .cor_method String parameter passed to stats::cor.
#' Defaults to "pearson"
#' @param .cor_limit Numeric designating upper bound of varaible correlation
#' until variable is removed because of collinearity. None removed if NA.
#' Choice of variable removed is first listed.
#' @param .fa_nfactors Integer vector specifying number of factors to extract OR
#' string "ncol" if want to be equal to number of columns used in EFA
#' @param .fa_rotate String of rotations type to use in psych::fa
#' @param .response_list List of character vectors specifying response variables
#' @param .file_path Path of output excel extract
#' @param .file_name Name of output file extract
#' @param ... Additional arguments passed to psych::fa
#' @return Primarily run for the side effect of creating an excel workbook of
#' efa extracts for each .response_list passed. Also returns a list of output
#' @export
#' @examples
#' efa_extract(bfi, .variable_names = c("A1", "A2", "A3", "A4", "A5"))
#' efa_extract(bfi, .variable_names = c("A1", "A2", "A3", "A4", "A5"),
#'            .file_path = "./", .file_name = "set_a_multi_factor")
efa_extract <- function(
    .data,
    .variable_names,
    .variable_labels = NULL,
    .cor_use = c("pairwise.complete.obs",
                 "complete.obs",
                 "all.obs",
                 "everything",
                 "na.or.complete"),
    .cor_method = c("pearson", "spearman"),
    .cor_limit = .8,
    .fa_nfactors = 1L,
    .fa_rotate = c("varimax", "none", "quartimax", "bentlerT",
                   "equamax", "varimin", "geominT", "bifactor",
                   "Promax", "promax", "oblimin", "simplimax",
                   "betlerQ", "geominQ", "biquartimin", "cluster"),
    .file_path,
    .file_name,
    ...) {

  # Parameter Set-Up

  # # Checking variable_labels length equivalent to names
  if (!is.null(.variable_labels)) {
    if (length(.variable_names) != length(.variable_labels)) {
      stop(".variable_labels must match .variable_names length")
    }
  }
  # # Passing arguments for stats::cor
  .cor_use <- match.arg(.cor_use)
  .cor_method <- match.arg(.cor_method)
  .cor_limit <- as.numeric(.cor_limit)
  # # Arguments for psych::fa
  if (!(all(.fa_nfactors == "ncol") | is.numeric(.fa_nfactors))) {
    stop(".fa_nfactors must be either 'ncol' or an integer vector")
  }
  .fa_rotate <- match.arg(.fa_rotate)
  # # Making n_factors into a named vector for output naming
  names(.fa_nfactors) <- .fa_nfactors
  # Other ellipsis will be passed to psych::fa
  fa_args <- list(...)
  # Initializing Excel File
  wb <- openxlsx::createWorkbook()

  # # Create Styling Used for EFA Tables
  style1 <- openxlsx::createStyle(fgFill = "#cae3eb", fontColour = "white")
  style2 <- openxlsx::createStyle(fgFill = "#63b6cf", fontColour = "white")
  style3 <- openxlsx::createStyle(fgFill = "#00558C", fontColour = "white")
  style4 <- openxlsx::createStyle(fgFill = "#002B49", fontColour = "white")
  style_bold <- openxlsx::createStyle(textDecoration = c("bold", "underline"))

  # Generating data with only relevant vectors
  .df <- .data %>%
    dplyr::select({{ .variable_names }})

  # Create sheet for each element in nfactors
  efa_list <- purrr::map(.fa_nfactors, function(.fa_nfactor) {
    # Printing number of factors for diagnostic purposes
    cat(.fa_nfactor, "number of factor solutions\n")
    # Generating sheet name from .fa_nfactor
    .sheet_name <- paste(.fa_nfactor, "Factors")
    # Adding Sheet with name
    openxlsx::addWorksheet(wb, .sheet_name)


    # Generating Correlation matrix to check for collinearity
    .matrix_corr <- stats::cor(.df,  use = .cor_use,  method = .cor_method)
    # Generating factors that should probably be reversed
    .df_principal <- psych::principal(.df, scores = FALSE)
    # Generating alphas dataframes
    .df_alpha <- psych::alpha(.df, warnings = FALSE)
    .df_alpha_total <- .df_alpha$total %>%
      dplyr::select(raw_alpha, std.alpha, average_r, ase, mean, sd)
    .df_alpha_drop <- .df_alpha$alpha.drop %>%
      tibble::as_tibble(rownames = "Item") %>%
      dplyr::select(Item, raw_alpha, std.alpha, average_r, `alpha se`)
    # # Adding Labels if possible
    if (!is.null(.variable_labels)) {
      .df_alpha_drop <- .df_alpha_drop %>%
        dplyr::mutate(Label = .variable_labels, .after = "Item")
    }
    # Determining collinear variables
    if (!is.na(.cor_limit)) {
      # Create working matrix based on correlations
      collinearity_matrix <- .matrix_corr

      # Setting all correlations at or above diagonal to 0
      collinearity_matrix[upper.tri(collinearity_matrix, diag = T)] <- 0
      # Converting to named list
      collinearity_list <- tibble::as_tibble(collinearity_matrix) %>%
        as.list()
      # Identifying collinear pairs
      collinear_pairs <- collinearity_list %>%
        # Taking all rows that are colinear
        purrr::map(~.variable_names[abs(.x) > .cor_limit]) %>%
        purrr::keep(~length(.x) > 0)
      # Take collinear column names
      # Creating dataframe without collinear columns
      if (length(collinear_pairs) > 0) {
        # Removing labels as well if provided
        if ((!is.null(.variable_labels))) {
          .variable_labels <- .variable_labels %>%
            magrittr::extract(!.variable_names %in% names(collinear_pairs))
        }
        .df <- .df %>%
          dplyr::select(-dplyr::all_of(names(collinear_pairs)))
        # removing any .fa_nfactors that exceed new maximal number
        if (is.numeric(.fa_nfactor) & .fa_nfactor > ncol(.df)) {
          warning(paste("Attempting to fit too many factor solutions because",
                        paste(names(collinear_pairs), collapse = ", "),
                        "removed because of collinearity"))
          efa_list <- list("Hindrance_Alphas" = .df_alpha_total,
                           "Drop_Alphas" = .df_alpha_drop)
          openxlsx::writeData(
            wb, sheet = .sheet_name,
            x = paste("Variable removed due to collinearity. Insufficient for",
                      .fa_nfactor, "solution"))
          .row_pointer <- 2L
          openxlsx::writeData(wb, sheet = .sheet_name,
                              x = "Variables removed for collinearity:",
                              startRow = .row_pointer)
          openxlsx::writeData(wb, sheet = .sheet_name,
                              x = "Collinear with:",
                              startRow = .row_pointer, startCol = 2)
          .row_pointer <- .row_pointer + 1L
          for (collinear_colindex in seq_along(collinear_pairs)) {
            openxlsx::writeData(wb, sheet = .sheet_name,
                                x = names(collinear_pairs)[collinear_colindex],
                                startRow = .row_pointer)
            openxlsx::writeData(wb, sheet = .sheet_name,
                                x = paste(collinear_pairs[collinear_colindex],
                                          collapse = ", "),
                                startRow = .row_pointer, startCol = 2)
            .row_pointer <- .row_pointer + 1L
          }
          return(efa_list)
        }
      }
    } else{
      collinear_pairs <- list()
    }

    # Exploratory Factory Analysis, list to all of for multiple nfactors
    if (all(.fa_nfactor == "ncol")) {
      # Using do.call to pass other ... arguments from initial function call
      .df_fa <- do.call(psych::fa,
                        args = c(list(
                          "r" = .df,
                          "nfactors" = ncol(.df),
                          "rotate" = .fa_rotate),
                          fa_args))
      # # Generating Scree diagram dataframe
      .df_scree <- dplyr::tibble("Factor_n" = seq_len(length(.df_fa$e.values)),
                                 "Eigenvalues" = .df_fa$e.values)
    } else {
      .df_fa <- do.call(psych::fa,
                        args = c(list(
                          "r" = .df,
                          "nfactors" = .fa_nfactor,
                          "rotate" = .fa_rotate),
                          fa_args))
      # Performing efa on all columns for scree plot
      .df_fa_scree <- do.call(psych::fa,
                              args = c(list(
                                "r" = .df,
                                "nfactors" = ncol(.df),
                                "rotate" = .fa_rotate),
                                fa_args))
      # # Generating Scree diagram dataframe
      .df_scree <- dplyr::tibble(
        "Factor_n" = seq_len(length(.df_fa_scree$e.values)),
        "Eigenvalues" = .df_fa$e.values)
    }

    # Creating Extracts in Excel Sheet
    .row_pointer <- 1L

    # # If any collinear columns, print removed columns
    if (length(collinear_pairs) > 0) {
      openxlsx::writeData(wb, sheet = .sheet_name,
                          x = "Variables removed for collinearity:",
                          startRow = .row_pointer)
      openxlsx::writeData(wb, sheet = .sheet_name,
                          x = "Collinear with:",
                          startRow = .row_pointer, startCol = 2)
      .row_pointer <- .row_pointer + 1L
      for (collinear_colindex in seq_along(collinear_pairs)) {
        openxlsx::writeData(wb, sheet = .sheet_name,
                          x = names(collinear_pairs)[collinear_colindex],
                          startRow = .row_pointer)
        openxlsx::writeData(wb, sheet = .sheet_name,
                            x = paste(collinear_pairs[collinear_colindex],
                                      collapse = ", "),
                            startRow = .row_pointer, startCol = 2)
        .row_pointer <- .row_pointer + 1L
      }
      .row_pointer <- .row_pointer + 2L
    }
    # # Display warning if any factors should be reversed
    if (any(.df_principal$loadings < 0)) {
      warning("Some items (",
              rownames(.df_principal$loadings)[(.df_principal$loadings <  0)],
              ") were negatively correlated with the total scale and",
              "probably should be reversed.")
      openxlsx::writeData(
        wb, sheet = .sheet_name,
        x = paste("Some items (",
                  rownames(.df_principal$loadings)[.df_principal$loadings < 0],
                  ") were negatively correlated with the total scale and",
                  "probably should be reversed."),
        startRow = .row_pointer)
      .row_pointer <- .row_pointer + 2L
    }
    # # Printing Alphas
    # # # Hindrances to Scale Alpha
    # # # # Subtitle
    openxlsx::writeData(wb, sheet = .sheet_name,
                        x = "Hindrances to Scale Alpha",
                        startRow = .row_pointer)
    openxlsx::addStyle(wb, sheet = .sheet_name, style_bold,
                       rows = .row_pointer, cols = 1)
    .row_pointer <- .row_pointer + 2L
    # # # # Table
    openxlsx::writeDataTable(wb, sheet = .sheet_name,
                             x = .df_alpha_total, startRow = .row_pointer)
    .row_pointer <- .row_pointer + nrow(.df_alpha_total) + 3L
    # # # Alpha Drop Table
    # # # # Subtitle
    openxlsx::writeData(wb, sheet = .sheet_name,
                        x = "Alpha Drop Table", startRow = .row_pointer)
    openxlsx::addStyle(wb, sheet = .sheet_name, style_bold,
                       rows = .row_pointer, cols = 1)
    .row_pointer <- .row_pointer + 2L
    # # # # Table
    openxlsx::writeDataTable(wb, sheet = .sheet_name,
                             x = .df_alpha_drop, startRow = .row_pointer)
    .row_pointer <- .row_pointer + nrow(.df_alpha_drop) + 3L
    # # Printing Scree Plot and Eigenvalues
    # # # Eigenvalues
    openxlsx::writeData(wb, sheet = .sheet_name,
                        x = "Scree Plot", startRow = .row_pointer)
    openxlsx::addStyle(wb, sheet = .sheet_name, style_bold,
                       rows = .row_pointer, cols = 1)
    .row_pointer <- .row_pointer + 2L
    openxlsx::writeDataTable(wb, sheet = .sheet_name,
                             x = .df_scree, startRow = .row_pointer)
    # # # Scree Plot
    scree_plot <- ggplot2::ggplot(
      .df_scree,
      ggplot2::aes(x = Factor_n, y = Eigenvalues, group = 1)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::xlab("Number of Factors") +
      ggplot2::ylab("Initial Eigenvalue") +
      ggplot2::labs(title = "Scree Plot")
    print(scree_plot)
    grDevices::recordPlot()
    openxlsx::insertPlot(wb, sheet = .sheet_name,
                         startRow = .row_pointer - 10, startCol = 10)
    .row_pointer <- .row_pointer + nrow(.df_scree) + 2L
    # # Printing Each factor Analysis Produced

    # # # Extacting and formatting loadings
    .df_fa_loading <- .df_fa$loadings %>%
      unclass() %>%
      tibble::as_tibble(rownames = "item") %>%
      dplyr::mutate(dplyr::across(tidyselect:::where(is.numeric), ~round(., 2)))
    if (!is.null(.variable_labels)) {
      .df_fa_loading <- .df_fa_loading %>%
        dplyr::mutate(Label = .variable_labels, .after = 1)
    }
    # # # Subtitle
    openxlsx::writeData(wb, sheet = .sheet_name,
                        x = paste("Factor Analysis"),
                        startRow = .row_pointer)
    openxlsx::addStyle(wb, sheet = .sheet_name, style_bold,
                       rows = .row_pointer, cols = 1)
    .row_pointer <- .row_pointer + 2L

     # # Identifying if any loadings should be reversed  because negative
    # # # Print Loadings
    openxlsx::writeDataTable(wb, sheet = .sheet_name,
                             x = .df_fa_loading,
                             startRow = .row_pointer)
    # # # # Applying Formatting
    # # # # # Mutating loading tibble to help with tracking rows
    .df_fa_loading_format <- .df_fa_loading %>%
      dplyr::select(-item) %>%
      dplyr::mutate(rn = dplyr::row_number()) %>%
      dplyr::relocate(rn)
    if (is.null(.variable_labels)) {
      names(.df_fa_loading_format) <- c(
        "rn", seq_len(ncol(.df_fa_loading_format) - 1L) + 1L)
    } else{
      names(.df_fa_loading_format) <- c(
        "rn", "Label", seq_len(ncol(.df_fa_loading_format) - 2L) + 2L)
    }
    purrr::pwalk(.df_fa_loading_format, function(rn, ...) {
      row_values <- list(...) %>%
        magrittr::extract(names(.) != "Label") %>%
        # Changing arguments into a named numeric vector
        unlist()
      names(row_values)[row_values < 0.2] %>%
        purrr::walk(~openxlsx::addStyle(wb, sheet = .sheet_name, style = style1,
                                        rows = as.integer(rn) + .row_pointer,
                                        cols = as.integer(.)))
      names(row_values)[row_values >= 0.2 & row_values < 0.4] %>%
        purrr::walk(~openxlsx::addStyle(wb, sheet = .sheet_name, style = style2,
                                        rows = as.integer(rn) + .row_pointer,
                                        cols = as.integer(.)))
      names(row_values)[row_values >= 0.4 & row_values < 0.6] %>%
        purrr::walk(~openxlsx::addStyle(wb, sheet = .sheet_name, style = style3,
                                        rows = as.integer(rn) + .row_pointer,
                                        cols = as.integer(.)))
      names(row_values)[row_values >= 0.6] %>%
        purrr::walk(~openxlsx::addStyle(wb, sheet = .sheet_name, style = style4,
                                        rows = as.integer(rn) + .row_pointer,
                                        cols = as.integer(.)))

    })
    # # # Print FA Plot
    psych::fa.diagram(.df_fa)
    grDevices::recordPlot()
    openxlsx::insertPlot(wb, sheet = .sheet_name,
                         height = max(4, nrow(.df_fa_loading_format) * 0.3),
                         width = max(4, nrow(.df_fa_loading_format) * 0.3),
                         startRow = .row_pointer,
                         startCol = max(12, ncol(.df_fa_loading_format) + 5))
    .row_pointer <- .row_pointer + nrow(.df_fa_loading) + 3L
    # # # Printing what accounts for variance
    .df_fa_vaccounted <- .df_fa$Vaccounted %>%
      tibble::as_tibble(rownames = "Variance Accounted")
    openxlsx::writeDataTable(wb, sheet = .sheet_name,
                             x = .df_fa_vaccounted,
                             startRow = .row_pointer)
    .row_pointer <- .row_pointer + nrow(.df_fa_vaccounted) + 1L
    # # # Printing Factor Analysis Summary
    # # # # Extracting from pre-formatted print of fa object
    .df_fa_captured <- invisible(utils::capture.output(.df_fa))
    # # # # Extracting strings desired for printing
    .df_fa_summary_extract <- c(
      # Description of call, manaully generated for a more intuitive print
      paste("Variables:", paste(.variable_names, collapse = ", ")),
      paste("Rotation:", .fa_rotate),
      paste("N Factor:", .fa_nfactor),
      paste("Other Arguments:",
            paste(purrr::imap_chr(fa_args, ~paste0(.y, ": ", .x)),
                  collapse = ", ")),
      # Other information
      .df_fa_captured[
        seq(stringr::str_which(.df_fa_captured, "Mean item complexity"),
            length(.df_fa_captured))]
    )
    # # # # Adding to Excel, each string is a separate row
    for (string in .df_fa_summary_extract) {
      openxlsx::writeData(wb, sheet = .sheet_name,
                          x = string,
                          startRow = .row_pointer)
      .row_pointer <- .row_pointer + 1L
    }
    # Line break for console
    cat("\n")
    # Generating elements to return
    efa_list <- list("Hindrance_Alphas" = .df_alpha_total,
                     "Drop_Alphas" = .df_alpha_drop, "Eigenvalues" = .df_scree,
                     "EFA" = .df_fa_loading)
  })

  # Generating File name to write to file
  if (!missing(.file_path) & !missing(.file_name)) {
    .file_name <- dplyr::if_else(stringr::str_detect(.file_name, "\\.xlsx$"),
                                 .file_name, paste0(.file_name, ".xlsx"))
    openxlsx::saveWorkbook(wb, file = file.path(.file_path, .file_name),
                           overwrite = T)
  }

  # Invisibly returning list of output
  invisible(efa_list)
}
