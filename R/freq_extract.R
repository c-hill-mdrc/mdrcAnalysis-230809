#' Generates frequency tables from a given dataset
#'
#' @param .dataset Required dataframe, Input dataset.
#' @param .freq_vars Character vector of variables to process. Defaults to all
#' variables with fewer than .max_lvls unique values
#' @param .subgroup Optional string, Variable name to group frequency by.
#' @param .wt_var Optional string, Variable name to weight observations by.
#' @param .max_lvls Optional integer, Used in automatic .freq_vars selection
#' @param .na_rm Logical, Indicates whether NA values are dropped from
#' cumulative calculations.
#' @param .round Optional integer, Decimal places to round to.
#' @param .output_path Optional string, File path of output excel sheet.
#' @param .output_file Optional string, File name of output excel sheet.
#'
#' @returns Dataframe of frequency table based on input specifications.
#'
#' @author Audrey Yu
#'
#' @export freq_extract
#'
#' @examples
#' freq_extract(iris, c("Species", "Sepal.Length"))
#' freq_extract(iris)
#' freq_extract(mtcars, .subgroup = "cyl", .max_lvls = 5, .na_rm = TRUE)
freq_extract <- function(.dataset,
                         .freq_vars,
                         .max_lvls = 24,
                         .subgroup,
                         .wt_var,
                         .na_rm = FALSE,
                         .round = 2,
                         .output_path = NA_character_,
                         .output_file = NA_character_) {
  # Extracting arguments as a list to pass
  args_list <-
    match.call.defaults() %>%
    as.list()

  # Setting .freq_vars to be all variables with fewer than .max_lvls unique values
  .freq_vars <-
    rlang::maybe_missing(.freq_vars,
                         .dataset %>%
                           dplyr::select(
                             tidyselect::where(~dplyr::n_distinct(.)
                                               <= .max_lvls))
                         %>%
                           names()
                         )

  # Excluding .subgroup and .wt_var if not missing
  if (!rlang::is_missing(.subgroup)) {
    .freq_vars <- .freq_vars[.freq_vars != .subgroup]
  }

  if (!rlang::is_missing(.wt_var)) {
    .freq_vars <- .freq_vars[.freq_vars != .wt_var]
  }

  # # Ensuring that .round is positive integer
  .round <- abs(floor(.round))

  # Writing function to iterate over
  # As a package, we could possibly write this as an internal "non-exported"
  # function and keep it in the package environment
  # this would keep it defined outside of the function
  # Rigorous argument checks not needed internally
  freq_extract_singleton <- function(.var, .dataset, .subgroup, .wt_var, .na_rm){

    # Grouping with subgroups if applies
    if (missing(.subgroup)) {
      data <- .dataset %>%
        dplyr::group_by(dplyr::across(tidyselect::all_of(.var)))
    } else{
      data <- .dataset %>%
        dplyr::group_by(dplyr::across(tidyselect::all_of(c(.subgroup, .var))))
    }

    # Adding weights if applies
    if (!missing(.wt_var)) {
      data <- data %>%
        dplyr::summarise(Freq = dplyr::n(),
                         Wt_Freq = sum(!!ensym(.wt_var)),
                         .groups = "drop")
    } else {
      data <- data %>%
        dplyr::summarise(Freq = dplyr::n(),
                         .groups = "drop")
    }

    # Making variable value and subgroup column, also unique row id
    data <- data %>%
      dplyr::mutate(Variable = .var,
                    Value = !!ensym(.var),
                    .before = 1)
    if (!missing(.subgroup)) {
      data <- data %>%
        dplyr::mutate(Subgrp_Var   = .subgroup,
                      Subgrp_Value = !!ensym(.subgroup),
                      .before = 1) %>%
        dplyr::mutate(Unique_ID  = paste(Subgrp_Value,
                                         Variable,
                                         Value,
                                         sep = ":"),
                      .before = 1)
    } else {
      data <- data %>%
        dplyr::mutate(Unique_ID  = paste(Variable, Value, sep = ":"),
                      .before = 1)
    }

    # Removing Unneeded columns
    if (missing(.subgroup)) {
      data <- data %>%
        dplyr::select(-(tidyselect::all_of(.var)))
    } else {
      data <- data %>%
        dplyr::select(-tidyselect::all_of(c(.var, .subgroup)))
    }

    # Performing other frequency summaries
    if (missing(.subgroup)) {
      data <-
        data %>%
        dplyr::mutate(
        ##non-weight variables
        #when .na_rm=no
        CFreq = cumsum(Freq),
        Pct = ((Freq) / sum(Freq)) * 100,
        CPct  = (CFreq / max(CFreq)) * 100
      )
    } else {
      data <-
        data %>%
        dplyr::group_by(Subgrp_Value) %>%
        dplyr::mutate(
          ## non-weight variables
          #  when .na_rm=no
          CFreq = cumsum(Freq),
          Pct = ((Freq) / sum(Freq)) * 100,
          CPct  = (CFreq / max(CFreq)) * 100
          )
    }

    # If .na_rm true, want to ignore where missing
    if (.na_rm) {
      data <- data %>%
        dplyr::mutate(
          Nm     = ifelse(is.na(Value) & .na_rm, Freq, 0L),
          Nm_Freq = Freq - Nm,
          Nm_CFreq  = cumsum(Nm_Freq),
          Nm_Pct = (Nm_Freq / sum(Nm_Freq)) * 100,
          Nm_CPct =  (Nm_CFreq / max(Nm_CFreq)) * 100
        )
    }

    # Weight Variables if na_rm true
    if (!missing(.wt_var)) {
      if(.na_rm) {
        data <- data %>%
          dplyr::mutate(
            ##weight variables
            #when .na_rm=no
            Wt_CFreq = cumsum(Wt_Freq),
            Wt_Pct = ((Wt_Freq) / sum(Wt_Freq)) * 100,
            Wt_CPct  = (Wt_CFreq / max(Wt_CFreq)) * 100,
            Wt_Nm     = ifelse(is.na(Value) & .na_rm, Wt_Freq, 0L),
            Wt_Nm_Freq = Wt_Freq - Wt_Nm,
            Wt_Nm_CFreq  = cumsum(Wt_Nm_Freq),
            Wt_Nm_Pct = ( Wt_Nm_Freq / sum(Wt_Nm_Freq)) * 100,
            Wt_Nm_CPct = (Wt_Nm_CFreq / max(Wt_Nm_CFreq)) * 100
          )
      } else {
        data <- data %>%
          dplyr::mutate(
            ##weight variables
            #when .na_rm=no
            Wt_CFreq = cumsum(Wt_Freq),
            Wt_Pct = ((Wt_Freq) / sum(Wt_Freq)) * 100,
            Wt_CPct  = (Wt_CFreq / max(Wt_CFreq)) * 100
          )
      }
    }
    # Cleaning columns
    data <- data %>%
      # Making all non-numeric into character
      dplyr::mutate(
        dplyr::across(dplyr::any_of(c("Value", "Subgrp_Value")), as.character),
        dplyr::across(tidyselect::where(is.numeric), ~round(., .round))
      )
    return(data)
  }

  # Take matching arguments
  args_list_matched <- args_list %>%
    magrittr::extract(names(.) %in% rlang::fn_fmls_names(freq_extract_singleton))
  # Iterating through each .freq_vars
  frequency_output <- .freq_vars %>%
    purrr::map_dfr(function(.var, .args_list_matched){
      args_list_matched <- c(list(".var" = .var), .args_list_matched)
      do.call(freq_extract_singleton, args_list_matched)
    }, args_list_matched)

  if (!missing(.output_path) & !missing(.output_file)) {
    create_excel(.output_path = .output_path,
                 .output_file = .output_file, frequency_output)
  }

  return(frequency_output)
}



