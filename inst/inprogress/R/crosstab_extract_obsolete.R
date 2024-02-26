#' Crosstab main function
#'
#' @param .dataset Input dataframe
#' @param .col_var Required, string of "columns" of cross tab
#' @param .subgroup String of column to use as by variable
#' @param .missing Logical, keep missing in data
#' @param .l_row_var list of Variables to cross by with the single ColVar
#' @param .max_lvls maximum number of levels for a categorical variable
#' @param .l_row_label optional list of Row Variables to go by
#' @param .output_path path to output location
#' @param .output_file filename for output
#'
#' @returns A tibble of all the cross-tabulations

#' @export crosstab_extract
#'
#' @examples
#' crosstab_extract(
#' .dataset = bif_final,
#' .col_var = "COHORT",
#' .l_row_var = c("blmale", "bldiplomas_as")
#' )
#'
crosstab_extract <- function(
  .dataset,
  .col_var,
  .l_row_var = NULL,
  .l_row_label = NULL,
  .subgroup = NULL,
  .missing = FALSE,
  .max_lvls = 24,
  .output_path = NA_character_,
  .output_file = NA_character_
){

  # convert everything to "NA"
  if(.missing){
  .dataset <- .dataset %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),
                                ~ifelse(is.na(.), "NA", .)))
  }

  # if no list of row variables is provided,
  # and if subgroups are not specified,
  # pick all the variables that is not specified as col var to cross tab by
  # set a max level limit too. Otherwise, the function won't get out of a loop!
  if (missing(.l_row_var)) {

    if (is.null(.subgroup)) {
      .l_row_var <- names(.dataset)[!names(.dataset) %in% .col_var]

      # limiting the maximum levels
      .l_row_var <- .dataset %>%
        dplyr::select(all_of(.l_row_var)) %>%
        dplyr::select(tidyselect::where(
          ~dplyr::n_distinct(., na.rm = !.missing) <= .max_lvls)) %>%
        names()
    } else {
      pre_selected <- c(.col_var, .subgroup)
      .l_row_var <- names(.dataset)[!names(.dataset) %in% pre_selected]


      # limiting the maximum level
      .l_row_var <- .dataset %>%
        dplyr::select(all_of(.l_row_var)) %>%
        dplyr::select(tidyselect::where(
          ~dplyr::n_distinct(., na.rm = !.missing) <= .max_lvls)) %>%
        names()

    }
  }

  # .subgroup set, extract all unique values
  if (!missing(.subgroup) & !is.null(.subgroup)) {
     .subgroup_values <- .dataset %>%
       dplyr::pull({{.subgroup}}) %>%
       unique()
  } else{
     .subgroup_values <- NULL
  }

  # grid of all possible parameters
    parameter_grid <- tidyr::expand_grid(
      ".row_var" = .l_row_var,
      ".col_var" = .col_var,
      ".subgroup_values"= .subgroup_values

    )

    # we do not need to iterate through a case where subgroup_values are missing
    # or any of the variant

    if(.missing == FALSE){

      parameter_grid <- na.omit(parameter_grid)

    } else {

      parameter_grid <-
        parameter_grid %>%
        # Coercing Actual NA values to the string "NA", making all character
        dplyr::mutate(dplyr::across(dplyr::everything(),
                                    ~ifelse(is.na(.), "NA", .)))

    }

  nLoops <- nrow(parameter_grid)
  cTabsL <- list()

  # for loop and looping through etc!
  for (i in 1:nLoops){

    if (!missing(.subgroup) & !is.null(.subgroup)){

      cTabsL[[i]] <-
        crosstab_single(.dataset = .dataset,
                        .col_var = .col_var,
                        .subgroup = .subgroup,
                        .missing = .missing,
                        .row_var = parameter_grid$.row_var[[i]],
                        .subgroup_values = parameter_grid$.subgroup_values[[i]])


    } else {

      cTabsL[[i]] <-
        crosstab_single(.dataset = .dataset,
                        .col_var = .col_var,
                        .row_var = parameter_grid$.row_var[[i]],
                        .subgroup = .subgroup,
                        .subgroup_values = NA_character_,
                        .missing = .missing)

    }

  }


cTabsCombined <- tibble::as_tibble(dplyr::bind_rows(cTabsL))

# Adding labels if it exists!
if (!missing(.l_row_label)) {

  # if labels provided as a character vector
  if (is.character(.l_row_label)) {
    labeldb <- tibble::tibble(
      COLUMN_NAME = .l_row_var,
      COLUMN_LABEL = .l_row_label)
  } else {
    labeldb <- .l_row_label
  }

  cTabsCombined <- cTabsCombined %>%
    dplyr::left_join(labeldb, by = c("Variable" = "COLUMN_NAME")) %>%
    dplyr::rename(Variable_Label = "COLUMN_LABEL") %>%
    dplyr::relocate(Variable, Variable_Label)
}

# Outputting to an excel spreadsheet
if (!missing(.output_path) & !missing(.output_file)) {

  create_excel(
    .output_path = .output_path,
    .output_file = .output_file,
    .x = cTabsCombined
  )
}

return(cTabsCombined)
}

#' crosstab single
#'
#' @param .dataset data to supply
#'
#' @param .col_var the cross tab variable on the column end
#' @param .row_var variables to cross again on the X axis!
#' @param .subgroup 3 variable to enhance from 2 way to 3 way
#' @param .missing Logical, keep missing in data
#' @param .max_lvls maximum number of levels for a categorical variable
#' @param .row_label optional row label
#' @param .subgroup subgroup variable
#' @param .subgroup_values subgroup values
#'
#' @return tibble of single variable crosstab
#'
crosstab_single <- function(
    .dataset,
    .col_var,
    .row_var,
    .row_label = NULL,
    .subgroup = NULL,
    .subgroup_values = NULL,
    .missing = FALSE,
    .max_lvls = 24) {


  # cleaning the data set
  indata <- crosstab_cleaning(
    .dataset = .dataset,
    .col_var = .col_var,
    .row_var = .row_var,
    .missing = .missing,
    .subgroup = .subgroup,
    .subgroup_values = .subgroup_values
  )

  # creating crosstab
  if(.missing){

    crossData <- janitor::tabyl(indata,
                              !!sym(.col_var),
                              !!sym(.row_var),
                              show_na = TRUE,
                              show_missing_levels = TRUE)

  } else {

    crossData <- janitor::tabyl(indata,
                           !!sym(.col_var),
                           !!sym(.row_var),
                           show_na = FALSE,
                           show_missing_levels = FALSE)

  }
  # calculating chisquare estimation
  chisq_results <- crosstab_chisq(.dataset = crossData)
  chisq_output <- chisq_results$chisq_output
  ChiSqWarning <- chisq_results$ChiSqWarning

  # combing all results including percentage calculation into a single table
  crosstab_results <- crosstab_table(
    .dataset = crossData,
    .chisq_output = chisq_output,
    .chisq_warning = ChiSqWarning,
    .row_var = .row_var,
    .col_var = .col_var,
    .subgroup = .subgroup,
    .subgroup_values = .subgroup_values
  )

  return(crosstab_results)

}

#' Crosstab cleaning
#'
#' @param .dataset DataFrame that needs to be cleaned
#' @param .col_var ColVar to cross tab by
#' @param .row_var Row var to cross tab by
#' @param .subgroup 3 way table
#' @param .subgroup_values subgroup values
#' @param .missing Missing
#'
#' @return tibble with only necessary variables and rows
crosstab_cleaning <- function(.dataset,
                              .col_var,
                              .row_var,
                              .missing = FALSE,
                              .subgroup = NULL,
                              .subgroup_values = NULL){

  # we must do it for each of the
  if (.missing) {

      indata <- .dataset %>%
        dplyr::select({{.col_var}}, {{.row_var}}) %>%
        #select relevant vars from data
        # Coercing Actual NA values to the string "NA", making all character
        dplyr::mutate(dplyr::across(dplyr::everything(),
                                    ~ifelse(is.na(.), "NA", .)))


  } else {

      # If Missing set to FALSE, remove observations where there are any missing
      indata <- .dataset %>%
        dplyr::select({{.col_var}}, {{.row_var}}) %>%
        tidyr::drop_na()

  }

  # .subgroup - makes the two way table 3 ways!
  #  use of .subgroup.
  if (!missing(.subgroup) & !is.null(.subgroup)) {

    indata <- .dataset %>%
      dplyr::filter(!!sym(.subgroup) == {{.subgroup_values}})

  } else {

    indata <- .dataset

  }

  return(indata)

}

#' Crosstab Chisquare
#'
#' @param .dataset Input dataframe
#'
#' @return list with chi-square and any warnings
crosstab_chisq <- function(.dataset){

  # Capturing chisq message
  ChiSqWarning <- tryCatch(
    janitor::chisq.test(.dataset),
    error = function(e) {
      return(paste("Chi-squared failed to calculate.",
                   "Value is likely scalar."))
    },
    warning = function(w) {
      return("Chi-squared approximation may be incorrect")
    }
  )
  if (inherits(ChiSqWarning, "htest")) {
    ChiSqWarning <- "No error or warning."
  }
  # Creating chisq object
  if (ChiSqWarning == paste("Chi-squared failed to calculate.",
                            "Value is likely scalar.")) {
    # If chisq cannot be calculate, synthetically create list to ref
    chisq_output <- list(
      statistic = NA_real_,
      p.value = NA_real_
    )
  } else {
    chisq_output <- suppressWarnings(janitor::chisq.test(.dataset))
  }

  return(list(chisq_output = chisq_output,
              ChiSqWarning = ChiSqWarning))
}

#' Crosstab tabling
#'
#' @param .dataset input data
#' @param .chisq_warning Chisquare warning
#' @param .chisq_output Chisquare output
#' @param .row_var RowVar
#' @param .subgroup  ByVal
#' @param .row_label optional RowLabel
#' @param .subgroup_values subgroup values
#' @param .col_var ColVar
#'
#' @return tibble with
crosstab_table <- function(.dataset,
                           .chisq_output,
                           .chisq_warning,
                           .row_var,
                           .col_var,
                           .row_label,
                           .subgroup,
                           .subgroup_values){

  # ZH: Insert header
  crosstab_final <- .dataset %>%
    janitor::adorn_totals(c("row","col"),
                          na.rm = FALSE)
  crosstab_final <-  crosstab_final %>%
    dplyr::mutate(
      ChiSqValue = .chisq_output$statistic,
      ChiSqProb = .chisq_output$p.value
    ) %>% #add in chisq vals
    # extract freq and perc from tabyl in new cols, drop the old
    dplyr::mutate(
      dplyr::across(!c(1, "ChiSqValue", "ChiSqProb"), list(
        Frequency = ~as.numeric(stringr::str_extract(.x, "\\d+"))
      ), .names = "{.fn}_{.col}")) %>%
    dplyr::mutate(ChiSqWarning = .chisq_warning) %>%
    dplyr::select(1, ChiSqValue, ChiSqProb, ChiSqWarning,
                  dplyr::contains(c("ColPct", "Frequency"))) %>%
    # create cols specifying row vars and values
    dplyr::mutate(Variable = .row_var, .before = 1)%>%
    dplyr::rename_with(~return("Value"), dplyr::matches(.col_var))

  if (!missing(.subgroup) & !is.null(.subgroup)) {
    crosstab_final <- crosstab_final %>%
      dplyr::mutate(!!.subgroup  := {{.subgroup_values}},
                    UniqueID = paste0(Variable, "_", {{ .subgroup_values}}),
                    .after = 1)
  }

  return(crosstab_final)

}
