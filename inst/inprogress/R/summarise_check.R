#' Generates dataframe comparing raw data and calculaated measures
#'
#' Takes in a "wide" dataset containing measure and component raw data as
#' columns. Calculates a "check" by pivoting the wide dataset and running
#' \code{dplyr::summarise} using a purrr style ~ function
#'
#' @param .data Dataset containing both the measure and the raw data used to
#' generate the measure. Must be a wide dataset, i.e. each variable specified
#' is a separate column
#' @param .f \code{purrr} style ~ function the describes how raw data should be
#' transformed into the measure
#' @param .measure Single variable to check against. Supports tidyselection
#' @param ... Names of relevant raw data used to generate measure. Supports all
#' tidyselection verbs
#' @param .nonmising Maximal number of cases to display where no data is
#' missing and measure check matches with measure variable
#' @param .na_drop Boolean specifying if cases where all specified variables
#' are missing should be ignored
#' @param .tolerance Maximal difference allowed between checked measure and
#' actual measure. used to account for rounding/datatype discrepancies. Set to
#' \code{"exact"} no tolerance desired
#' especially when \code{.data} originates from a non-R data source
#' @returns Dataframe of original raw data, original measure, checked measure,
#' and descriptives such as the number of occurences of that variable
#' combination and number of missing raw data variables. Displays all variable
#' combinations where discrepancy found, where some raw data is missing, and
#' up to \code{.nonmissing} other measure/variable combinations
#' @export
summarise_check <- function(.data,
                            .f,
                            .measure,
                            ...,
                            .nonmissing = 10,
                            .na_drop = TRUE,
                            .tolerance = 0.00001) {

  # Selecting relevant variables
  # # Defusing selection expression
  measure_expr <- rlang::enquo(.measure)
  vars_expr <- rlang::expr(c(...))
  # # Identifying position of vars
  measure_pos <- tidyselect::eval_select(measure_expr, data = .data)
  vars_pos <- tidyselect::eval_select(vars_expr, data = .data)
  # # Creating new, selected dataframe
  data <- rlang::set_names(.data[c(measure_pos, vars_pos)],
                           c(names(measure_pos), names(vars_pos)))
  # if .na_drop is true, filter out where variable missing
  if (.na_drop == T) {
    data <- tidyr::drop_na(data)
  }

  # Generating unique, per-row id for data
  data <- dplyr::mutate(data, ID = dplyr::row_number())
  # Calculating Measure according to .f specified
  # note calculation will be done on on a long format
  func <- purrr::as_mapper(.f)
  check_data <- data[c("ID", names(vars_pos))] %>%
    tidyr::pivot_longer(
      cols = names(vars_pos),
      names_to = "variable_name",
      values_to = "NumVal"
    ) %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(
      !!sym(paste0(names(measure_pos), "_check")) := func(NumVal,
                                                          variable_name),
      "N Missing" = sum(is.na(NumVal)))
  new_names <- c(names(check_data),
                 purrr::discard(names(data), ~. == "ID"))

  ordered_data <- dplyr::full_join(data, check_data, by = "ID") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(new_names))) %>%
    dplyr::summarise(Frequency = n(), .groups = "drop") %>%
    dplyr::arrange(Frequency)

  # Order and filter displayed cases
  check_dataframe <- dplyr::bind_rows(
    # Functionalized if else to use in bind_rows
    (if (.tolerance == "exact") {
      dplyr::filter(
        ordered_data,
        !!rlang::sym(names(measure_pos)) != !!rlang::sym(
          paste0(names(measure_pos), "_check"))
      )} else {
        dplyr::filter(
          ordered_data,
          abs(magrttr::subtract(
            !!rlang::sym(names(measure_pos)),
            !!rlang::sym(paste0(names(measure_pos), "_check")))) > .tolerance
        )}),
    # Display cases with missing data
    dplyr::filter(ordered_data, `N Missing` > 0),
    # Display at most .nonmissing other cases
    dplyr::filter(
      ordered_data,
      !!rlang::sym(names(measure_pos)) == !!rlang::sym(paste0(
        names(measure_pos), "_check"))) %>%
      dplyr::filter(`N Missing` == 0) %>%
      dplyr::filter(dplyr::row_number() <= .nonmissing)) %>%
    dplyr::distinct()
  return(check_dataframe)
}
