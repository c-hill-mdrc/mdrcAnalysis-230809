#' Applies flag to duplicate observations and generate summary tables
#' @param .data Dataframe of data to check for duplicates
#' @param .group_vars Variables to check for duplicates across
#' Defaults to all
#' @return Original Dataframe with "DUPFLAG" column attached if detected as
#' a duplicate
#' @export
check_duplicate <- function(.data, .group_vars){
  # Setting .group_vars to all if none specified
  if(missing(.group_vars)){
    .group_vars <- names(.data)
  }

  # Creating a dummy variable on .data indicating duplicate
  dup_check <- .data %>%
    dplyr::group_by(dplyr::across(tidyselect::all_of(.group_vars))) %>%
    dplyr::arrange(.by_group = T) %>%
    dplyr::mutate(DUPFLAG = dplyr::if_else(dplyr::n() == 1, 0, 1)) %>%
    dplyr::ungroup()

  # Printing table counting duplicate observations
  dup_check %>%
    dplyr::count(DUPFLAG, name = "Number") %>%
    # Making flag verbose
    dplyr::mutate(
      Condition = dplyr::if_else(
        DUPFLAG == 0,
        paste("Unique: Input observations with",
              "identifier on only one observation"),
        paste("Duplicates: Input observations with",
              "identifier on more than one observation")
      )
    ) %>%
    # Adding Total row count
    dplyr::add_row(.,  Condition = "Total Input Observations",
                   dplyr::summarise(., Number = sum(Number))) %>%
    dplyr::mutate(Percent = round((Number / max(Number)) * 100, 2)) %>%
    dplyr::arrange(dplyr::desc(Percent)) %>%
    dplyr::select(Condition, Number, Percent) %>%
    # # Keeping just as dataframe
    # flextable::flextable() %>%
    # flextable::set_caption("Count of duplicated observations") %>%
    print()

  # Printing table summarizing combinations
  .data %>%
    dplyr::group_by(dplyr::across(tidyselect::all_of(.group_vars))) %>%
    dplyr::summarise(
      Condition = dplyr::if_else(
        n() == 1,
        paste("Unique combination of group variables"),
        paste("Non-unique combiantion of group variables")
      ),
      Number = n(),
      .groups = "drop") %>%
    dplyr::add_row(., Condition = "Total combinations of group variables",
                   dplyr::summarise(., Number = sum(Number))) %>%
    dplyr::mutate(Percent = round((Number / max(Number)) * 100, 2)) %>%
    dplyr::arrange(dplyr::desc(Percent)) %>%
    # # Keeping just as dataframe
    # flextable::flextable() %>%
    # flextable::set_caption("Summary of combinations of group variables.") %>%
    print()

  # Return original dataframe but with duplicate dummy
  invisible(dup_check)
}
