#' crosstab_extract main function
#'
#' For creating crosstabs for an entire data frame. Options include creating
#' crosstabs by subgroup, including missing values, and including chi-square
#' test results.
#'
#' @param .dataset required, data frame with .col_var and .row_var
#' @param .col_var required, variable in columns of cross-tabulation
#' @param .row_vars optional, character vector of variables in rows of cross-
#'  tabulation. All variables are included if not provided.
#' @param .subgroup optional, character vector of variables to group the data by
#' @param .chisq optional, Boolean, include chi-square test, default is FALSE
#' @param .missing optional, Boolean, include missing values, default is FALSE
#' @param .max_lvls optional, maximum number of levels for a variable to be
#'  included in the cross-tabulations. Default is 24.
#' @param .row_labels optional, data frame with variable labels. Data frame
#'  should
#' @param .output_path optional, path to output location
#' @param .output_file optional, filename for output file
#'
#' @returns A tibble of all the cross-tabulations
#'
#' @export crosstab_extract
#'
#' @examples
#' crosstab_extract(
#'   .dataset  = bif_final,
#'   .col_var  = "COHORT",
#'   .row_vars = c("blmale", "bldiplomas_as")
#'   )
#'
crosstab_extract <-
  function(.dataset
          ,.col_var
          ,.row_vars    = NA_character_
          ,.subgroup    = NA_character_
          ,.chisq       = FALSE
          ,.missing     = FALSE
          ,.max_lvls    = 24
          ,.row_labels  = NA_character_
          ,.output_path = NA_character_
          ,.output_file = NA_character_
          )
    {

    # drop unnecessary variables

    ## if no .row_vars were provided, include all variables
    if (missing(.row_vars)) {
      .row_vars <- setdiff(names(.dataset)
                          ,c(.col_var, .subgroup)
                           )
    } # end if no .row_vars provided

    ## check all .row_vars for .max_levels
    ### store all requested row variables
    old_rvs <- .row_vars

    ### create list of row variables that are less than max_lvls
    .row_vars <-
      .dataset |>
      dplyr::select(tidyselect::all_of({{.row_vars}})) |>
      dplyr::select(
        tidyselect::where(
          ~dplyr::n_distinct(., na.rm = !.missing) <= .max_lvls)) |>
      names()

    ### List any removed variables and report
    rmvd_rvs <- setdiff(old_rvs, .row_vars)

    if (length(rmvd_rvs) > 0) {
      print(paste("The following requested variables exceed the maximum number"
                  ,"of levels and were removed from analysis:"
                  ,paste(rmvd_rvs, collapse = ", ")))
    }

    ## create vector of necessary variables
    cols <- c(.col_var, .row_vars, .subgroup)

    ## create vector of necessary variables available in input data frame
    keep <- intersect(cols, names(.dataset))

    ## create data frame with only necessary variables
    keepdf <-
      .dataset |>
      dplyr::select(tidyselect::all_of(keep))

    # run helpers on each combination of .subgroup and .row_var
    if (missing(.subgroup)) {

      # no subgroups, iterate over row variables
      outdf <-
        purrr::map(.row_vars, \(x) crosstab_single(.dataset = keepdf
                                                  ,.col_var = {{.col_var}}
                                                  ,.row_var = x
                                                  ,.missing = {{.missing}}
                                                  ,.chisq   = {{.chisq}}
                                                  )
                   ) |>
        purrr::list_rbind()

    } else {

      # subgroups, iterate over groups and row variables

      # first create split data frames
      groupdf <-
        keepdf |>
        dplyr::group_by(dplyr::pick(.subgroup)) |>
        dplyr::glimpse()

      # iterate over row variables and subgroups
      outdf <-
        purrr::map(.row_vars, \(y) {
          dplyr::group_map(groupdf, ~
            crosstab_single(.dataset = .x
                           ,.col_var = {{.col_var}}
                           ,.row_var = y
                           ,.missing = {{.missing}}
                           ,.chisq   = {{.chisq}}
            )
          ) |>
            purrr::list_rbind(names_to = "Subgroup") |>
            dplyr::mutate(Subgroup = paste0(unique(.subgroup)
                                          ,"_"
                                          ,unique(.dataset[[.subgroup]])[Subgroup])
            )
        }) |>
        purrr::list_rbind()

    } # end looping over subgroups and row variables

    # add row variable labels if provided
    if (!is.na(.row_labels)) {
      outdf <-
        dplyr::left_join(outdf, {{.row_labels}}, by = "Variable") |>
        dplyr::relocate(Label, .after = Variable)
    } # end adding row labels

    # add trail
    outdf <-
      outdf |>
      dplyr::bind_cols(trail(call = match.call.defaults()))

    # Outputting to an excel spreadsheet
    if (!missing(.output_path) & !missing(.output_file)) {
      create_excel(.output_path = .output_path,
                   .output_file = .output_file,
                   .x = outdf
                   )
      } # end exporting to Excel if provided

  return(outdf)
}

#' crosstab_extract helper function: runs crosstab_clean and crosstab_single
#'
#' Consolidates both steps for use in main function
#'
#' @param .dataset required, data frame with .col_var and .row_var
#' @param .col_var required, variable in columns of cross-tabulation
#' @param .row_var required, variable in rows of cross-tabulation
#' @param .missing optional, Boolean, include missing values, default is FALSE
#' @param .chisq optional, Boolean, request chi-square test, default is FALSE
#'
#' @return clean tibble with all requested statistics
#'
crosstab_single <-
  function(.dataset
           ,.col_var
           ,.row_var
           ,.missing  = FALSE
           ,.chisq    = FALSE
  ){

    # Run both functions
    singledf <-
      .dataset |>
      crosstab_clean(.col_var = {{.col_var}}
                    ,.row_var = {{.row_var}}
                    ,.missing = {{.missing}}
                    ) |>
      crosstab_stats(.col_var = {{.col_var}}
                    ,.row_var = {{.row_var}}
                    ,.chisq   = {{.chisq}}
                    )

     return(singledf)
  }

#' crosstab_extract helper function to create cross-tab data frame
#'
#' Takes input data frame, subsets for relevant .col_var and .row_var. Checks
#' for missing values and converts to string, if requested. Outputs data
#' frame with minimally necessary data for cross-tab.
#'
#' @param .dataset required, data frame with .col_var and .row_var
#' @param .col_var required, variable in columns of cross-tabulation
#' @param .row_var required, variable in rows of cross-tabulation
#' @param .missing optional, Boolean, include missing values, default is FALSE
#'
#' @return tibble with necessary variables and rows
#'
crosstab_clean <-
  function(.dataset
           ,.col_var
           ,.row_var
           ,.missing  = FALSE
           ){

    # Remove missing values, if requested
    if (.missing) {
      # if .missing, convert all variables to character and NA values to "NA"
      cleandf <-
        .dataset |>
        dplyr::select(tidyselect::all_of({{.col_var}})
                     ,tidyselect::all_of({{.row_var}})) |>
        # making variables character, coercing NA values to "NA"
        dplyr::mutate(dplyr::across(dplyr::everything(),
                                    ~ifelse(is.na(.), "NA", .)))
      } else {
        # if .missing = FALSE, remove observations with NA values
        cleandf <-
          .dataset |>
          dplyr::select(tidyselect::all_of({{.col_var}})
                       ,tidyselect::all_of({{.row_var}})) |>
          tidyr::drop_na()
      }

    # Output data frame
    return(cleandf)
}


#' crosstab_extract helper which runs the crosstab for a pair of variables
#'
#' Takes input data frame, creates cross-tabulation, calculates desired
#' frequencies and percentages. Also runs chi-square test, if requested.
#'
#' @param .dataset required, data frame with .col_var and .row_var
#' @param .col_var required, variable in columns of cross-tabulation
#' @param .row_var required, variable in rows of cross-tabulation
#' @param .chisq optional, Boolean, request chi-square test, default is FALSE
#'
#' @return tibble of single variable cross-tab with chi-square, if requested
#'
crosstab_stats <-
  function(.dataset
           ,.col_var
           ,.row_var
           ,.chisq     = FALSE
           ) {

    # Calculate frequences and percentages
    crosstab <-
      .dataset |>
      dplyr::count(dplyr::pick(.col_var)
                  ,dplyr::pick(tidyselect::all_of(.row_var))) |>
      dplyr::mutate(UniqueID = paste0({{.row_var}},
                               "_",
                               .data[[{{.row_var}}]]
                               )) |>
      dplyr::relocate(UniqueID,
                      .before = tidyselect::everything()) |>
      dplyr::group_by(dplyr::pick({{.col_var}})) |>
      dplyr::mutate(CumFreq = cumsum(n),
                    Pct     = n/max(CumFreq)*100,
                    CumPct  = cumsum(Pct)) |>
      dplyr::ungroup() |>
      tidyr::pivot_longer(cols = {{.row_var}},
                          names_to = "Variable") |>
      dplyr::mutate(Value = as.character(.data[["value"]]), .keep = "unused") |>
      dplyr::rename(Freq = n) |>
      tidyr::pivot_wider(names_from = tidyselect::all_of(.col_var)
                        ,values_from = c(Freq, CumFreq, Pct, CumPct)
                        ,values_fill = 0
                        )

    # If requested, add chisq
    if (.chisq) {
      # checking if .col_var and .row_var have at least 2 levels each
      clvls <- dplyr::n_distinct(.dataset[[.col_var]])
      rlvls <- dplyr::n_distinct(.dataset[[.row_var]])

      if (clvls >= 2 & rlvls >= 2) {
        # calculating chisquare estimation
        chisqr <-
          broom::tidy(
            suppressWarnings(
              stats::chisq.test(x = dplyr::pull(.dataset, var = {{.col_var}}),
                                y = dplyr::pull(.dataset, var = {{.row_var}}))
            )
          ) |>
          dplyr::rename(ChiSqStatistic = statistic)
      } else {
        print(paste("Variables"
                   ,.col_var
                   ,"or"
                   ,.row_var
                   ,"have fewer than 2 levels. Cannot calculate chi-square."))

        # create empty tibble
        chisqr <-
          tibble::tibble(ChiSqStatistic = NA
                        ,p.value        = NA
                        ,parameter      = NA
                        ,method         = NA_character_
                        )
      }

      # Combine results
      crosstab <-
        dplyr::bind_cols(crosstab, chisqr)
    }

    return(crosstab)
}

