#' mean extract function
#'
#' @param .dataset required, the data set with necessary variables
#' @param .subgroup optional, character vector of variables to group the data by
#' @param .quantile optional, request quantiles, default is FALSE
#' @param .inc_sample RTU parameter for sample size row inclusion, default is TRUE
#' @param .inc_trail  RTU parameter for trail inclusion, default is TRUE
#'
#' @return tibble with summary statistics for all variables
#' @export
#'
#' @examples
#' mean_extract(
#'             .dataset = sim_data_robust_reg,
#'             .subgroup = "treatment",
#'             .quantile = TRUE
#'             )
#'
#'
mean_extract <- function(
    .dataset,
    .subgroup = NA_character_,
    .quantile = FALSE,
    .inc_sample = TRUE,
    .inc_trail = TRUE
){


  # check if the subgroup is not empty
  if (!is.na(.subgroup)){

    mean_final <-
      .dataset %>%
        dplyr::group_split(!!rlang::sym(.subgroup), .keep = TRUE) %>%
        purrr::map(~ mean_extract_helper(.x,
                                         .quantile = .quantile,
                                         .subgroup = .subgroup)) %>%
        dplyr::bind_rows()


      if (.inc_trail == TRUE) {

        ## Capturing the function parameters
        call <- match.call()

        mean_final <-
          mean_final %>%
            dplyr::relocate(Variable, !!rlang::sym(.subgroup)) %>%
            dplyr::bind_cols(trail(call))

      } else {

        mean_final <-
          mean_final %>%
          dplyr::relocate(Variable, !!rlang::sym(.subgroup))

      }


    } else {

      if (.inc_trail == TRUE) {

        ## Capturing the function parameters
        funcParams <- match.call()

        # return mean calculation without the subgroups
        mean_final <- mean_extract_helper(
          .dataset = .dataset,
          .quantile = .quantile
        )

        ## Capturing the function parameters
        mean_final <- mean_final %>%
          dplyr::bind_cols(trail(funcParams))

      } else {

        # return mean calculation without the subgroups
        mean_final <- mean_extract_helper(
          .dataset = .dataset,
          .quantile = .quantile
        )

      }
    }

  # pull out the sample counts
  sampleCounts <-
    sample_size(.dataset   = .dataset
               ,.treatment = .subgroup
    )

  # rename the sample count function first column to the first column of the Tib
  names(sampleCounts)[1] <- names(mean_final)[1]

  # Adding the sample size function to the larger table
  if(.inc_sample == TRUE) {

    mean_final <-
      mean_final %>%
      dplyr::bind_rows(sampleCounts)

  }

  return(mean_final)

  } # mean_extract

mean_extract_helper <- function(
  .dataset,
  .quantile = FALSE,
  .subgroup = NA_character_
){

  # pull out all the variables
  .allvariables <- names(.dataset)

  # estimate the mean, standard deviation and max for all the dependent variables
  tmp_msdminmax <-
      .dataset %>%
      dplyr::select_if(is.numeric) %>%
      dplyr::summarise(dplyr::across(all_of(names(.)), list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE),
                                                            min = ~min(.x, na.rm = TRUE), max = ~max(.x, na.rm = TRUE)),
                                     .names ="{col}.{fn}")) %>%
      tidyr::pivot_longer(everything(), names_to = c("Variable", "function"), names_sep = "\\.") %>%
      tidyr::pivot_wider(id_cols = "Variable", names_from = "function") %>%
      dplyr::rename(Mean = mean,
                    Stddev = sd,
                    Min = min,
                    Max = max)

  # estimate the N values
  tmp_n <-
    .dataset %>%
    dplyr::select_if(is.numeric) %>%
    dplyr::summarise(dplyr::across(all_of(names(.)), list(n = ~dplyr::n()), .names = "{col}.{fn}")) %>%
    tidyr::pivot_longer(everything(), names_to = c("Variable", "function"), names_sep = "\\.") %>%
    tidyr::pivot_wider(id_cols = "Variable", names_from = "function") %>%
    dplyr::rename(n_total = n)

  # create missing measure
  tmp_na <-
        .dataset %>%
      dplyr::select_if(is.numeric) %>%
      dplyr::summarise(dplyr::across(all_of(names(.)), list(missing = ~sum(is.na(.))), .names = "{col}.{fn}")) %>%
      tidyr::pivot_longer(everything(), names_to = c("Variable", "function"), names_sep = "\\.") %>%
      tidyr::pivot_wider(id_cols = "Variable", names_from = "function") %>%
      dplyr::rename(NMiss = missing)

  # mean estimate final joins
  mean_final <- tmp_n %>%
    dplyr::inner_join(tmp_msdminmax, by = c("Variable")) %>%
    dplyr::inner_join(tmp_na, by = c("Variable"))


  # when quantiles are being called upon
  if(.quantile == TRUE) {

    tmp_quantiles <-
      .dataset %>%
        dplyr::select_if(is.numeric) %>%
        dplyr::summarise(dplyr::across(all_of(names(.)), list(q25 = ~quantile(., probs = 0.25),
                                                              q50 = ~quantile(., probs = 0.5),
                                                              q75 = ~quantile(., probs = 0.75)),
                                       .names = "{col}.{fn}"
        )) %>%
        tidyr::pivot_longer(everything(), names_to = c("Variable", "function"), names_sep = "\\.") %>%
        tidyr::pivot_wider(id_cols = "Variable", names_from = "function") %>%
        dplyr::rename(Quantile25 = q25,
                      Quantile50 = q50,
                      Quantile75 = q75)

    # Update the mean final
    mean_final <- mean_final %>%
      dplyr::inner_join(tmp_quantiles, by = c("Variable"))

  }

  # if subgroup is part of the function
  if(!is.na(.subgroup)){

    mean_final <- mean_final %>%
      dplyr::mutate({{.subgroup}} := unique(.dataset[[.subgroup]]))

    mean_final <- mean_final %>%
      dplyr::filter(!Variable %in% .subgroup)

  } #is.na(.subgroup)

  return(mean_final)

} # mean_extract_helper



