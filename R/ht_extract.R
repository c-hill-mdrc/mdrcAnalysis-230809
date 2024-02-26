#' Homogeneity Test Extract
#'
#' The HT statistic is a test of homogeneity of effect size.
#' The HT statistic is used to check whether impacts for subgroups are
#' statistically different.
#'
#' @param .dataset Input data
#' @param .subgroup String specifying column encoding subgroup rows belong to
#' @param .dependents Character vector of dependent variables to run models on
#' @param .treatment String specifying treatment variable
#' @param .covariates Character Vector of covariates to run model on
#' @param .asfactorcovs Character vector of covariates to be run as
#'                      factors/categorical variables
#' @param .confintalpha decimal value of desired alpha level for
#'                      confidence intervals. Default is 0.95
#'
#' @param .inc_trail   RTU parameter to include trail in output
#' @param .inc_sample  RTU parameter to include sample size row in output
#'
#' @param .output_path file path to output location
#' @param .output_file output file name
#'
#' @return Tibble displaying model outputs over all combinations of
#' dependents and subgroup values
#'
#'
#' @export
ht_extract <- function(.dataset,
                       .subgroup,
                       .dependents,
                       .treatment,
                       .covariates,
                       .asfactorcovs = NA_character_,
                       .confintalpha = 0.95,
                       .inc_trail = TRUE,
                       .inc_sample = TRUE,
                       .output_path = NA_character_,
                       .output_file = NA_character_) {

  # Parameter checks
  if (!is.data.frame(.dataset)) {
    stop("ERROR: lm_extract expects a data frame as the first argument.")
  }
  if (!rlang::is_scalar_character(.subgroup)) {
    stop("ERROR: ht_extract expects a string")
  }
  if (!is.character(.dependents)) {
    stop("ERROR: lm_extract expects a character vector of dependents.")}
  if (!is.character(.covariates)) {
    stop("ERROR: lm_extract expects a character vector of covariates.")}
  if (!is.character(.asfactorcovs)) {
    stop("ERROR: lm_extract expects a character vector of factor covariates")}
  if (!is.character(.treatment)) {
    stop("ERROR: lm_extract expects a character vector of treatment.")}
  if (dplyr::is.grouped_df(.data)) {
    # Ungrouping if necessary
    .data <- dplyr::ungroup(.data)
  }

  # estimate treatment effects for each dependent variable by subgroup

  # Generating a dataframe by row for each unique combination of
  # .dependents and .subgroup. Since we are using the lm_extract custom f
  # unction, cannot perform this on a grouped dataframe

  # Generating all combos of .dependents and dataframes filtered by subgroup
  # subgroup splits here (add explicit in splits!)
  parameter_combinations <- tidyr::expand_grid(
    ".dependent" = .dependents,
    ".data" = split(.dataset,
                    dplyr::pull(.dataset, {{.subgroup}}))
  )

  #browser()

  # Loop over each combination of parameters and run lm_extract
  lm_mdls <- parameter_combinations %>%
    purrr::pmap_dfr(function(.dependent, .data){

      SubGrpDesc <- .data %>%
        dplyr::select(SubGrpVal = !!rlang::sym(.subgroup)) %>%
        dplyr::distinct() %>%
        dplyr::mutate(SubGrpVar = !!rlang::quo(.subgroup))
      # Attaching results of lm_extract
      SubGrpDesc <- SubGrpDesc %>%
        dplyr::bind_cols(
          lm_extract(.data,
                     .dependent,
                     .treatment,
                     .covariates,
                     .asfactorcovs,
                     .confintalpha,
                     .inc_trail = FALSE,
                     .inc_sample = FALSE,
                     .output_path = NA_character_,
                     .output_file = NA_character_)
        )

      # remove the additional sample line
      SubGrpDesc <-
        SubGrpDesc %>%
          dplyr::filter(dependent != "sample_size")


      # Diagnostic if any std errors are equal to 0
      # ZH: should the standard errors be small or large? (do not expect a significant impact so std errors would be large compared to when there is an impact!)
      if(SubGrpDesc$std_error == 0){
        warning(paste("Standard error of 0 detected where treatment is ",
                      .treatment,
                      ", dependent is ",
                      .dependent,
                      " and subgroup is ", .data[[1, .subgroup]]))
      }
      return(SubGrpDesc)
    })

  # calculate the sum of impacts (expressed in squared se) divided by
  # the sum of groups (expressed in squared se terms) ~ a weighted mean
  # impact across groups
  # ZH: if you can send me the paper on htextract underlying subgroup paper!
  theta_den <- lm_mdls %>%
    dplyr::mutate(isesq = (std_error^2),
                  numtheta = impact/isesq,
                  dentheta = 1/isesq) %>%
    dplyr::group_by(dependent) %>%
    dplyr::summarise(Num_Den = sum(numtheta)/sum(dentheta),
                     .groups = "drop")

  # calculate the homogeneity of Impacts as
  # the sum of each groups deviation from the weighted mean impact across groups
  # (group impact minus the weighted mean impact across groups in
  # squared standard errors) divided by the squared standard error
  HTtest <- lm_mdls %>%
    dplyr::left_join(theta_den, by = "dependent") %>%
    dplyr::mutate(HTi = ((impact-Num_Den)^2) / (std_error^2)) %>%
    dplyr::group_by(dependent) %>%
    dplyr::summarise(HT=sum(HTi),
                     .groups = "drop")

  # c01alculate the probability of the homogeneity of impacts to be random
  # p value
  ProbHTtest <- lm_mdls %>%
    dplyr::left_join(HTtest, by = "dependent") %>%
    dplyr::group_by(dependent,HT) %>%
    dplyr::summarise(ProbHT = 1 - stats::pchisq(unique(HT),dplyr::n()-1),
                     .groups = "drop") %>%
    dplyr::mutate(HT_Stars = create_stars(ProbHT))

  # stack the HT statistics onto the subgroup statistics and reorder the results
  lm_mdls <- lm_mdls %>%
    dplyr::mutate(UniqueId = paste0(
      dependent,"_", SubGrpVar, "_", SubGrpVal )) %>%
    dplyr::select(-SubGrpVal , -SubGrpVar ) %>%
    dplyr::left_join(ProbHTtest, by = "dependent") %>%
    dplyr::relocate(UniqueId,
                    dependent,
                    starts_with("estimate"),
                    impact,
                    p_value,
                    stars,
                    std_error,
                    HT,
                    ProbHT,
                    HT_Stars)

  # pull out the sample counts
  sampleCounts <- sample_size(
    .dataset = .dataset,
    .treatment = .treatment
  )

  # Adding the sample size function to the larger table
  if(.inc_sample == TRUE) {

    names(sampleCounts)[1] <- names(lm_mdls)[1]

    sampleCounts <-
      sampleCounts %>%
        dplyr::select(-n_na)

    lm_mdls <-
      lm_mdls %>%
        dplyr::bind_rows(sampleCounts)

  } else {

    lm_mdls

  }

  # Adding the trail function
  if(.inc_trail == TRUE) {

    ## Capturing the function parameters
    call <- match.call()

    lm_mdls <-
      lm_mdls %>%
      dplyr::bind_cols(trail(call))

  } else {

    lm_mdls

  }

  # Outputting to an excel spreadsheet
  if (!is.na(.output_path) & !is.na(.output_file)) {

    create_excel(
      .output_path = .output_path,
      .output_file = .output_file,
      .x = lm_mdls
    )
  }

  return(lm_mdls)
}
