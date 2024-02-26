#### All Row Variable Testing
testthat::test_that("LM Extract Unit Testing",{

  # load ht_extract_sas output
  extract_sas <-
    readxl::read_excel(
      path = testthat::test_path("sas_tests", "lmglmcomparison.xlsx"))

  # Convert COHORT to a factor for analysis
  finalBif <- bif_final %>%
    dplyr::mutate(CohortFactor = factor(COHORT))

  # Create vector of dependent variables
  finalBif_dependents <- names(finalBif[setdiff(names(finalBif), c("SAMPLEID", "RA_CODE", "CohortFactor",
                                                                   "COHORT", "RA_DATE", "AGE",
                                                                   "DOB", "blcurhrs_73plus"))])
  finalBif_treatment  <- "RA_CODE"
  finalBif_covariates <- "CohortFactor"

  extract_r <- lm_extract(
    .dataset = finalBif,
    .dependents = finalBif_dependents,
    .treatment= finalBif_treatment,
    .covariates = finalBif_covariates,
    .confintalpha = .9)  # confidence interval - p value is .1 here.


  extract_sas$Dependent <- tolower(extract_sas$Dependent)
  colnames(extract_sas) <- tolower(colnames(extract_sas))

  extract_sas <- extract_sas %>%
    dplyr::select_if(~!any(is.na(.)))

  extract_sas <- extract_sas[extract_sas$dependent %in% finalBif_dependents, ]

  extract_r_matching <-
    extract_r %>%
    dplyr::select(
      dependent, estimate_p, estimate_c, impact, n_p, n_c, n_total, std_error,
      r_squared, unadj_impact, dep_var_mean, dep_var_sd
    ) %>%
    dplyr::filter(dependent != "sample_size")

  extract_sas_matching <-
    extract_sas %>%
    dplyr::select(
      dependent, adjmean_p, adjmean_c,
      impact_p_c, nobs_p, nobs_c, nobs,
      stderr_p_c, rsquare, unadj_impact_p_c, depvarmean,
      depvarstddev
    ) %>%
    dplyr::rename(
      dependent = dependent,
      estimate_p = adjmean_p,
      estimate_c = adjmean_c,
      impact = impact_p_c,
      n_p = nobs_p,
      n_c = nobs_c,
      n_total = nobs,
      std_error = stderr_p_c,
      r_squared = rsquare,
      unadj_impact = unadj_impact_p_c,
      dep_var_mean = depvarmean,
      dep_var_sd = depvarstddev
    )

  testthat::expect_equal(object = extract_r_matching,
                         expected = extract_sas_matching,
                         tolerance = 0.01)

})

