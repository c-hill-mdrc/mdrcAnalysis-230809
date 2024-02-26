# Checking for corr_extract functions

test_that("Error if non-numeric column passed to .vars", {
  expect_error(corr_extract(sim_data, "individuals"))
})

test_that("Warnings occur if items reversed for alphas need to be reversed", {
  # Add dplyr
  withr::local_package("dplyr")

  suppressWarnings(expect_warning(corr_extract(sim_data, .warnings = TRUE)))
})

test_that(".stats defaults to mean only", {
  # Add purrr
  withr::local_package("purrr")

  expect_identical(
    purrr::keep(names(corr_extract(sim_data)), ~. %in% eval(formals(corr_extract)$.stats)),
    "mean"
  )
})

test_that("Correlation missing if sd is 0", {
  # Add dplyr
  withr::local_package("dplyr")

  expect_true(all(sim_data$pre_income_raw_25pct == sim_data$pre_income_raw_25pct[1]))
  expect_true(all(is.na(
    corr_extract(sim_data,
                 .with_var = names(dplyr::select(sim_data,
                                                 tidyselect:::where(is.numeric),
                                                 -pre_income_raw_25pct)),
                 .use = "everything",
                 .prob = FALSE,
                 .nobs = FALSE,
                 .alpha = FALSE)$pre_income_raw_25pct))
  )
  expect_true(all(is.na(
    corr_extract(sim_data,
                 .with_var = names(dplyr::select(sim_data,
                                                 tidyselect:::where(is.numeric),
                                                 -pre_income_raw_25pct)),
                 .use = "all",
                 .prob = FALSE,
                 .nobs = FALSE,
                 .alpha = FALSE)$pre_income_raw_25pct))
  )
  expect_true(all(is.na(
    corr_extract(sim_data,
                 .with_var = names(dplyr::select(sim_data,
                                                 tidyselect:::where(is.numeric),
                                                 -pre_income_raw_25pct)),
                 .use = "pairwise",
                 .prob = FALSE,
                 .nobs = FALSE,
                 .alpha = FALSE)$pre_income_raw_25pct))
  )
  expect_true(all(is.na(
    corr_extract(sim_data,
                 .with_var = names(dplyr::select(sim_data,
                                                 tidyselect:::where(is.numeric),
                                                 -pre_income_raw_25pct)),
                 .use = "na",
                 .prob = FALSE,
                 .nobs = FALSE,
                 .alpha = FALSE)$pre_income_raw_25pct))
  )
})

test_that("Extracts identical between .use if .impute is TRUE", {
  # Add dplyr
  withr::local_package("dplyr")

  expect_identical(corr_extract(sim_data,
                                .use = "pairwise",
                                .impute = TRUE),
                   corr_extract(sim_data,
                                .use = "all",
                                .impute = TRUE))
  expect_identical(corr_extract(sim_data,
                                .use = "pairwise",
                                .impute = TRUE),
                   corr_extract(sim_data,
                                .use = "everything",
                                .impute = TRUE))
  expect_identical(corr_extract(sim_data,
                                .use = "pairwise",
                                .impute = TRUE),
                   corr_extract(sim_data,
                                .use = "complete",
                                .impute = TRUE))
  expect_identical(corr_extract(sim_data,
                                .use = "pairwise",
                                .impute = TRUE),
                   corr_extract(sim_data,
                                .use = "na",
                                .impute = TRUE))
})

test_that("pearson calculated properly", {
  # Add dplyr
  withr::local_package("dplyr")

  expect_equal(
    corr_extract(sim_data,
                 .vars = "cities",
                 .with_var = "population",
                 .method = "pearson",
                 .use = "pairwise",
                 .prob = FALSE,
                 .nobs = FALSE,
                 .alpha = FALSE,
                 .round = NULL)$cities,
    stats::cor(sim_data$cities, sim_data$population,
                    use = "pairwise", method = "pearson"))
})

test_that("spearman calculated properly", {
  # Add dplyr
  withr::local_package("dplyr")

  expect_equal(
    corr_extract(sim_data,
                 .vars = "cities",
                 .with_var = "population",
                 .method = "spearman",
                 .use = "pairwise",
                 .prob = FALSE,
                 .nobs = FALSE,
                 .alpha = FALSE,
                 .round = NULL)$cities,
    stats::cor(sim_data$cities, sim_data$population,
                    use = "pairwise", method = "spearman"))
})

test_that("kendall calculated properly", {
  # Add dplyr
  withr::local_package("dplyr")

  expect_equal(
    corr_extract(sim_data,
                 .vars = "cities",
                 .with_var = "population",
                 .method = "kendall",
                 .use = "pairwise",
                 .prob = FALSE,
                 .nobs = FALSE,
                 .alpha = FALSE,
                 .round = NULL)$cities,
    stats::cor(sim_data$cities, sim_data$population,
                    use = "pairwise", method = "kendall"))
})

#### All Row Variable Testing
testthat::test_that("Correlation for Cities and Individuals Part 1",{
  # Add readxl
  withr::local_package("readxl")

  # load data set
  data("sim_data_robust_reg")

  # load ht_extract_sas output
  corr_extract_sas <-
    readxl::read_excel(
      path = testthat::test_path("sas_tests", "corrextractcomp1.xlsx"),
      sheet = "ALL")

  # clean up the column names in correlation
  corr_extract_sas <-
    corr_extract_sas %>%
      janitor::clean_names()

  # Convert COHORT to a factor for analysis
  finalBif <- bif_final %>%
    dplyr::mutate(CohortFactor = factor(COHORT))

  # converting to numeric from character to run correlation
  sim_data$individuals <- as.numeric(sim_data$individuals)
  sim_data_robust_reg$individuals <- as.numeric(sim_data_robust_reg$individuals)

  # corr extract r
  corr_extract_r <-
    corr_extract(sim_data_robust_reg, .vars = c("cities", "individuals"))

  # rename corr_extract_sas to align with
  corr_extract_sas <-
    corr_extract_sas %>%
      dplyr::rename(
        "Variable" = "variable",
        "Type" = "type",
        "Total_Raw_Alpha" = "raw_alpha_with_deleted_variable",
        "Std_Alpha" = "std_corr_with_total",
        "Total_Std_Alpha" = "std_alpha_with_deleted_variable",
        "mean" = "mean",
        "cities" = "cities",
        "individuals" = "individuals"
      ) %>%
      dplyr::select(
        Variable,
        Type,
        Std_Alpha,
        Total_Raw_Alpha,
        Total_Std_Alpha,
        mean,
        cities,
        individuals)

  # dropping one column from corr_extract_r (Raw_Alpha) to align
  corr_extract_r <-
    corr_extract_r %>%
      dplyr::select(-Raw_Alpha)

  # Restricting corr_extract to 6 columns
  corr_extract_r_first_part <-
    corr_extract_r %>%
      dplyr::select(Variable, Type, Std_Alpha,
                    Total_Raw_Alpha, Total_Std_Alpha,
                    mean)

  # Restricting corr_extract to 6 columns
  corr_extract_sas_first_part <-
    corr_extract_sas %>%
    dplyr::select(Variable, Type, Std_Alpha,
                  Total_Raw_Alpha, Total_Std_Alpha,
                  mean)


  # comparing ttest for no group by object
  testthat::expect_equal(object = corr_extract_sas_first_part,
                         expected = corr_extract_r_first_part,
                         tolerance = 0.01,
                         ignore_attr = TRUE)

})


#### All Row Variable Testing
testthat::test_that("Correlation for Cities and Individuals Part 2",{
  # Add readxl
  withr::local_package("readxl")

  # load ht_extract_sas output
  corr_extract_sas <-
    readxl::read_excel(
      path = testthat::test_path("sas_tests", "corrextractcomp1.xlsx"),
      sheet = "ALL")

  # clean up the column names in correlation
  corr_extract_sas <-
    corr_extract_sas %>%
    janitor::clean_names()

  # Convert COHORT to a factor for analysis
  finalBif <- bif_final %>%
    dplyr::mutate(CohortFactor = factor(COHORT))

  # converting to numeric from character to run correlation
  sim_data$individuals <- as.numeric(sim_data$individuals)
  sim_data_robust_reg$individuals <- as.numeric(sim_data_robust_reg$individuals)

  # corr extract r
  corr_extract_r <-
    corr_extract(sim_data_robust_reg, .vars = c("cities", "individuals"))

  # rename corr_extract_sas to align with
  corr_extract_sas <-
    corr_extract_sas %>%
    dplyr::rename(
      "Variable" = "variable",
      "Type" = "type",
      "Total_Raw_Alpha" = "raw_alpha_with_deleted_variable",
      "Std_Alpha" = "std_corr_with_total",
      "Total_Std_Alpha" = "std_alpha_with_deleted_variable",
      "mean" = "mean",
      "cities" = "cities",
      "individuals" = "individuals"
    ) %>%
    dplyr::select(
      Variable,
      Type,
      Std_Alpha,
      Total_Raw_Alpha,
      Total_Std_Alpha,
      mean,
      cities,
      individuals)

  # dropping one column from corr_extract_r (Raw_Alpha) to align
  corr_extract_r <-
    corr_extract_r %>%
    dplyr::select(-Raw_Alpha)

  # Restricting corr_extract to 6 columns
  corr_extract_r_second_part <-
    corr_extract_r$cities[4:6]

  # Restricting corr_extract to 6 columns
  corr_extract_sas_second_part <-
    as.numeric(corr_extract_sas$cities[4:6])

  # comparing ttest for no group by object
  testthat::expect_equal(object = corr_extract_sas_second_part,
                         expected = corr_extract_r_second_part,
                         tolerance = 0.01,
                         ignore_attr = TRUE)

})


#### All Row Variable Testing
testthat::test_that("Correlation for different age groups and education levels Part 1",{

  # Add readxl
  withr::local_package("readxl")

  # load corr_extract_sas output
  corr_extract_sas <-
    readxl::read_excel(
      path = testthat::test_path("sas_tests", "corrextractcomp2.xlsx"),
      sheet = "ALL")

  # clean up the column names in correlation
  corr_extract_sas <-
    corr_extract_sas %>%
    janitor::clean_names()

  # Convert COHORT to a factor for analysis
  finalBif <- bif_final %>%
    dplyr::mutate(CohortFactor = factor(COHORT))

  # converting to numeric from character to run correlation
  sim_data$individuals <- as.numeric(sim_data$individuals)
  sim_data_robust_reg$individuals <- as.numeric(sim_data_robust_reg$individuals)

  # corr extract r
  corr_extract_r <-
    corr_extract(sim_data_robust_reg, .vars = c("age_raw", "eduLevel"))

  # rename corr_extract_sas to align with
  corr_extract_sas <-
    corr_extract_sas %>%
    dplyr::rename(
      "Variable" = "variable",
      "Type" = "type",
      "Total_Raw_Alpha" = "raw_alpha",
      "Std_Alpha" = "std_corr_with_total",
      "Total_Std_Alpha" = "std_alpha",
      "mean" = "mean",
      "age_raw" = "age_raw",
      "eduLevel" = "edu_level"
    ) %>%
    dplyr::select(
      Variable,
      Type,
      Std_Alpha,
      Total_Raw_Alpha,
      Total_Std_Alpha,
      mean,
      age_raw,
      eduLevel)


  # filtering sas and r
  corr_extract_sas <-
    corr_extract_sas %>%
      dplyr::select(Variable, Type, Std_Alpha,
                  Total_Raw_Alpha, Total_Std_Alpha,
                  mean)

  corr_extract_r <-
    corr_extract_r %>%
      dplyr::select(Variable, Type, Std_Alpha,
                    Total_Raw_Alpha, Total_Std_Alpha,
                    mean)

  # comparing ttest for no group by object
  testthat::expect_equal(object = corr_extract_sas,
                         expected = corr_extract_r,
                         tolerance = 0.01,
                         ignore_attr = TRUE)

})


#### All Row Variable Testing
testthat::test_that("Correlation for different age groups and education levels Part 2",{

  # Add readxl
  withr::local_package("readxl")

  # load data set
  data("sim_data_robust_reg")

  # load ht_extract_sas output
  corr_extract_sas <-
    readxl::read_excel(
      path = testthat::test_path("sas_tests", "corrextractcomp2.xlsx"),
      sheet = "ALL")

  # clean up the column names in correlation
  corr_extract_sas <-
    corr_extract_sas %>%
    janitor::clean_names()

  # Convert COHORT to a factor for analysis
  finalBif <- bif_final %>%
    dplyr::mutate(CohortFactor = factor(COHORT))

  # converting to numeric from character to run correlation
  sim_data$individuals <- as.numeric(sim_data$individuals)
  sim_data_robust_reg$individuals <- as.numeric(sim_data_robust_reg$individuals)

  # corr extract r
  corr_extract_r <-
    corr_extract(sim_data_robust_reg, .vars = c("age_raw", "eduLevel"))

  # rename corr_extract_sas to align with
  corr_extract_sas <-
    corr_extract_sas %>%
    dplyr::rename(
      "Variable" = "variable",
      "Type" = "type",
      "Total_Raw_Alpha" = "raw_alpha",
      "Std_Alpha" = "std_corr_with_total",
      "Total_Std_Alpha" = "std_alpha",
      "mean" = "mean",
      "age_raw" = "age_raw",
      "eduLevel" = "edu_level"
    ) %>%
    dplyr::select(
      Variable,
      Type,
      Std_Alpha,
      Total_Raw_Alpha,
      Total_Std_Alpha,
      mean,
      age_raw,
      eduLevel)


  # Restricting corr_extract to 6 columns
  corr_extract_r_second_part <-
    corr_extract_r$age_raw[4:6]

  # Restricting corr_extract to 6 columns
  corr_extract_sas_second_part <-
    as.numeric(corr_extract_sas$age_raw[4:6])

  # comparing ttest for no group by object
  testthat::expect_equal(object = corr_extract_sas_second_part,
                         expected = corr_extract_r_second_part,
                         tolerance = 0.01,
                         ignore_attr = TRUE)

})


#### All Row Variable Testing
testthat::test_that("Correlation for different age groups and education levels by
                    within gender and number of children",{

  # Add readxl
  withr::local_package("readxl")

  # load data set
  data("sim_data_robust_reg")

  # load ht_extract_sas output
  corr_extract_sas <-
    readxl::read_excel(
      path = testthat::test_path("sas_tests", "corrextractcomp3.xlsx"),
      sheet = "ALL")

  # clean up the column names in correlation
  corr_extract_sas <-
    corr_extract_sas %>%
    janitor::clean_names()

  # Convert COHORT to a factor for analysis
  finalBif <- bif_final %>%
    dplyr::mutate(CohortFactor = factor(COHORT))

  # converting to numeric from character to run correlation
  sim_data$individuals <- as.numeric(sim_data$individuals)
  sim_data_robust_reg$individuals <- as.numeric(sim_data_robust_reg$individuals)

  # corr extract r
  corr_extract_r <-
    corr_extract(sim_data_robust_reg, .vars = c("age_raw", "eduLevel"),
                 .with_vars = c("gender", "numChild"))

  # dropping one column from corr_extract_r (Raw_Alpha) to align
  corr_extract_r <-
    corr_extract_r %>%
    dplyr::select(Variable,
                  mean,
                  Type,
                  age_raw,
                  eduLevel) %>%
    dplyr::arrange(Variable, Type) %>%
    janitor::clean_names()

  # rearrange sas to get the comparison close to correct
  corr_extract_sas <-
    corr_extract_sas %>%
      dplyr::arrange(variable, type)

  # comparing ttest for no group by object
  testthat::expect_equal(object = corr_extract_sas,
                         expected = corr_extract_r,
                         tolerance = 0.01,
                         ignore_attr = TRUE)

})

#### All Row Variable Testing
testthat::test_that("Correlation for different age groups and education levels by
                    within gender and number of children",{
  # Add readxl
  withr::local_package("readxl")

  # load data set
  data("sim_data_robust_reg")

  # load ht_extract_sas output
  corr_extract_sas <-
    readxl::read_excel(
      path = testthat::test_path("sas_tests", "corrextractcomp4.xlsx"),
      sheet = "0")

  # clean up the column names in correlation
  corr_extract_sas <-
    corr_extract_sas %>%
    janitor::clean_names()

  # Convert COHORT to a factor for analysis
  finalBif <- bif_final %>%
    dplyr::mutate(CohortFactor = factor(COHORT))

  # converting to numeric from character to run correlation
  sim_data$individuals <- as.numeric(sim_data$individuals)
  sim_data_robust_reg$individuals <- as.numeric(sim_data_robust_reg$individuals)

  # corr extract r
  corr_extract_r <-
    corr_extract(sim_data_robust_reg,
                 .vars = c("age_raw", "eduLevel"),
                 .with_vars = c("gender"),
                 .by_x = "treatment")

  # dropping one column from corr_extract_r (Raw_Alpha) to align
  corr_extract_r <-
    corr_extract_r %>%
    dplyr::select(Variable,
                  mean,
                  Type,
                  age_raw,
                  eduLevel) %>%
    dplyr::arrange(Variable, Type) %>%
    janitor::clean_names()

  # rearrange sas to get the comparison close to correct
  corr_extract_sas <-
    corr_extract_sas %>%
    dplyr::arrange(variable, type)

  # comparing ttest for no group by object
  testthat::expect_equal(object = corr_extract_sas,
                         expected = corr_extract_r,
                         tolerance = 0.01,
                         ignore_attr = TRUE)

  })





































