testthat::test_that("Distribution of 2 Categorical Variables (Unit test)",{

  twoCatsDs <- freq_extract(
    .dataset = sim_data, # data set you have loaded above
    .freq_vars = c("gender", "numChild"), # categorical variables you want to calculate a frequency distribution
    .max_lvls = 10, # maximum number of expected levels
    .round = 0) # rounding for percentage estimation

  #-------------#
  firstCat <- "gender"
  secondCat <- "numChild"
  #-------------#
  testSet01 <-
    sim_data %>%
    dplyr::group_by(gender) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Value" = "gender",
      "Freq" = "n") %>%
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 0),
      CPct = cumsum(Pct),
      Variable = "gender",
      Unique_ID = paste0(Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Variable,
      Value,
      Freq,
      CFreq,
      Pct,
      CPct
    )
  testSet02 <-
    sim_data %>%
    dplyr::group_by(numChild) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Value" = "numChild",
      "Freq" = "n") %>%
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 0),
      CPct = cumsum(Pct),
      Variable = "numChild",
      Unique_ID = paste0(Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Variable,
      Value,
      Freq,
      CFreq,
      Pct,
      CPct
    )
  testSet <- dplyr::bind_rows(testSet01,
                              testSet02)
  testSet$Value <- as.character(testSet$Value)

  testthat::expect_equal(testSet,
                         twoCatsDs)


})

testthat::test_that("Distribution of 2 Categorical Variables (Decimals)",{

  #-------------#
  firstCat <- "gender"
  secondCat <- "numChild"
  #-------------#
  testSet01 <-
    sim_data %>%
    dplyr::group_by(gender) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Value" = "gender",
      "Freq" = "n") %>%
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 5),
      CPct = cumsum(Pct),
      Variable = "gender",
      Unique_ID = paste0(Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Variable,
      Value,
      Freq,
      CFreq,
      Pct,
      CPct
    )
  testSet02 <-
    sim_data %>%
    dplyr::group_by(numChild) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Value" = "numChild",
      "Freq" = "n") %>%
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 5),
      CPct = cumsum(Pct),
      Variable = "numChild",
      Unique_ID = paste0(Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Variable,
      Value,
      Freq,
      CFreq,
      Pct,
      CPct
    )
  testSet <- dplyr::bind_rows(testSet01,
                              testSet02)
  testSet$Value <- as.character(testSet$Value)

  twoCatsDsDecis <- freq_extract(
    .dataset = sim_data, # data set you have loaded above
    .freq_vars = c("gender", "numChild"), # categorical variables you want to calculate a frequency distribution
    .max_lvls = 10, # maximum number of expected levels
    .round = 2) # rounding for percentage estimation

  testthat::expect_equal(object = testSet,
                         expected = twoCatsDsDecis,
                         tolerance = 0.05)

})

testthat::test_that("Distribution of 2 Categorical Variables within a group (i.e. Cities)",{

  twoCatsWGroupDsDecis <- freq_extract(
    .dataset = sim_data,                  # data set you have loaded above
    .freq_vars = c("gender", "numChild"), # categorical variables you want to calculate a frequency distribution
    .subgroup = c("cities"),              # the subgroup by which we want to divide by
    .max_lvls = 10,                       # maximum number of expected levels
    .round = 2)                           # rounding for percentage estimation

  twoCatsWGroupDsDecis <-
    twoCatsWGroupDsDecis %>%
    dplyr::arrange(Unique_ID) %>%
    dplyr::select(-Subgrp_Var)

  testSet01 <-
    sim_data %>%
    dplyr::group_by(cities, gender) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Subgrp_Value" = "cities",
      "Value" = "gender",
      "Freq" = "n") %>%
    #dplyr::ungroup() %>% # Critical
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 5),
      CPct = cumsum(Pct),
      Variable = "gender",
      Unique_ID = paste0(Subgrp_Value, ":", Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Subgrp_Value,
      Variable,
      Value,
      Freq,
      CFreq,
      Pct,
      CPct
    )

  testSet02 <-
    sim_data %>%
    dplyr::group_by(cities, numChild) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Subgrp_Value" = "cities",
      "Value" = "numChild",
      "Freq" = "n") %>%
    #dplyr::ungroup() %>% # Critical
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 5),
      CPct = cumsum(Pct),
      Variable = "numChild",
      Unique_ID = paste0(Subgrp_Value, ":", Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Subgrp_Value,
      Variable,
      Value,
      Freq,
      CFreq,
      Pct,
      CPct
    )

  testSet <- dplyr::bind_rows(testSet01,
                              testSet02)

  testSet$Value <- as.character(testSet$Value)
  testSet <- testSet %>%
    dplyr::arrange(Unique_ID)

  testthat::expect_equal(object = testSet,
                         expected = twoCatsWGroupDsDecis,
                         tolerance = 0.05)

})

testthat::test_that("Distribution of a Continuous variable - Stage 1 test",{

  # Example 1
  contDsDecis01 <- freq_extract(
    .dataset = sim_data, # data set you have loaded above
    .freq_vars = c("population", "gender", "numChild"), # categorical variables you want to calculate a frequency distribution
    .max_lvls = 10, # maximum number of expected levels
    .round = 2) # rounding for percentage estimation
  contDsDecis01

  # Example 2
  contDsDecis03 <- freq_extract(
    .dataset = sim_data, # data set you have loaded above
    .max_lvls = 4, # maximum number of expected levels
    .round = 2) # rounding for percentage estimation
  contDsDecis03

  #### Does the function filter to the right number of categorical variables?
  colLevels <- as.data.frame(t(as.data.frame(lapply(apply(sim_data, 2, unique), length))))
  names(colLevels) <- "levelCount"
  colLevels

  colLevelsFilter <-
    colLevels %>%
    dplyr::arrange(levelCount) %>%
    dplyr::filter(levelCount <= 4)
  testSet <- base::rownames(colLevelsFilter)
  testSet

  functionSet <- unique(contDsDecis03$Variable)
  functionSet

  testthat::expect_equal(object = testSet ,
                         expected = functionSet)

})

testthat::test_that("Distribution of 2 Categorical Variables with Weights - Stage 2 test",{

  # Example 1
  contDsDecis01 <- freq_extract(
    .dataset = sim_data, # data set you have loaded above
    .freq_vars = c("population", "gender", "numChild"), # categorical variables you want to calculate a frequency distribution
    .max_lvls = 10, # maximum number of expected levels
    .round = 2) # rounding for percentage estimation
  contDsDecis01

  # Example 2
  contDsDecis03 <- freq_extract(
    .dataset = sim_data, # data set you have loaded above
    .max_lvls = 4, # maximum number of expected levels
    .round = 2) # rounding for percentage estimation
  contDsDecis03

  # Filtered categorical variables
  #-------------#
  firstCat <- "pre_income_raw_25pct"
  secondCat <- "tanf_receipt"
  thirdCat <- "snap_receipt"
  fourthCat <- "gender"
  fifthCat <- "maritalStatus"
  sixthCat <- "treatment"
  seventhCat <- "employed_01"
  eighthCat <- "employed_02"
  #-------------#
  testSet01 <-
    sim_data %>%
    dplyr::group_by(pre_income_raw_25pct) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Value" = "pre_income_raw_25pct",
      "Freq" = "n") %>%
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 5),
      CPct = cumsum(Pct),
      Variable = "pre_income_raw_25pct",
      Unique_ID = paste0(Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Variable,
      Value,
      Freq,
      CFreq,
      Pct,
      CPct
    )
  testSet02 <-
    sim_data %>%
    dplyr::group_by(tanf_receipt) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Value" = "tanf_receipt",
      "Freq" = "n") %>%
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 5),
      CPct = cumsum(Pct),
      Variable = "tanf_receipt",
      Unique_ID = paste0(Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Variable,
      Value,
      Freq,
      CFreq,
      Pct,
      CPct
    )
  testSet03 <-
    sim_data %>%
    dplyr::group_by(snap_receipt) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Value" = "snap_receipt",
      "Freq" = "n") %>%
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 5),
      CPct = cumsum(Pct),
      Variable = "snap_receipt",
      Unique_ID = paste0(Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Variable,
      Value,
      Freq,
      CFreq,
      Pct,
      CPct
    )
  testSet04 <-
    sim_data %>%
    dplyr::group_by(gender) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Value" = "gender",
      "Freq" = "n") %>%
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 5),
      CPct = cumsum(Pct),
      Variable = "gender",
      Unique_ID = paste0(Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Variable,
      Value,
      Freq,
      CFreq,
      Pct,
      CPct
    )
  testSet05 <-
    sim_data %>%
    dplyr::group_by(maritalStatus) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Value" = "maritalStatus",
      "Freq" = "n") %>%
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 5),
      CPct = cumsum(Pct),
      Variable = "maritalStatus",
      Unique_ID = paste0(Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Variable,
      Value,
      Freq,
      CFreq,
      Pct,
      CPct
    )
  testSet06 <-
    sim_data %>%
    dplyr::group_by(treatment) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Value" = "treatment",
      "Freq" = "n") %>%
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 5),
      CPct = cumsum(Pct),
      Variable = "treatment",
      Unique_ID = paste0(Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Variable,
      Value,
      Freq,
      CFreq,
      Pct,
      CPct
    )
  testSet07 <-
    sim_data %>%
    dplyr::group_by(employed_01) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Value" = "employed_01",
      "Freq" = "n") %>%
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 5),
      CPct = cumsum(Pct),
      Variable = "employed_01",
      Unique_ID = paste0(Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Variable,
      Value,
      Freq,
      CFreq,
      Pct,
      CPct
    )
  testSet08 <-
    sim_data %>%
    dplyr::group_by(employed_02) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Value" = "employed_02",
      "Freq" = "n") %>%
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 5),
      CPct = cumsum(Pct),
      Variable = "employed_02",
      Unique_ID = paste0(Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Variable,
      Value,
      Freq,
      CFreq,
      Pct,
      CPct
    )
  testSet <- dplyr::bind_rows(testSet01,
                              testSet02,
                              testSet03,
                              testSet04,
                              testSet05,
                              testSet06,
                              testSet07,
                              testSet08)
  testSet$Value <- as.character(testSet$Value)

  testthat::expect_equal(object = testSet ,
                         expected = contDsDecis03,
                         tolerance = 0.05)

})

testthat::test_that("Distribution of 2 Categorical Variables with Weights",{

  sim_data_weighted <- sim_data %>%
    dplyr::mutate(genWeight = ifelse(gender %in% c(0), 2, 1))

  # two weighted categories
  twoCatsDsWeighted <- freq_extract(
    .dataset = sim_data_weighted, # data set you have loaded above
    .freq_vars = c("gender", "numChild"), # categorical variables you want to calculate a frequency distribution
    .wt_var = "genWeight",
    .max_lvls = 10, # maximum number of expected levels
    .round = 2) # rounding for percentage estimation

  # Weight for Gender(0/1) - 2 : 1
  #-------------#
  firstCat <- "gender"
  secondCat <- "numChild"
  #-------------#
  testSet01 <-
    sim_data_weighted %>%
    dplyr::group_by(gender) %>%
    dplyr::summarise(
      Freq = dplyr::n(),
      Wt_Freq = sum(genWeight)) %>%
    dplyr::rename(
      "Value" = "gender",
      "Freq" = "Freq",
      "Wt_Freq" = "Wt_Freq") %>%
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 2),
      CPct = cumsum(Pct),
      Wt_CFreq = cumsum(Wt_Freq),
      Wt_Pct = round(100 * (Wt_Freq/sum(Wt_Freq)),2),
      Wt_CPct = cumsum(Wt_Pct),
      Variable = "gender",
      Unique_ID = paste0(Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Variable,
      Value,
      Freq,
      Wt_Freq,
      CFreq,
      Pct,
      CPct,
      Wt_CFreq,
      Wt_Pct,
      Wt_CPct
    )

  testSet02 <-
    sim_data_weighted %>%
    dplyr::group_by(numChild) %>%
    dplyr::summarise(
      Freq = dplyr::n(),
      Wt_Freq = sum(genWeight)) %>%
    dplyr::rename(
      "Value" = "numChild",
      "Freq" = "Freq",
      "Wt_Freq" = "Wt_Freq") %>%
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 2),
      CPct = cumsum(Pct),
      Wt_CFreq = cumsum(Wt_Freq),
      Wt_Pct = round(100 * (Wt_Freq/sum(Wt_Freq)),2),
      Wt_CPct = cumsum(Wt_Pct),
      Variable = "numChild",
      Unique_ID = paste0(Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Variable,
      Value,
      Freq,
      Wt_Freq,
      CFreq,
      Pct,
      CPct,
      Wt_CFreq,
      Wt_Pct,
      Wt_CPct
    )

  testSet <- dplyr::bind_rows(testSet01,
                              testSet02)
  testSet$Value <- as.character(testSet$Value)

  # the count multiplication
  testthat::expect_equal(object = testSet ,
                         expected = twoCatsDsWeighted,
                         tolerance = 0.05)

})

testthat::test_that("Missing Values",{

  missingDataSet <-
    freq_extract(.dataset = sim_three_levels_missing, # data set you have loaded above
                 .freq_vars = c("numChild"),
                 .max_lvls = 10,
                 .na_rm = TRUE,
                 .round = 2)

  testSet01 <-
    sim_three_levels_missing %>%
    dplyr::group_by(numChild) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Value" = "numChild",
      "Freq" = "n") %>%
    dplyr::mutate(
      CFreq = cumsum(Freq),
      Pct = round(100 * (Freq/sum(Freq)), 5),
      CPct = cumsum(Pct),
      Variable = "numChild",
      Unique_ID = paste0(Variable, ":", Value)
    ) %>%
    dplyr::select(
      Unique_ID,
      Variable,
      Value,
      Freq,
      CFreq,
      Pct,
      CPct
    )

  testSet02 <-
    sim_three_levels_missing %>%
    dplyr::group_by(numChild) %>%
    dplyr::tally() %>%
    dplyr::rename(
      "Value" = "numChild",
      "Freq" = "n") %>%
    dplyr::mutate(
      Variable = "numChild",
      Unique_ID = paste0(Variable, ":", Value)) %>%
    dplyr::mutate(
      Nm = ifelse(is.na(Value), Freq, 0),
      Nm_Freq = Freq - Nm,
      Nm_CFreq = cumsum(Nm_Freq),
      Nm_Pct = round(100 * (Nm_Freq/sum(Nm_Freq)), 5),
      Nm_CPct = cumsum(Nm_Pct),
      Variable = "numChild",
    ) %>%
    dplyr::select(
      Nm,
      Nm_Freq,
      Nm_CFreq,
      Nm_Pct,
      Nm_CPct
    )

  testSet <- cbind(testSet01,
                   testSet02)
  testSet$Value <- as.character(testSet$Value)

  missingDataSet <- as.data.frame(missingDataSet)

  # the count multiplication
  testthat::expect_equal(object = testSet ,
                         expected = missingDataSet,
                         tolerance = 0.05)
})



testthat::test_that("Testing",{

  # sim_data weighted
  sim_data_weighted <- sim_data %>%
    dplyr::mutate(genWeight = ifelse(gender %in% c(0), 2, 1))

  #summarize sim_data_weighted
  sim_data_weighted_sum_genW <-
    sum(sim_data_weighted$genWeight)

  # two weighted categories
  # the count multiplication
  testthat::expect_equal(object = sim_data_weighted_sum_genW ,
                         expected = sum(sim_data_weighted$genWeight),
                         tolerance = 0.05)

})
