#' BIF_Final
#'
#' A sample baseline information form data set
#'
#' @format ## `bif_final`
#' A data frame with 2,761 rows and 34 columns:
#' \describe{
#'   \item{SAMPLEID}{Sample ID}
#'   \item{COHORT}{Cohort}
#'   \item{RA_CODE}{Program or Control Status}
#'   \item{RA_DATE}{Date of Random Assignment}
#'   \item{AGE}{Age at random assignment}
#'   \item{DOB}{Date of Birth}
#'   \item{blmale, blfemale}{Flag for gender}
#'   \item{blmsmarried_lws, blmsmarried_lwos, blmsunmarried_lwp, blmsunmmarried_lwop}{Flags for marital status}
#'   \item{blworkingyes, blworkingno}{Flags for working status}
#'   \item{blcurhrs_1to19, blcurhrs_20to29, blcurhrs_30to34, blcurhrs_35to39, blcurhrs_40to72, blcurhrs_73plus}{Flags for current hours worked}
#'   \item{blhavechild_lt18}{Flag for children under 18}
#'   \item{blhgrade_8ol, blhgrade_9,blhgrade_10,blhgrade_11,blhgrade_12}{Flags for final year of school completed}
#'   \item{bldiplomas_hs, bldiplomas_ged, bldiplomas_tec, bldiplomas_as, bldiplomas_4yr, bldiplomas_md, bldiplomas_non}{Flags for highest level of education}
#'   \item{blotherlang}{Flag for other language spoken at home}
#' }
#' @source Sample data created internally
"bif_final"

#' sim_data
#'
#' A simulated data set
#'
#' @format ## `sim_data`
#' A data frame with 323 rows and 20 columns:
#' \describe{
#'   \item{cities}{cities ID}
#'   \item{population}{population of city}
#'   \item{sites}{Site ID}
#'   \item{individuals}{ID}
#'   \item{age_raw, age_scaled}{Age}
#'   \item{eduLevel}{Flag for education level}
#'   \item{pre_income_raw, pre_income_scaled, pre_income_raw_25pct}{income}
#'   \item{tanf_receipt, snap_receipt}{Benefits}
#'   \item{gender}{gender}
#'   \item{numChild}{children under 18}
#'   \item{maritalStatus}{Marital Status}
#'   \item{treatment}{Program or Control Status}
#'   \item{probs_01, probs_02}{probabilities}
#'   \item{employed_01, employed_02}{employment status}
#' }
#' @source Sample data created internally
"sim_data"

#' sim_data_robust_reg
#'
#' A simulated data set
#'
#' @format ## `sim_data`
#' A data frame with 323 rows and 21 columns:
#' \describe{
#'   \item{cities}{cities ID}
#'   \item{population}{population of city}
#'   \item{sites}{Site ID}
#'   \item{individuals}{ID}
#'   \item{age_raw, age_scaled}{Age}
#'   \item{eduLevel}{Flag for education level}
#'   \item{pre_income_raw, pre_income_scaled, pre_income_raw_25pct}{income}
#'   \item{tanf_receipt, snap_receipt}{Benefits}
#'   \item{gender}{gender}
#'   \item{numChild}{children under 18}
#'   \item{maritalStatus}{Marital Status}
#'   \item{treatment}{Program or Control Status}
#'   \item{probs_01, probs_02}{probabilities}
#'   \item{employed_01, employed_02}{employment status}
#'   \item{favAnimals}{favorite animals}
#' }
#' @source Sample data created internally
"sim_data_robust_reg"

#' sim_data_weighted
#'
#' A simulated data set
#'
#' @format ## `sim_data`
#' A data frame with 323 rows and 21 columns:
#' \describe{
#'   \item{cities}{cities ID}
#'   \item{population}{population of city}
#'   \item{sites}{Site ID}
#'   \item{individuals}{ID}
#'   \item{age_raw, age_scaled}{Age}
#'   \item{eduLevel}{Flag for education level}
#'   \item{pre_income_raw, pre_income_scaled, pre_income_raw_25pct}{income}
#'   \item{tanf_receipt, snap_receipt}{Benefits}
#'   \item{gender}{gender}
#'   \item{numChild}{children under 18}
#'   \item{maritalStatus}{Marital Status}
#'   \item{treatment}{Program or Control Status}
#'   \item{probs_01, probs_02}{probabilities}
#'   \item{employed_01, employed_02}{employment status}
#'   \item{genWeight}{generated Weights}
#' }
#' @source Sample data created internally
"sim_data_weighted"

#' sim_three_levels_missing
#'
#' A simulated data set with missing data
#'
#' @format ## `sim_data`
#' A data frame with 323 rows and 21 columns:
#' \describe{
#'   \item{cities}{cities ID}
#'   \item{population}{population of city}
#'   \item{sites}{Site ID}
#'   \item{individuals}{ID}
#'   \item{age_raw, age_scaled}{Age}
#'   \item{eduLevel}{Flag for education level}
#'   \item{pre_income_raw, pre_income_scaled, pre_income_raw_25pct}{income}
#'   \item{tanf_receipt, snap_receipt}{Benefits}
#'   \item{gender}{gender}
#'   \item{numChild}{children under 18}
#'   \item{maritalStatus}{Marital Status}
#'   \item{treatment}{Program or Control Status}
#'   \item{probs_01, probs_02}{probabilities}
#'   \item{employed_01, employed_02}{employment status}
#'   \item{favAnimals}{favorite Animals}
#' }
#' @source Sample data created internally
"sim_three_levels_missing"


