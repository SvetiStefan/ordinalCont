#' ordinalCont.
#'
#' @name ordinalCont.
#' @docType package
NULL

#' VAS scores marked by 84 patients at three time marks.
#'
#' A dataset containing the VAS scores of 84 patients at three time marks: baseline, the end of the course of treatment and one month after completion of the treatment. The variables are as follows:
#'
#' \itemize{
#'   \item time. the time mark at which the VAS score is recorded.
#'   \item ID. The ID number of the patient.
#'   \item sex Sex of patient.
#'   \item age Age of patient at entry in the study.
#'   \item group Trial group (2002 or 2003).
#'   \item duration duration of ...
#'   \item laterali uni- or bi-lateral pain.
#'   \item laser 1 if the laser was active, 0 otherwise.
#'   \item localisa diffuse or lacalized pain.
#'   \item vas. VAS score (0--1).
#'   \item lasert1. dummy variable. 1 if the laser was active when the VAS score has been recorded at baseline, 0 otherwise.
#'   \item lasert2. dummy variable. 1 if the laser was active when the VAS score has been recorded at the end of the course of treatment, 0 otherwise.
#'   \item lasert3. dummy variable. 1 if the laser was active when the VAS score has been recorded one month after completion of the treatment, 0 otherwise.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name pain
#' @usage data(pain)
#' @format A data frame with 252 rows and 13 variables
NULL

#' QoL scores marked by 80 patients at five time marks (only baseline included).
#'
#' A dataset containing the QoL scores of 80 patients at five time marks: baseline, 2, 4, 6, 12 months after a stroke event. The variables are as follows:
#'
#' \itemize{
#'   \item ID. The ID number of the patient.
#'   \item Age Age of patient at entry in the study.
#'   \item Gender Sex of patient.
#'   \item Height
#'   \item Weight
#'   \item BMI
#'   \item Side.Affected
#'   \item Time...stroke duration of ...
#'   \item Sensation ...
#'   \item Neglect ...
#'   \item Spasticity ...
#'   \item Contracture ...
#'   \item M1.EQoL.health QoL score (from 0 to 100)
#' }
#'
#' @docType data
#' @keywords datasets
#' @name qol
#' @usage data(qol)
#' @format A data frame with 80 rows and 13 variables
NULL

