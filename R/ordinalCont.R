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

#' @title ANZ0001 trial
#' The ANZ0001 trial, conducted by the ANZ Breast Cancer Trials Group, is an unblinded, multi-centre, randomised trial with three chemotherapy treatment arms concluded in 2005 (Stockler, Sour- jina, Grimison, Gebski, Byrne, Harvey, Francis, Nowak, Hazel, Forbes, and Group, 2007). 
#' Health-related quality of life measures (Overall quality of life, Physical Well-Being, Mood, Pain, Nausea and Vomiting, Appetite) are assessed at each chemotherapy treatment cycle, from randomization until disease progression, when treatment is interrupted. 
#' The treatments Intermittent Capecitabine (IC) and Continuous Capecitabine (CC) are compared with the standard combination treatment CMF, each with its own protocol. 
#' There is no maximum duration of treatment, but it is interrupted on disease progression, or when patient intolerance or unacceptable toxicity are recorded.
#' The data set is extracted from the ANZ0001 trial and contains information from 292 patients with complete quality of life measurements.
#'
#' The variables are as follows:
#'
#' \itemize{
#'   \item randno. The ID number of the patient.
#'   \item cycleno The chemo cycle number.
#'   \item age Age of patient at entry in the study.
#'   \item bsa Patient Body Surface Area (m2).
#'   \item treatment The treatment received by the patient.
#'   \item overall Overall quality of life as recorded by the patient on a LASA scale, here normalised in the interval (0, 1).
#'   \item phys Physical Well-Being as recorded by the patient on a LASA scale, here normalised in the interval (0, 1).
#'   \item mood Mood as recorded by the patient on a LASA scale, here normalised in the interval (0, 1).
#'   \item pain Pain as recorded by the patient on a LASA scale, here normalised in the interval (0, 1).
#'   \item nausvom Nausea and Vomiting as recorded by the patient on a LASA scale, here normalised in the interval (0, 1).
#'   \item appetite Appetite as recorded by the patient on a LASA scale, here normalised in the interval (0, 1).
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ANZ0001
#' @usage data(ANZ0001)
#' @format A data frame with 3473 rows and 11 variables
NULL

