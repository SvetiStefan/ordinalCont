#' ordinalCont.
#'
#' @name ordinalCont.
#' @docType package
NULL

#' @title Pain data 
#' 
#' @details A dataset containing the VAS scores in a clinical trial for the effectiveness
#' of laser treatment for chronic neck pain. Patients (n=84) were randomly assigned to receive 
#' either active laser treatment or placebo. VAS scores were recorded by the patients
#'  at three time points: baseline, 
#' the end of the course of treatment and one month after completion of the treatment. 
#' 
#' The variables are as follows:
#' \tabular{rlll}{
#'  \tab \bold{Name}\tab \bold{Type}\tab \bold{Description}\cr
#'    \tab time\tab numeric\tab  time  at which the VAS score is recorded: 1=baseline; 2=end of treatment; \cr
#'    \tab\tab\tab 3=one month after end of  treatment\cr
#'   \tab ID\tab numeric\tab The ID number of the patient\cr
#'   \tab sex\tab binary\tab 1=male; 2=female\cr
#'   \tab age\tab numeric\tab age of patient at entry in the study\cr
#'   \tab group\tab  numeric \tab trial group (2002 or 2003)\cr
#'   \tab duration\tab numeric \tab duration of ...\cr
#'   \tab laterali \tab binary\tab uni- or bi-lateral pain\cr
#'   \tab laser \tab binary\tab 1 if the laser was active, 0 otherwise\cr
#'   \tab localisa\tab ??\tab diffuse or localized pain\cr
#'   \tab vas\tab numeric\tab VAS score (0--1)\cr
#'   \tab lasert1\tab binary\tab 1 if the laser was active and the VAS score was recorded at baseline, 0 otherwise\cr
#'   \tab lasert2\tab binary\tab 1 if the laser was active and the VAS score was recorded at the end of the course\cr
#'   \tab\tab\tab of treatment, 0 otherwise\cr
#'   \tab lasert3\tab binary\tab 1 if the laser was active and the VAS score was recorded one month after completion\cr
#'   \tab\tab\tab of the treatment, 0 otherwise\cr
#'  }
#'  
#'  @references Chow RT, Heller GZ, Barnsley L (2006). The effect of 300mW, 830nm laser on 
#'  chronic neck pain: A double-blind, randomized, placebo-controlled study, \emph{Pain}, 
#'  124: 201-210.

#'
#' @docType data
#' @keywords datasets
#' @name pain
#' @usage data(pain)
#' @format A data frame with 252 rows and 13 variables
NULL

#' @title ANZ0001 trial
#' 
#' @details  The ANZ0001 trial, conducted by the ANZ Breast Cancer Trials Group, is an unblinded, multi-centre, randomised trial with three chemotherapy treatment arms concluded in 2005 (Stockler, Sour- jina, Grimison, Gebski, Byrne, Harvey, Francis, Nowak, Hazel, Forbes, and Group, 2007). 
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

