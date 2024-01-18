#' Forstmann et al.'s data
#'
#' A dataset containing the speed or accuracy manipulation for a Random Dot
#' Motion experiment.
#'
#' Details on the dataset can be found in the following paper:
#'
#' \strong{Striatum and pre-SMA facilitate decision-making under time pressure}
#'
#' Birte U. Forstmann, Gilles Dutilh, Scott Brown, Jane Neumann,
#' D. Yves von Cramon, K. Richard Ridderinkhof, Eric-Jan Wagenmakers.
#'
#' \emph{Proceedings of the National Academy of Sciences Nov 2008, 105 (45)
#' 17538-17542; DOI: 10.1073/pnas.0805903105}
#'
#' @format A data frame with 15818 rows and 5 variables:
#' \describe{
#'   \item{subject}{integer ID for each subject}
#'   \item{rt}{reaction time for each trial as a double}
#'   \item{condition}{Factor with 3 levels for Speed, Accuracy and
#'     Neutral}
#'   \item{stim}{Factor with 2 levels for Left and Right trials}
#'   \item{resp}{Factor with 2 levels for Left and Right responses}
#' }
#' @source \url{https://www.pnas.org/content/105/45/17538}
"forstmann"
