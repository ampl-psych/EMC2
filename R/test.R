donothing <- function(i) repeat{}

pdonothing <- function(i)
  parallel:mclapply(1:3,donothing,mc.cores=3)

#' Title
#'
#' @returns
#' @export
#'
#' @examples
test_parallel <- function() {
  ncores=3
  parallel:mclapply(1:ncores,pdonothing,mc.cores=ncores)
}
