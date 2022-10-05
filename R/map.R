get_pars <- function(p_vector,dadm) {
  # Add constants, transform p_vector, map to design, transform mapped parameters
  # to the natural scale, and create trial-dependent parameters
  attr(dadm,"model")$Ttransform(
    attr(dadm,"model")$Ntransform(
      map_p(
        attr(dadm,"model")$transform(add_constants(p_vector,attr(dadm,"constants"))),
        dadm)),
    dadm)
}

add_constants <- function(p,constants)
  # augments parameter matrix or vector p with constant parameters (also used in data)
{
  if (is.null(constants)) return(p)
  if (is.matrix(p)) {
    nams <- c(dimnames(p)[[2]],names(constants))
    p <- cbind(p,matrix(rep(constants,each=dim(p)[1]),nrow=dim(p)[1]))
    dimnames(p)[[2]] <- nams
    p
  } else c(p,constants)
}
