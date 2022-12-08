#' Summarize deep sequenced data into overall infection status by well
#' @name any_DVL
#' @param assay Assay data, with rows representing the distinct viral lineages (DVL) and columns representing the wells.
#' @return A vector of length \code{nrow(assay)}
#'
any_DVL = function(assay)  {
  assay = apply(X = assay, MARGIN = 2, FUN = max, na.rm = FALSE)
  assay[is.na(assay)] = 1
  return(assay)
}