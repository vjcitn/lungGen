#' use MAPPED_TRAIT_URI to find hits after mapping text to EBI ontology tag
#' @param interm character(1) term of interest in English
#' @param gwc instance of gwaswloc
#' @examples
#' data(lgenGWC_17)
#' ph = phenoToHits("pulmonary", lgenGWC_17)
#' ph
#' table(ph$MAPPED_TRAIT)
#' @export
phenoToHits = function (interm, gwc) 
{
    data(taggedPhenoDF)
    tag = (taggedPhenoDF %>% dplyr::filter(grepl(interm, term)) %>% dplyr::select(tag))[[1]]
    stopifnot(length(tag)>0)
    tag = colon2und(tag)
    inds = unlist(lapply(tag, function(x) grep(x, gwc$MAPPED_TRAIT_URI)))
    gwc[unlist(inds)]
}
