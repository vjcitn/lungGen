#' mapping from MAPPED_TRAIT_URI tags to natural language terms
#' @importFrom utils data
#' @docType data
#' @format data.frame
#' @examples
#' lungGen::taggedPhenoDF
"taggedPhenoDF"
#' illustrative excerpt from EBI GWAS catalog for jan 2019; limit to chr17
#' @format gwaswloc instance as defined in gwascat package, GRCh38 addresses
#' @examples
#' lungGen::lgenGWC_17
"lgenGWC_17"
#' location, rsid and predicted probability of phenorelevance based on jan 30 2019 model
#' @format GRanges instance, GRCh38 addresses
#' @examples
#' head(predGR)
"predGR"
#' liftOver support for erma data in hg19
#' @format chain imported using rtracklayer
"ch19to38"
