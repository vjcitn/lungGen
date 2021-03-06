#' @importFrom GenomicFiles files
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb seqnames
statesOLD = function(snpgr, ermaset, genome="GRCh38") {
  rowRanges(ermaset) = snpgr
  efil = files(ermaset)
  sts = bplapply(efil, function(x) {
   imp = import(x, which=rowRanges(ermaset), genome=genome) # states overlapped by SNP
   chkback = subsetByOverlaps(as(rowRanges(ermaset), "GRanges"), imp)
   ans = data.frame(
     state = 
       as.character(imp$name),
     cell = 
       as.character(short_celltype[cellTypes(ermaset)[match(x, efil)]]),
     rsid = names(chkback),
     snploc = start(chkback)
     ) 
    rownames(ans) = NULL
    ans$start.state = start(imp)
    ans$end.state = end(imp)
    ans$seqnames = as.character(GenomeInfoDb::seqnames(imp))
    ans
    }
   )
  ans = do.call(rbind, sts)
  rownames(ans) = NULL
  ans
  }

states = function (snpgr, ermaset, genome = "GRCh38") 
{
    cellTypes2 = function(x) colData(x)[,"Standardized.Epigenome.name"]
    rowRanges(ermaset) = snpgr
    efil = files(ermaset)
    sts = bplapply(efil, function(x) {
        imp = import(x, which = rowRanges(ermaset), genome = genome)
        chkback = subsetByOverlaps(as(rowRanges(ermaset), "GRanges"), 
            imp)
        ans = data.frame(state = as.character(imp$name), cell = as.character(short_celltype[cellTypes2(ermaset)[match(x, 
            efil)]]), rsid = names(chkback), snploc = start(chkback))
        rownames(ans) = NULL
        ans$start.state = start(imp)
        ans$end.state = end(imp)
        ans$seqnames = as.character(GenomeInfoDb::seqnames(imp))
        ans
    })
    ans = do.call(rbind, sts)
    rownames(ans) = NULL
    ans
}

