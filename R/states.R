states = function(snpgr, ermaset, genome="GRCh38") {
  rowRanges(ermaset) = snpgr
  efil = files(ermaset)
  sts = bplapply(efil, function(x) {
   imp = import(x, which=rowRanges(ermaset), genome=genome) # states overlapped by SNP
   chkback = subsetByOverlaps(as(rowRanges(ermaset), "GRanges"), imp)
   ans = data.frame(
     state = 
       #as.character(erma:::liberalImport(x, which=rowRanges(ermaset), genome=genome)$name),
       as.character(imp$name),
     cell = 
       as.character(short_celltype[cellTypes(ermaset)[match(x, efil)]]),
     rsid = names(chkback),
     snploc = start(chkback)
     ) 
    rownames(ans) = NULL
    ans$start.state = start(imp)
    ans$end.state = end(imp)
    ans$seqnames = as.character(seqnames(imp))
    ans
    }
   )
  ans = do.call(rbind, sts)
  rownames(ans) = NULL
  ans
  }
