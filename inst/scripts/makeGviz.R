library(Gviz)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(lungGen)
data(predGR)
ugenes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)
id = ugenes$gene_id
ss = select(org.Hs.eg.db, keys=id, keytype="ENTREZID", columns="SYMBOL")$SYMBOL
ugenes = ugenes[-which(is.na(ss))]
ugenes$symbol = ss[-which(is.na(ss))]
gs = subsetByOverlaps(predGR, GRanges("chr17", IRanges(39.9e6,40e6)))
genome(gs) = "hg38"
subsetByOverlaps(ugenes, range(gs)) -> hh
plotTracks(list(GenomeAxisTrack(name="chr17"), 
  DataTrack(gs, name="phrProb", ylim=c(.7,.87)), 
  GeneRegionTrack(hh, showId=TRUE, name="asLoc", shape="arrow")))
