#

library(lungGen)
library(ggplot2)
library(plotly)
library(gggmvis)
regexpr="lung|asthma|pulmonary"
bpp=BiocParallel::SerialParam()
BiocParallel::register(bpp)
data(taggedPhenoDF)
data(lgenGWC_17)
GenomeInfoDb::seqlevelsStyle(lgenGWC_17) = "UCSC"
locErmaSet = makeEpig17()
cellTypes2 = function(x) colData(x)[,"Standardized.Epigenome.name"]
ct = cellTypes2(locErmaSet)
data(short_celltype, package="lungGen")
st = short_celltype[ct]
#


#' @importFrom GenomicFiles files
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb seqnames
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


server = function(input, output) {
  cellTypes2 = function(x) colData(x)[,"Standardized.Epigenome.name"]
  colnames(locErmaSet) = short_celltype[cellTypes2(locErmaSet)]
  output$cells = renderPrint( input$cellpicks )
  output$pms = renderPrint({   # some lists are cached 1/21/2019
   if (input$pheno == "asthma") {
     data(asthmaGWASpubs)
     ans = asthmaGWASpubs
     }
   else if (input$pheno == "chronic obstructive pulmonary disease") {
     data(copdGWASpubs)
     ans = copdGWASpubs
     }
   else if (input$pheno == "idiopathic pulmonary fibrosis") {
     data(idioGWASpubs)
     ans = idioGWASpubs
     }
   else {
     ids = unique(phenoToHits(input$pheno, lgenGWC_17)$PUBMEDID)
     ans = lapply(as.character(ids), pmid2MIAME)
     }
   ans
   })
  output$snps = DT::renderDataTable({
   hitgr = phenoToHits(input$pheno, lgenGWC_17)
   names(hitgr) = hitgr$SNPS
   dfr = as.data.frame(mcols(phenoToHits(input$pheno, lgenGWC_17))[,c("CHR_ID", "CHR_POS", "SNPS", "MAPPED_GENE", "PUBMEDID")])
   chkdup = apply(dfr, 1, paste0, collapse=":")
   todr = which(duplicated(chkdup))
   if (length(todr)>0) dfr = dfr[-todr,]
   dfr$CHR_POS = as.integer(dfr$CHR_POS)
   dfr$PUBMEDID = paste0("<a href=https://www.ncbi.nlm.nih.gov/pubmed/?term=", dfr$PUBMEDID, " target='_blank'>PMID ", dfr$PUBMEDID, "</a>")
   o = order(dfr$CHR_POS)
   DT::datatable(dfr[o,], escape=FALSE)
   })
  output$states = DT::renderDataTable({
    xx = phenoToHits(input$pheno, lgenGWC_17)
    names(xx) = xx$SNPS
    DT::datatable(states(xx, locErmaSet[, input$cellpicks]))
    })
  output$genesel = renderUI({
   dfr = as.data.frame(mcols(phenoToHits(input$pheno, lgenGWC_17))[,c("CHR_ID", "CHR_POS", "SNPS", "MAPPED_GENE", "PUBMEDID", "REPORTED GENE(S)")])
   chkdup = apply(dfr, 1, paste0, collapse=":")
   todr = which(duplicated(chkdup))
   if (length(todr)>0) dfr = dfr[-todr,]
   gns = as.character(dfr$REPORTED.GENE.S.)
   gns = sort(unique(unlist(strsplit(gns, ", "))))
   selectInput("gsel", "mapped gene", choices=gns, selected=gns[1])
   })
  output$txplot = renderPlotly({
    vis = try(ggvisForSymbol(input$gsel, arrmm=2.3))
    validate(need(vis, "can't build gene model; obsolete symbol?"))
#
# FIXME -- assuming v79 (hg38)
#
    g1 = genes(EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79, filter=AnnotationFilter::SymbolFilter(input$gsel))+1000
    GenomeInfoDb::seqlevelsStyle(g1) = "UCSC"
    preds = subsetByOverlaps(predGR, g1)
    stt = states(preds, locErmaSet[, input$cellpicks])
    save(stt, file="stt.rda")
    preddf = data.frame(start=start(preds), pred=preds$pred+1, rsid=names(preds), state=stt$state)
    preddf = preddf[preddf$pred > 1.03,] # FIXME
    ggplotly(vis + geom_point(data=preddf, 
        aes(x=start, y=pred, text=rsid, colour=state)) + 
        theme(axis.text.y=element_blank()) + ylab(" "))
    })
  } # end server
