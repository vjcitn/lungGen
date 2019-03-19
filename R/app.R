#' interface to txregnet resource for lung disease genomics
#' @import magrittr annotate gwascat shiny ggplot2
#' @importFrom dplyr filter select
#' @importFrom GenomicFeatures genes
#' @importFrom S4Vectors mcols
#' @importFrom IRanges subsetByOverlaps
#' @importFrom plotly plotlyOutput renderPlotly ggplotly
#' @importFrom BiocParallel register SerialParam bplapply
#' @importFrom AnnotationFilter SymbolFilter
#' @importFrom erma cellTypes
#' @importFrom gggmvis ggvisForSymbol
#' @importFrom SummarizedExperiment "rowRanges<-"
#' @importFrom DelayedArray "rowRanges"
# #' @param regexpr character(1) will be used to grep taggedPhenoDF$term
# #' @param bpp a BiocParallel bpparam instance
#' @export
lungGen = function() {
if (!requireNamespace("shiny")) stop("install shiny to use this function")
setwd(system.file("app", package="lungGen"))
shiny::runApp()
}

oldlungGen = function(regexpr="lung|asthma|pulmonary", bpp=BiocParallel::SerialParam()) {
register(bpp)
data(taggedPhenoDF)
data(lgenGWC_17)
GenomeInfoDb::seqlevelsStyle(lgenGWC_17) = "UCSC"
data(locErmaSet)
ct = cellTypes(locErmaSet)
data(short_celltype, package="erma")
st = short_celltype[ct]
ui = fluidPage(
  sidebarLayout(
   sidebarPanel(width=3,
    helpText("lungGen: annotation support for lung disease SNP"),
    selectInput("pheno", "phenotype", 
         choices=grep(regexpr, taggedPhenoDF$term, value=TRUE), selected = "asthma"),
    helpText("epigenomes available:"),
    selectInput("cellpicks", "cells", choices=st, selected=st[1], multiple=TRUE),
    uiOutput("genesel")
    ),
   mainPanel(
    tabsetPanel(
     tabPanel("pubmed",
      helpText("EBI/EMBL GWAS catalog entries 1/15/2019"),
      verbatimTextOutput("pms")
      ),
     tabPanel("loci",
      DT::dataTableOutput("snps")
      ),
     tabPanel("cell",
      verbatimTextOutput("cells")
      ),
     tabPanel("states",
      DT::dataTableOutput("states")
      ),
     tabPanel("txmodel",
      helpText("Gene model: red lines denote exons, each transcript vertically displaced.  Black points are 1+predicted probability of being within R^2 .8 of a GRASP SNP."),
      plotlyOutput("txplot")
      )
     )
    )
   )
  )
server = function(input, output) {
  colnames(locErmaSet) = short_celltype[cellTypes(locErmaSet)]
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
   print(head(dfr))
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
    preddf = data.frame(start=start(preds), pred=preds$pred+1, rsid=names(preds))
    preddf = preddf[preddf$pred > 1.03,] # FIXME
    ggplotly(vis + geom_point(data=preddf, 
        aes(x=start, y=pred, text=rsid)) + 
        theme(axis.text.y=element_blank()) + ylab(" "))
    })
  } # end server
runApp(list(ui=ui, server=server))
}

miameList = function(srchstring="asthma") {
     ids = unique(phenoToHits(srchstring, lgenGWC_17)$PUBMEDID)
     lapply(as.character(ids), annotate::pmid2MIAME)
}
