#' interface to txregnet resource for lung disease genomics
#' @import dplyr magrittr annotate gwascat shiny
#' @param regexpr character(1) will be used to grep taggedPhenoDF$term
#' @export
lungGen = function(regexpr="lung|asthma|pulmonary") {
data(taggedPhenoDF)
data(lgenGWC_17)
seqlevelsStyle(lgenGWC_17) = "UCSC"
data(locErmaSet)
ct = cellTypes(locErmaSet)
data(short_celltype)
st = short_celltype[ct]
ui = fluidPage(
  sidebarLayout(
   sidebarPanel(width=3,
    helpText("lungGen: annotation support for lung disease SNP"),
    selectInput("pheno", "phenotype", 
         choices=grep(regexpr, taggedPhenoDF$term, value=TRUE), selected = "asthma"),
    helpText("epigenomes available:"),
    checkboxGroupInput("cellpicks", "cells", choices=st, selected=st[1])
    ),
   mainPanel(
    tabsetPanel(
     tabPanel("pubmed",
      helpText("EBI/EMBL GWAS catalog entries 1/15/2019"),
      verbatimTextOutput("pms")
      ),
     tabPanel("loci",
      dataTableOutput("snps")
      ),
     tabPanel("cell",
      verbatimTextOutput("cells")
      ),
     tabPanel("states",
      dataTableOutput("states")
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
  output$snps = renderDataTable({
   hitgr = phenoToHits(input$pheno, lgenGWC_17)
   names(hitgr) = hitgr$SNPS
   dfr = as.data.frame(mcols(phenoToHits(input$pheno, lgenGWC_17))[,c("CHR_ID", "CHR_POS", "SNPS", "PUBMEDID")])
   chkdup = apply(dfr, 1, paste0, collapse=":")
   todr = which(duplicated(chkdup))
   if (length(todr)>0) dfr = dfr[-todr,]
   dfr$CHR_POS = as.integer(dfr$CHR_POS)
   dfr$PUBMEDID = paste0("<a href=https://www.ncbi.nlm.nih.gov/pubmed/?term=", dfr$PUBMEDID, " target='_blank'>PMID ", dfr$PUBMEDID, "</a>")
   o = order(dfr$CHR_POS)
   datatable(dfr[o,], escape=FALSE)
   })
  output$states = renderDataTable({
    xx = phenoToHits(input$pheno, lgenGWC_17)
    names(xx) = xx$SNPS
    states(xx, locErmaSet[, input$cellpicks])
    })
  }
runApp(list(ui=ui, server=server))
}

miameList = function(srchstring="asthma") {
     ids = unique(phenoToHits(srchstring, lgenGWC_17)$PUBMEDID)
     lapply(as.character(ids), annotate::pmid2MIAME)
}
