#
# UI PREP
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
#data(locErmaSet)
library(GenomicFiles)
locErmaSet = makeEpig17()
cellTypes2 = function(x) colData(x)[,"Standardized.Epigenome.name"]
ct = cellTypes2(locErmaSet)
data(short_celltype, package="lungGen")
st = short_celltype[ct]
#
#
#

ui = fluidPage(
  sidebarLayout(
   sidebarPanel(width=3,
    helpText(h3("lungGen: annotation support for lung disease SNP")),
    selectInput("pheno", "phenotype", 
         choices=grep(regexpr, taggedPhenoDF$term, value=TRUE), selected = "asthma"),
    helpText(h3("epigenomes available:")),
    selectInput("cellpicks", "cells for epigenomic state exploration", choices=st, selected=st[1], multiple=FALSE),
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
     tabPanel("states",
      DT::dataTableOutput("states")
      ),
     tabPanel("txmodel",
      helpText("Gene model: red lines denote exons, each transcript vertically displaced.  Colored points are 1+predicted probability of being within R^2 .8 of a GRASP SNP.  Hover over for details.  Colors are based on ChromHMM states for selected cell type."),
      plotlyOutput("txplot")
      )
     )
    )
   )
  )


