library(shiny)
library(shinyjs)
library(shinythemes)
library(shinydashboard)
library(shinyWidgets)
library(shinybusy)
library(shinycssloaders)
library(shinycustomloader)
library(DT)
library(rjson)
library(GenomicFeatures)
library(GenomicAlignments)
library(GenomicRanges)
library(fastqcr)
library(Rfastp)
library(ngsReports)
library(QuasR)
library(Rhisat2)
library(rtracklayer)
library(biomartr)
library(Gviz)

fastqc_plot <- function(x, plot.type) {
  switch(plot.type,
         "Summary" = plotSummary(x),
         "Per Base Sequence Quality" = plotBaseQuals(x),
         "Mean Sequence Quality Per Read" = plotSeqQuals(x),
         "Per Base Sequence Content" = plotSeqContent(x),
         "GC Content" = plotGcContent(x),
         "N Content" = plotNContent(x),
         "Sequence Duplication Levels" = plotDupLevels(x),
         "Adapter Content" = plotAdapterContent(x),
         "Overrepresented Sequences" = plotOverrep(x))
}

multiqc.plot.type <- c("Summary", "Per Base Sequence Quality", "Mean Sequence Quality Per Read",
                       "Per Base Sequence Content", "GC Content", "N Content", "Sequence Duplication Levels",
                       "Adapter Content", "Overrepresented Sequences")

singleqc.plot.type <- c("Per Base Sequence Quality", "Mean Sequence Quality Per Read",
                        "Per Base Sequence Content", "GC Content", "N Content", "Sequence Duplication Levels",
                        "Adapter Content", "Overrepresented Sequences")


ensembl_table <- data.frame(getENSEMBLInfo())

align_plot <- function(qc_proj, type) {
  switch (type,
          "Quality Score" = QuasR:::plotQualByCycle(qc_proj$raw$qa),
          "Nucleotide Content by Cycle" = QuasR:::plotNuclByCycle(qc_proj$raw$qa),
          "Duplicate Level" = QuasR:::plotDuplicated(qc_proj$raw$qa),
          "Mapping Statistics" = QuasR:::plotMappings(qc_proj$raw$mapdata, a4layout = TRUE),
          "Library Complexity" = QuasR:::plotUniqueness(qc_proj$raw$unique, a4layout = TRUE),
          "Mismatch Frequency" = QuasR:::plotErrorsByCycle(qc_proj$raw$mm),
          "Mismatch Types" = QuasR:::plotMismatchTypes(qc_proj$raw$mm),
          "Fragment Size" = QuasR:::plotFragmentDistribution(qc_proj$raw$frag)
  )
}






