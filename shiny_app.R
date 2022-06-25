
test <- function () 
{
  refseq_kingdoms <- c("invertebrate")
  all_refseqOrgs <- as.vector(unlist(lapply(refseq_kingdoms, 
                                            function(kingdom) strsplit(RCurl::getURL(paste0("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/", 
                                                                                            kingdom, "/"), ftp.use.epsv = FALSE, dirlistonly = TRUE), 
                                                                       "\n"))))
  all_refseqOrgs <- stringr::str_replace(all_refseqOrgs, "_", 
                                         " ")
  return(all_refseqOrgs)
}


# Get available database
library(dplyr)
library(biomartr)
options <- getKingdoms(db = "refseq") 
refseq_table <- data.frame(listGenomes(db = "refseq", type = "kingdom", details = T))
kg_options <- unique(refseq_table$kingdoms)
refseq_sel <- refseq_table %>%
  filter(kingdoms == "Eukaryota")


# Rshiny able to load sample files seperately
library(shiny)
library(shinyWidgets)
library(DT)
library(shinythemes)

ui <- fluidPage(
  
)


server <- function(input, output, server) {
  
}


shinyApp(ui, server)

shiny::insertUI()
