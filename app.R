################################# Download packages & dependencies #################################
if (length(find.package(package = "shiny", quiet = T)) > 0) {
  library(shiny)
} else {
  print("Package shiny not installed")
  install.packages("shiny")
  print("Package shiny installed")
  library(shiny)
}


if (length(find.package(package = "shinyWidgets", quiet = T)) > 0) {
  library(shinyWidgets)
} else {
  print("Package shinyWidgets not installed")
  install.packages("shinyWidgets")
  print("Package shinyWidgets installed")
  library(shinyWidgets)
}


if (length(find.package(package = "shinythemes", quiet = T)) > 0) {
  library(shinythemes)
} else {
  print("Package shinythemes not installed")
  install.packages("shinythemes")
  print("Package shinythemes installed")
  library(shinythemes)
}


if (length(find.package(package = "shinydashboard", quiet = T)) > 0) {
  library(shinydashboard)
} else {
  print("Package shinydashboard not installed")
  install.packages("shinydashboard")
  print("Package shinydashboard installed")
  library(shinydashboard)
}


if (length(find.package(package = "shinycustomloader", quiet = T)) > 0) {
  library(shinycustomloader)
} else {
  print("Package shinycustomloader not installed")
  install.packages("shinycustomloader")
  print("Package shinycustomloader installed")
  library(shinycustomloader)
}


if (length(find.package(package = "shinycssloaders", quiet = T)) > 0) {
  library(shinycssloaders)
} else {
  print("Package shinycssloaders not installed")
  install.packages("shinycssloaders")
  print("Package shinycssloaders installed")
  library(shinycssloaders)
}


if (length(find.package(package = "shinyBS", quiet = T)) > 0) {
  library(shinyBS)
} else {
  print("Package shinyBS not installed")
  install.packages("shinyBS")
  print("Package shinyBS installed")
  library(shinyBS)
}


if (length(find.package(package = "shinyvalidate", quiet = T)) > 0) {
  library(shinyvalidate)
} else {
  print("Package shinyvalidate not installed")
  install.packages("shinyvalidate")
  print("Package shinyvalidate installed")
  library(shinyvalidate)
}


if (length(find.package(package = "DT", quiet = T)) > 0) {
  library(DT)
} else {
  print("Package DT not installed")
  install.packages("DT")
  print("Package DT installed")
  library(DT)
}


if (length(find.package(package = "BiocManager", quiet = T)) > 0) {
  library(BiocManager)
} else {
  print("Package BiocManager not installed")
  install.packages("BiocManager")
  print("Package BiocManager installed")
  library(BiocManager)
}


if (length(find.package(package = "remotes", quiet = T)) > 0) {
  library(remotes)
} else {
  print("Package remotes not installed")
  BiocManager::install("remotes")
  print("Package remotes installed")
  library(remotes)
}


## Quality control: FASTQC
if (length(find.package(package = "fastqcr", quiet = T)) > 0) {
  library(fastqcr)
} else {
  print("Package fastqcr not installed")
  install.packages("fastqcr")
  print("Package fastqcr installed")
  library(fastqcr)
}


if (length(find.package(package = "ngsReports", quiet = T)) > 0) {
  library(ngsReports)
} else {
  print("Package ngsReports not installed")
  BiocManager::install("steveped/ngsReports")
  print("Package ngsReports installed")
  library(ngsReports)
}


## Quality control: Trimming
if (length(find.package(package = "threadr", quiet = T)) > 0) {
  library(threadr)
} else {
  print("Package threadr not installed")
  install.packages("threadr")
  print("Package threadr installed")
  library(threadr)
}


if (length(find.package(package = "Rfastp", quiet = T)) > 0) {
  library(Rfastp)
} else {
  print("Package Rfastp not installed")
  BiocManager::install("Rfastp")
  print("Package Rfastp installed")
  library(Rfastp)
}


if (length(find.package(package = "rjson", quiet = T)) > 0) {
  library(fastqcr)
} else {
  print("Package rjson not installed")
  install.packages("rjson")
  print("Package rjson installed")
  library(rjson)
}


## Align
if (length(find.package(package = "Rhisat2", quiet = T)) > 0) {
  library(Rhisat2)
} else {
  print("Package Rhisat2 not installed")
  BiocManager::install("Rhisat2")
  print("Package Rhisat2 installed")
  library(Rhisat2)
}


if (length(find.package(package = "Rbowtie", quiet = T)) > 0) {
  library(Rbowtie)
} else {
  print("Package Rbowtie not installed")
  BiocManager::install("Rbowtie")
  print("Package Rbowtie installed")
  library(Rbowtie)
}


if (length(find.package(package = "QuasR", quiet = T)) > 0) {
  library(QuasR)
} else {
  print("Package QuasR not installed")
  BiocManager::install("QuasR")
  print("Package QuasR installed")
  library(QuasR)
}


if (length(find.package(package = "Biostrings", quiet = T)) > 0) {
  library(Biostrings)
} else {
  print("Package Biostrings not installed")
  BiocManager::install("Biostrings")
  print("Package Biostrings installed")
  library(Biostrings)
}


if (length(find.package(package = "GenomicFeatures", quiet = T)) > 0) {
  library(GenomicFeatures)
} else {
  print("Package GenomicFeatures not installed")
  BiocManager::install("GenomicFeatures")
  print("Package GenomicFeatures installed")
  library(GenomicFeatures)
}


if (length(find.package(package = "threadr", quiet = T)) > 0) {
  library(threadr)
} else {
  print("Package threadr not installed")
  remotes::install_github("skgrange/threadr")
  print("Package threadr installed")
  library(threadr)
}


## Post alignment
if (length(find.package(package = "plotly", quiet = T)) > 0) {
  library(plotly)
} else {
  print("Package plotly not installed")
  install.packages("plotly")
  print("Package plotly installed")
  library(plotly)
}

if (length(find.package(package = "ramwas", quiet = T)) > 0) {
  library(ramwas)
} else {
  print("Package ramwas not installed")
  BiocManager::install("ramwas")
  print("Package ramwas installed")
  library(ramwas)
}


## Genome Browser
if (length(find.package(package = "Gviz", quiet = T)) > 0) {
  library(Gviz)
} else {
  print("Package Gviz not installed")
  BiocManager::install("Gviz")
  print("Package Gviz installed")
  library(Gviz)
}

options(shiny.maxRequestSize=10000*1024^2, ucscChromosomeNames=FALSE)
######################################################## Helper Functions ########################################################
## UI module: Sample upload
fastqUI <- function(id, label = "upload_fastq") {
  ns <- NS(id)
  tagList(
    wellPanel(
      shiny::tags$head(
        shiny::tags$style(type="text/css", "#inline label{ display: table-cell; text-align: center; vertical-align: middle; }
                #inline .form-group { display: table-row;}")
      ),

      shiny::tags$div(id = "inline", textInput(ns("fastq_sample"), "Sample Name:")),
      div(style = "margin-top: +10px"),
      conditionalPanel("input.read_type == 'single'",
                       fileInput(ns("fastq_se"), "Select Single Read FASTQ File"), accept = c(".zip")),
      conditionalPanel("input.read_type == 'paired'",
                       fileInput(ns("fastq_pe_r1"), "Select Read 1 FASTQ File", accept = c(".zip")),
                       div(style = "margin-top: -20px"),
                       fileInput(ns("fastq_pe_r2"), "Select Read 2 FASTQ File"), accept = c(".zip")),
      div(style = "margin-top: -10px"),
      shiny::tags$div(id = "inline", textInput(ns("fastq_suffix"), "FASTQ Suffix:")),
      div(style = "margin-top: +20px"),
      actionBttn(ns("unzip_fastq"), "Unzip FASTQ", icon = icon("file-import"), size = "sm")
    )
  )
}


## Server module: Sample upload
fastqServer <- function(id, read_type, folder_name) {
  moduleServer(
    id,

    function(input, output, session) {
      observeEvent(input$unzip_fastq, {
        if (read_type == "single") {
          unzip(input$fastq_se$datapath, exdir = paste0(folder_name, "/", input$fastq_sample))

          dir <- paste0(folder_name, "/", input$fastq_sample, "/", tools::file_path_sans_ext(input$fastq_se$name))
          FileName <- paste0(tools::file_path_sans_ext(input$fastq_se$name), "/", list.files(path = dir, full.names = F))
          SampleName <- gsub(pattern = input$fastq_suffix, replacement = "", x = basename(FileName))

          df <- data.frame(FileName, SampleName)
          names(df) <- c("FileName", "SampleName")

        } else {
          unzip(input$fastq_pe_r1$datapath, exdir = paste0(folder_name, "/", input$fastq_sample))
          unzip(input$fastq_pe_r2$datapath, exdir = paste0(folder_name, "/", input$fastq_sample))

          dir_r1 <- paste0(folder_name, "/", input$fastq_sample, "/", tools::file_path_sans_ext(input$fastq_pe_r1$name))
          FileName1 <- paste0(tools::file_path_sans_ext(input$fastq_pe_r1$name), "/", list.files(path = dir_r1, full.names = F))
          SampleName <- gsub(pattern = input$fastq_suffix, replacement = "", x = basename(FileName1))

          dir_r2 <- paste0(folder_name, "/", input$fastq_sample, "/", tools::file_path_sans_ext(input$fastq_pe_r2$name))
          FileName2 <- paste0(tools::file_path_sans_ext(input$fastq_pe_r2$name), "/", list.files(path = dir_r2, full.names = F))

          df <- data.frame(FileName1, FileName2, SampleName)
          names(df) <- c("FileName1", "FileName2", "SampleName")

        }

        write.table(df, file = paste0(folder_name, "/", input$fastq_sample, "/", "untrimmed.txt"), sep="\t", row.names=FALSE, quote = FALSE)

        updateActionButton(session, "unzip_fastq", label = "Unzipped!", icon = icon("check"))
      })
    }
  )
}


## FASTQC plots
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

## FASTQC plot type: Multiple Files
multiqc.plot.type <- c("Summary", "Per Base Sequence Quality", "Mean Sequence Quality Per Read",
                       "Per Base Sequence Content", "GC Content", "N Content", "Sequence Duplication Levels",
                       "Adapter Content", "Overrepresented Sequences")

## FASTQC plot type: Single File
singleqc.plot.type <- c("Per Base Sequence Quality", "Mean Sequence Quality Per Read",
                        "Per Base Sequence Content", "GC Content", "N Content", "Sequence Duplication Levels",
                        "Adapter Content", "Overrepresented Sequences")


## Annotation Files: NCBI Download
refgene_info <- read.delim2("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt", fill = T, header = F, sep = '\t', skip = 2)

names(refgene_info) <- c(
  "assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category",
  "taxid", "species_taxid", "organism_name", "infraspecific_name"	, "isolate",
  "version_status", "assembly_level", "release_type", "genome_rep", "seq_rel_date",
  "asm_name", "submitter", "gbrs_paired_asm", "paired_asm_comp", "ftp_path",
  "excluded_from_refseq", "relation_to_type_material"	, "asm_not_live_date"
)



############################################################# UI ############################################################# 
ui <- tagList(useShinydashboard(),
  navbarPage(id = "navbar", title = "",
    navbarMenu("Quality Control",

      ############################################### UI: FAST QC ###############################################
      tabPanel("FastQC",
        sidebarLayout(
          sidebarPanel(
            shiny::tags$style(".popover{ max-width: 100%;}"),
            radioButtons("read_type", "Choose Read Type", choices = c("single", "paired"), inline = T),
            bsPopover("read_type",
                      title = "Info: Fastq Suffix",
                      content = paste("Single-End: File Name (sample_1.fastq.gz) = Sample Name (sample_1) + Suffix (.fastq.gz)",
                                      "Paired-End: File Name (sample_1_R1.fastq.gz) = Sample Name (sample_1) + Suffix (_R1.fastq.gz)",
                                      sep = "</br>"),
                      placement = "right",
                      options=list(container="body")),
            actionButton('addSample', '', icon = icon('plus'))
          ),

          mainPanel(
            fluidRow(
              column(4, uiOutput("select_sample_qc")),
              column(4,
                     div(style = "margin-top: +25px"),
                     actionButton("runfastqc", "Run FASTQC")
              )
            ),

            tabsetPanel(
              tabPanel("About FASTQC",
                fluidRow(
                  h3(strong("Per Base Sequence Quality")),
                  column(width = 6,
                         HTML('<font color="#1ec73a"><h4><b>&#10004; GOOD</b></h4></font>'),
                         img(src = 'good_bq.png', align = "left", width = "300px", height = "300px")),
                  column(width = 6,
                         HTML('<font color="#de0914"><h4><b>&#10060; BAD</b></h4></font>'),
                         img(src = 'bad_bq.png', align = "left", width = "300px", height = "300px")),
                  h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                  HTML(
                    "<p>
                      <b>
                        <ul>
                          <li>This plot shows aggregated quality scores at each position along the read.</li>
                          <li>The blue line shows mean quality score while the red line shows median quality score at each base position/window.</li>
                          <li>Note that the average quality score will steadily drop over the length of the read, so it is common to see base calls<br>
                          falling into the orange area near towards the end.</li>
                          <li>With paired end reads the average quality scores for read 1 will almost always be higher than for read 2.</li>
                        </ul>
                      </b>
                    </p>"
                  ),

                  h4(span("What to look for:", style = "color:#73C6B6;font-weight: bold;")),
                  HTML(
                    "<p>
                      <b>
                        <ul>
                          <li>Good: Distribution of average read quality is fairly tight in the upper range.</li>
                          <li>Warning: Lower quartile for any base is less than 10, or if the median for any base is less than 25.</li>
                          <li>Failure: Lower quartile for any base is less than 5 or if the median for any base is less than 20.</li>
                        </ul>
                      </b>
                    </p>"
                  )
                ),
                
                div(style = "margin-top: +30px"),
                
                fluidRow(
                  h3(strong("Per Sequence Quality Scores")),
                  column(width = 6,
                         HTML('<font color="#1ec73a"><h4><b>&#10004; GOOD</b></h4></font>'),
                         img(src = 'good_sq.png', align = "left", width = "300px", height = "300px")),
                  column(width = 6,
                         HTML('<font color="#de0914"><h4><b>&#10060; BAD</b></h4></font>'),
                         img(src = 'bad_sq.png', align = "left", width = "300px", height = "300px")),
                  h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                  HTML(
                    "<p>
                      <b>
                        <ul>
                          <li>This plot shows the average quality score on the x-axis and the number of sequences with that average<br>
                          score on the y-axis. </li>
                          <li>Often, a subset of sequences will have universally poor quality, often because they are poorly imaged<br>
                          (eg. on the edge of the field of view).</li>
                        </ul>
                      </b>
                    </p>"
                  ),

                  h4(span("What to look for:", style = "color:#73C6B6;font-weight: bold;")),
                  HTML(
                    "<p>
                      <b>
                        <ul>
                          <li>Good: Majority of reads have a high average quality score with no large bumps at the lower quality values</li>
                          <li>Warning: Most frequently observed mean quality is below 27 - this equates to a 0.2% error rate.</li>
                          <li>Fail: Most frequently observed mean quality is below 20 - this equates to a 1% error rate.</li>
                        </ul>
                      </b>
                    </p>"
                  )
                ),
                
                div(style = "margin-top: +30px"),
                
                fluidRow(
                  h3(strong("Per Base Sequence Content")),
                  column(width = 6,
                         HTML('<font color="#1ec73a"><h4><b>&#10004; GOOD</b></h4></font>'),
                         img(src = 'good_sc.png', align = "left", width = "300px", height = "300px")),
                  column(width = 6,
                         HTML('<font color="#de0914"><h4><b>&#10060; BAD</b></h4></font>'),
                         img(src = 'bad_sc.png', align = "left", width = "300px", height = "300px")),
                  h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                  HTML(
                    "<p>
                      <b>
                        <ul>
                          <li>This plot shows the proportion of each base position in a file for which each of the four normal DNA bases has been called.</li>
                          <li>In a random library, there should be little to no difference between the bases in a run, so the lines should run parallel<br>
                          to each other.</li>
                          <li>If there is a significant and consistent difference in bases, this may indicate an overrepresented sequence which is <br>
                          contaminating the library or that there was a systematic problem during the sequencing of the library.</li>
                        </ul>
                      </b>
                    </p>"
                  ),

                  h4(span("What to look for:", style = "color:#73C6B6;font-weight: bold;")),
                  HTML(
                    "<p>
                      <b>
                        <ul>
                          <li>Good: 4 bases should remain relatively constant over the length of the read with %A=%T and %G=%C.</li>
                          <li>Warning: Difference between A and T, or G and C is greater than 10% in any position.</li>
                          <li>Failure: Difference between A and T, or G and C is greater than 20% in any position.</li>
                        </ul>
                      </b>
                    </p>"
                  )
                ),
                
                div(style = "margin-top: +30px"),
                
                fluidRow(
                  h3(strong("Per Sequence GC Content")),
                  column(width = 6,
                         HTML('<font color="#1ec73a"><h4><b>&#10004; GOOD</b></h4></font>'),
                         img(src = 'good_gc.png', align = "left", width = "300px", height = "300px")),
                  column(width = 6,
                         HTML('<font color="#de0914"><h4><b>&#10060; BAD</b></h4></font>'),
                         img(src = 'bad_gc.png', align = "left", width = "300px", height = "300px")),
                  h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                  HTML(
                    "<p>
                      <b>
                        <ul>
                          <li>This plots the GC content across the whole length of each sequence in a file and compares it to a modelled normal<br>
                          distribution of GC content.</li>
                          <li>Generally, the distribution should be normal unless there are over-represented sequences (sharp peaks on a normal<br>
                          distribution) or or contamination with another organism (broad peak).</li>
                          <li>A shifted normal distribution indicates some systematic bias which is independent of base position (will not be flagged).</li>
                          <li>Note: There are many situations where the observed distribution deviates too far from the theoretical, so a fail status<br>
                          can be ignored.</li>
                        </ul>
                      </b>
                    </p>"
                  ),

                  h4(span("What to look for:", style = "color:#73C6B6;font-weight: bold;")),
                  HTML(
                    "<p>
                      <b>
                        <ul>
                          <li>Good: GC content of all reads should form a normal distribution with the peak of the curve at the mean GC content for<br>
                          the organism sequenced</li>
                          <li>Warning: Sum of the deviations from the normal distribution represents more than 15% of the reads.</li>
                          <li>Fail:  Sum of the deviations from the normal distribution represents more than 30% of the reads.</li>
                        </ul>
                      </b>
                    </p>"
                  )
                ),
                
                div(style = "margin-top: +30px"),
                
                fluidRow(
                  h3(strong("Per base N content")),
                  column(width = 6,
                         HTML('<font color="#1ec73a"><h4><b>&#10004; GOOD</b></h4></font>'),
                         img(src = 'good_n.png', align = "left", width = "300px", height = "300px")),
                  column(width = 6,
                         HTML('<font color="#de0914"><h4><b>&#10060; BAD</b></h4></font>'),
                         img(src = 'bad_n.png', align = "left", width = "300px", height = "300px")),
                  h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                  HTML(
                    "<p>
                      <b>
                        Plots the percentage of bases at each position or bin with no base call, i.e. ‘N’.
                      </b>
                    </p>"
                  ),

                  h4(span("What to look for:", style = "color:#73C6B6;font-weight: bold;")),
                  HTML(
                    "<p>
                      <b>
                        <ul>
                          <li>Good: No point where this curve rises noticeably above zero.</li>
                          <li>Warning: Any position shows an N content of >5%.</li>
                          <li>Fail: Any position shows an N content of >20%.</li>
                        </ul>
                      </b>
                    </p>"
                  )
                ),
                
                div(style = "margin-top: +30px"),
                
                fluidRow(
                  h3(strong("Sequence Duplication Levels")),
                  column(width = 6,
                         HTML('<font color="#1ec73a"><h4><b>&#10004; GOOD</b></h4></font>'),
                         img(src = 'good_dup.png', align = "left", width = "300px", height = "300px")),
                  column(width = 6,
                         HTML('<font color="#de0914"><h4><b>&#10060; BAD</b></h4></font>'),
                         img(src = 'bad_dup.png', align = "left", width = "300px", height = "300px")),
                  h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                  HTML(
                    "<p>
                      <b>
                        <ul>
                          <li>This plots the percentage of reads of a given sequence and counts the degree of duplication for every sequence<br>
                          in the set.</li>
                          <li>Only sequences which occur in the first 200,000 sequences in each file are analysed and each sequence is tracked<br>
                          to the end of the file.</li>
                          <li>Any sequences with more than 10 duplicates are placed into the 10 duplicates category so it is not usual to see a<br>
                          small rise in this category.</li>
                          <li>Generally, there are two sources of duplicate reads: PCR duplication in which library fragments have been over<br>
                          represented due to biased PCR enrichment or truly over represented sequences such as very abundant transcripts in<br>
                          an RNA-Seq library.</li>
                          <li>PCR duplicates misrepresent the true proportion of sequences in your starting material and is a cause for concern.</li>
                        </ul>
                      </b>
                    </p>"
                  ),

                  h4(span("What to look for:", style = "color:#73C6B6;font-weight: bold;")),
                  HTML(
                    "<p>
                      <b>
                        <ul>
                          <li>Good: For whole genome shotgun data it is expected that nearly 100% of your reads will be unique (appearing only<br>
                           once), and indicates a highly diverse library that was not over sequenced.</li>
                          <li>Warning: Non-unique sequences make up more than 20% of the total.</li>
                          <li>Fail: Non-unique sequences make up more than 50% of the total.</li>
                        </ul>
                      </b>
                    </p>"
                  )
                ),
                
                div(style = "margin-top: +30px"),
                
                fluidRow(
                  h3(strong("Adapter Content")),
                  column(width = 6,
                         HTML('<font color="#1ec73a"><h4><b>&#10004; GOOD</b></h4></font>'),
                         img(src = 'good_ad.png', align = "left", width = "300px", height = "300px")),
                  column(width = 6,
                         HTML('<font color="#de0914"><h4><b>&#10060; BAD</b></h4></font>'),
                         img(src = 'bad_ad.png', align = "left", width = "300px", height = "300px")),
                  h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                  HTML(
                    "<p>
                      <b>
                        <ul>
                          <li>A cumulative plot of the fraction of reads where the sequence library adapter sequence is identified at the<br>
                          indicated base position.</li>
                          <li>Only adapters specific to the library type are searched.</li>
                          <li>Ideally Illumina sequence data should not have any adapter sequence present, however when using long read<br>
                          lengths it is possible that some of the library inserts are shorter than the read length resulting in read-through<br>
                          to the adapter at the 3’ end of the read.</li>
                        </ul>
                      </b>
                    </p>"
                  )
                ),
                
                div(style = "margin-top: +30px"),
                
                fluidRow(
                  h3(strong("Overrepresented Sequences")),
                  column(width = 12,
                         img(src = 'ovr.png', align = "left", width = "80%", height = "80%")),
                  h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                  HTML(
                    "<p>
                      <b>
                        <ul>
                          <li>This module lists all of the sequence which make up more than 0.1% of the total.</li>
                          <li>Only sequences which occur in the first 200,000 sequences in each file are analysed and each sequence is tracked<br>
                          to the end of the file.</li>
                          <li>As an exact sequence match is needed, any reads over 75bp in length are truncated to 50bp.</li>
                          <li>For each overrepresented sequence the program will look for matches in a database of common contaminants and<br>
                          will report the best hit it finds.</li>
                          <li>Note: Finding a hit to a contaminant may not be the exact source of contamination, but it should point you in the<br>
                          right direction.</li>
                        </ul>
                      </b>
                    </p>"
                  ),

                  h4(span("What to look for:", style = "color:#73C6B6;font-weight: bold;")),
                  HTML(
                    "<p>
                      <b>
                        <ul>
                          <li>Warning: Any sequence is found to represent more than 0.1% of the total.</li>
                          <li>Fail: Any sequence is found to represent more than 1% of the total.</li>
                        </ul>
                      </b>
                    </p>"
                  )
                )
                
              ),

              tabPanel("Multi-QC",
                fluidRow(
                  div(style = "margin-top: +20px"),
                  column(4, selectInput("multiqc_plot_type", "Type of Multi-QC Plot", choices = multiqc.plot.type)),
                ),

                withLoader(plotOutput("multiqc_plot"), type = "html", loader = "dnaspin")
              ),

              tabPanel("Single-QC",
                fluidRow(
                  div(style = "margin-top: +20px"),
                  column(4, uiOutput("select_qc_file")),
                  column(4, selectInput("singleqc_plot_type", "Type of Single-QC Plot", choices = singleqc.plot.type))
                ),

                withLoader(plotOutput("singleqc_plot"), type = "html", loader = "dnaspin")
              )
            )
          )
        )
      ),

      ############################################### UI: TRIMMING ###############################################
      tabPanel('Trimming',
        sidebarLayout(
          sidebarPanel(
            uiOutput("select_sample_tr"),

            h4(strong("Fastp Trim Settings")),
            div(style = "margin-top: +30px"),

            div(dataTableOutput("trim_settings"), style = "font-size:80%; height:400px; overflow-y: scroll"),
            div(style = "margin-top: +30px"),

            actionButton("runtrimming", "Run Trimming")
          ),

          mainPanel(
            tabsetPanel(
              tabPanel("Main Settings",
                dashboardSidebar(disable = TRUE),
                dashboardBody(
                  fluidRow(
                    box(width = 4,
                      title = "Adapter", status = "primary", solidHeader = TRUE,
                      materialSwitch("adapterTrimming", "Trim adaptor", status = 'success'),

                      conditionalPanel(condition = "input.adapterTrimming",

                                       radioButtons("adapter1", "Adapter detection (R1 for PE)", choices = c("auto", "manual"), inline = TRUE),

                                       conditionalPanel(condition = "input.adapter1 == 'manual'",
                                                        textInput("adapter1_manual", "Adapter Sequence (R1 for PE)")),

                                       conditionalPanel(condition = "input.read_type == 'paired'",
                                                        radioButtons("adapter2", "Adapter detection (R2)", choices = c("auto", "manual"), inline = TRUE),
                                                        conditionalPanel(condition = "input.adapter2 == 'manual'",
                                                                         textInput("adapter2_manual", "Adapter Sequence (R2)", value = ""),
                                                                         div(style = "margin-top: -15px"),
                                                                         helpText(HTML("<font size=-1>Leave blank if adapter sequence for R1 & R2 are the same.</font>")))))
                    ),

                    box(width = 4,
                      title = "Quality", status = "primary", solidHeader = TRUE,

                      materialSwitch("qualityFiltering", "Filter by quality", status = "success"),

                      conditionalPanel(condition = "input.qualityFiltering",
                                       numericInput("qualityFilterPhred", "Min base quality score", value = 15),
                                       div(style = "margin-top: -15px"),
                                       helpText(HTML("<font size=-1>Min quality phred score that a base is qualified.</font>")),

                                       div(style = "margin-top: +20px"),

                                       sliderInput("qualityFilterPercent", "Max % of bases < Min quality", min = 0, max = 100, value = 40, post = '%'),
                                       div(style = "margin-top: -15px"),
                                       helpText(HTML("<font size=-1>
                                                        Max percentage of bases that is allowed to be unqualified (0-100).
                                                        Default 40 means 40%.
                                                     </font>")),

                                       div(style = "margin-top: +20px"),

                                       numericInput("averageQualFilter", "Avg read quality", value = 0),
                                       div(style = "margin-top: -15px"),
                                       helpText(HTML("<font size=-1>
                                                        If a read's average quality score is less than this value, then this
                                                        read/pair is discarded.
                                                     </font>"))),

                      div(style = "margin-top: +20px"),

                      numericInput("maxNfilter", "Max N", value = 5),
                      div(style = "margin-top: -15px"),
                      helpText(HTML("<font size=-1>Max no. of N allowed in the sequence, else discarded.</font>"))

                    ),

                    tabBox(width = 4,
                      title = "Per Read Cutting", side = "right", selected = "Enable",

                      tabPanel("Settings", value = "Settings",
                               conditionalPanel(condition = "input.cutLowQualTail",
                                                sliderInput("cutTailWindowSize", "Window size: Cut tail", min = 1, max=1000, value = 4),
                                                numericInput("cutTailMeanQual", "Mean quality: Cut tail", value = 20)),

                               conditionalPanel(condition = "input.cutSlideWindowRight",
                                                sliderInput("cutSlideWindowSize", "Window size: Cut right", min = 1, max=1000, value = 4),
                                                numericInput("cutSlideWindowQual", "Mean quality: Cut right", value = 20)),

                               conditionalPanel(condition = "input.cutLowQualFront",
                                                sliderInput("cutFrontWindowSize", "Window size: Cut front", min = 1, max=1000, value = 4),
                                                numericInput("cutFrontMeanQual", "Mean quality: Cut front", value = 20))
                      ),

                      tabPanel("Enable", value = "Enable",
                        materialSwitch("cutLowQualTail", "Cut tail", status = "success"),
                        div(style = "margin-top: -15px"),
                        helpText(HTML("<font size=-1>
                                        Move a sliding window from tail (3') to the front, drop the bases in the window if its
                                        mean quality < threshold, stop otherwise. This is similar to Trimmomatic TRAILING method.
                                      </font>")),

                        div(style = "margin-top: +20px"),

                        materialSwitch("cutSlideWindowRight", "Cut right", status = "success"),
                        div(style = "margin-top: -15px"),
                        helpText(HTML("<font size=-1>
                                        Move a sliding window from front to tail, if meet one window with mean quality < threshold,
                                        drop the bases in the window and the right part, and then stop. This is simliar to Trimmomatic
                                        SLIDINGWINDOW method.
                                      </font>")),

                        div(style = "margin-top: +20px"),

                        materialSwitch("cutLowQualFront", "Cut front", status = "success"),
                        div(style = "margin-top: -15px"),
                        helpText(HTML("<font size=-1>
                                        Move a sliding window from front (5') to the tail, drop the bases in the window if its
                                        mean quality < threshold, stop otherwise. This is similar to Trimmomatic LEADING method.
                                      </font>"))
                      )
                    )
                  )
                )
              ),

              tabPanel("Advanced Settings",
                dashboardSidebar(disable = TRUE),
                dashboardBody(
                  fluidRow(
                    box(width = 4,
                        title = "Length", status = "success", solidHeader = TRUE,

                        materialSwitch("lengthFiltering", "Filter by length", status = "success"),
                        div(style = "margin-top: -15px"),
                        helpText(HTML("<font size=-1>Note that length filtering is applied as a last step.</font>")),

                        div(style = "margin-top: +20px"),

                        conditionalPanel(condition = "input.lengthFiltering",
                                         numericInput("minReadLength", "Min read length", value = 15),
                                         div(style = "margin-top: -15px"),
                                         helpText(HTML("<font size=-1>Reads shorter than min length will be discarded.</font>")),

                                         div(style = "margin-top: +20px"),

                                         numericInput("maxReadLength", "Max read length", value = 0),
                                         div(style = "margin-top: -15px"),
                                         helpText(HTML("<font size=-1>Reads shorter than min length will be discarded. 0 means no limitation.</font>")))
                    ),

                    box(width = 4,
                      title = "Global Trimming", status = "success", solidHeader = TRUE,

                      numericInput("trimFrontRead1", "Trim front R1", value = 0),
                      numericInput("trimTailRead1", "Trim tail R1", value = 0),

                      conditionalPanel(condition = "input.read_type == 'paired'",
                                       numericInput("trimFrontRead2", "Trim front R2", value = 0),
                                       numericInput("trimTailRead2", "Trim tail R2", value = 0)),
                      div(style = "margin-top: -15px"),
                      helpText(HTML("<font size=-1>No. of bases to hard trim from the front/tail.</font>")),

                      div(style = "margin-top: +20px"),

                      numericInput("maxLengthRead1", "Max length (R1 for PE)", value = 0),

                      conditionalPanel(condition = "input.read_type == 'paired'",
                                       numericInput("maxLengthRead2", "Max length R2", value = 0)),
                      div(style = "margin-top: -15px"),
                      helpText(HTML("<font size=-1>
                                      0 means no limitation. Note that the max length limitation wil only be applied as
                                      the last step.
                                    </font>"))
                    ),

                    box(width = 4,
                      title = "Tail Trimming", status = "success", solidHeader = TRUE,

                      materialSwitch("forceTrimPolyG", "Force trim polyG tail", status = "success"),
                      div(style = "margin-top: -15px"),
                      helpText(HTML("<font size=-1>
                                      <li>For Illumina NextSeq/NovaSeq data, polyG can happen in read tails since G means no signal in the Illumina two-color systems.</li>
                                      <li>NextSeq/NovaSeq data is detected by the machine ID in the FASTQ records and this feature is automatically enabled for such data.</li>
                                    </font>")),

                      div(style = "margin-top: +20px"),

                      materialSwitch("disableTrimPolyG", "Disable polyG tail trimming", status = "success"),
                      div(style = "margin-top: -15px"),
                      helpText(HTML("<font size=-1>Disable automatic polyG trimming for Illumina NextSeq/NovaSeq data.</font>")),

                      div(style = "margin-top: +20px"),

                      numericInput("minLengthPolyG", "Min polyG tail length", value = 10),
                      div(style = "margin-top: -15px"),
                      helpText(HTML("<font size=-1>Min length to detect polyG</font>")),

                      div(style = "margin-top: +20px"),

                      materialSwitch("trimPolyX",  "Force trim polyX tail", status = "success"),
                      numericInput("minLengthPolyX", "Min polyX tail length", value = 10),
                      div(style = "margin-top: -15px"),
                      helpText(HTML("<font size=-1>Min length to detect polyX</font>"))
                    )
                  )
                )
              ),

              tabPanel("Other Settings",
                dashboardSidebar(disable = TRUE),
                dashboardBody(
                  fluidRow(

                    box(width = 4,
                      title = "Overrepresentation Analysis", status = "warning", solidHeader = TRUE,

                      materialSwitch("overrepresentationAnalysis", "Enable overrepresentation analysis", status = "success"),

                      conditionalPanel(condition = "input.overrepresentationAnalysis",
                                       numericInput("overrepresentationSampling", "No. of reads for analysis", value = 20),
                                       div(style = "margin-top: -15px"),
                                       helpText(HTML("<font size=-1>1 in n reads is used for sequence counting. Note that if n=1, processing will be extremely slow.</font>")))
                    ),

                    box(width = 4,
                      title = "Base Correction", status = "warning", solidHeader = TRUE,

                      conditionalPanel(condition = "input.read_type == 'single'",
                                       helpText("Only for PE reads"),
                                       div(style = "margin-top: +20px")),

                      conditionalPanel(condition = "input.read_type == 'paired'",
                                       materialSwitch("correctionOverlap", "Run correction overlap", status = "success"),
                                       div(style = "margin-top: -15px"),
                                       helpText(HTML('<font size=-1>
                                                        If a proper overlap between a pair of reads is found, it can correct mismatched base pairs or if 1 base is with
                                                        high quality while the other is with ultra low quality. If a base is corrected, the quality of its paired base
                                                        will be assigned to it so that they will share the same quality.
                                                     </font>')),

                                       div(style = "margin-top: +20px"),

                                       conditionalPanel(condition = "input.correctionOverlap",
                                                        numericInput("minOverlapLength", "Min overlap length", value = 30),
                                                        div(style = "margin-top: -15px"),
                                                        helpText(HTML("<font size=-1>Min length to detect overlapped region of PE reads.</font>")),

                                                        div(style = "margin-top: +20px"),

                                                        numericInput("maxOverlapMismatch", "Max mismatched bases", value = 5),
                                                        div(style = "margin-top: -15px"),
                                                        helpText(HTML("<font size=-1>Max number of mismatched bases to detect overlapped region of PE reads.</font>")),

                                                        div(style = "margin-top: +20px"),

                                                        sliderInput("maxOverlapMismatchPercentage", "Max mismatch % in overlap", min = 0, max = 100, value = 20, post = '%'),
                                                        div(style = "margin-top: -15px"),
                                                        helpText(HTML("<font size=-1>Max percentage of mismatched bases to detect overlapped region of PE reads.</font>"))))

                    ),

                    box(width = 4,
                      title = "UMI", status = "warning", solidHeader = TRUE,
                      materialSwitch("umi", "Preprocess UMI", status = "success"),

                      conditionalPanel(condition = "input.umi",
                                       radioButtons("umiLoc", "Location of UMI", choices = c("read1", "read2", "per_read"), inline = T),
                                       div(style = "margin-top: -15px"),
                                       helpText(HTML("<font size=-1>
                                                        <li>read1: Head of read1 is used as UMI. If PE, used for both reads.</li>
                                                        <li>read2: Head of read2 is used as UMI. If PE, used for both reads.</li>
                                                        <li>For single-end reads, select read1</li>
                                                     </font>")),

                                       div(style = "margin-top: +20px"),

                                       numericInput("umiLength", "Length of UMI in read1/read2", value = 0),

                                       numericInput("umiSkipBaseLength", "No. of bases to skip", value = 0),
                                       div(style = "margin-top: -15px"),
                                       helpText(HTML("<font size=-1>Skip n bases after UMI to trim the UMI separator and A/T tailing in read1/read2.</font>")),

                                       div(style = "margin-top: +20px"),

                                       textInput("umiPrefix", "UMI Prefix", value = ""),
                                       div(style = "margin-top: -15px"),
                                       helpText(HTML("<font size=-1>If prefix is specified, an underline will be used to connect it and UMI.</br>
                                                        For example, UMI=AATTCCGG, prefix=UMI, then the final string presented in the name will be UMI_AATTCCGG.
                                                     </font>")),

                                       div(style = "margin-top: +20px"),

                                       materialSwitch("umiNoConnection", "Remove '_' between UMI prefix & UMI string", status = "success"))
                    )
                  )
                )
              ),

              ############################################### UI: TRIM RESULTS ###############################################
              tabPanel("Results",
                div(style = "margin-top: +20px"),

                fluidRow(
                  column(4, uiOutput("select_trim_fastq")),
                  column(4,
                         div(style = "margin-top: +25px"),
                         downloadButton("trim_report", "Download Trim Report")
                  ),
                ),

                h4(span("Summary:", style = "color:#73C6B6;font-weight: bold;")),
                withLoader(tableOutput("table_trim_summary"), type = "html", loader = "loader6"),

                div(style = "margin-top: +20px"),

                h4(span("Base Quality Comparison:", style = "color:#73C6B6;font-weight: bold;")),
                withLoader(plotOutput("trim_base_quality"), type = "html", loader = "dnaspin"),

                div(style = "margin-top: +20px"),

                h4(span("Base Content Comparison:", style = "color:#73C6B6;font-weight: bold;")),
                withLoader(plotOutput("trim_gc_content"), type = "html", loader = "dnaspin"),

                div(style = "margin-top: +20px"),

                h4(span("Over-Represented Sequences:", style = "color:#73C6B6;font-weight: bold;")),
                div(withLoader(dataTableOutput("trim_overrepresented"), type = "html", loader = "dnaspin"),
                    style = "overflow-y: scroll; overflow-x: scroll")
              )
            )
          )
        )
      )
    ),

    navbarMenu("Generate Counts",

      ############################################### UI: ALIGN ###############################################
      tabPanel("Align",
        sidebarLayout(
          sidebarPanel(width = 5,
            tabsetPanel(
              tabPanel("Upload",
                div(style = "margin-top: +20px"),
                materialSwitch("upload_how", "Upload Annotation Files?", status = "success"),

                conditionalPanel(condition = "input.upload_how",
                                 fileInput("upload_fasta", "Upload FASTA", accept = c(".fasta", ".fa", ".fna")),
                                 fileInput("upload_gtf", "Upload Gene Annotation File", accept = c(".gtf", ".gff", "gff3", ".gz"))),

                conditionalPanel(condition = "!input.upload_how",
                                 textInput("ncbi_accession", "Enter NCBI Accession No."),
                                 uiOutput("select_ncbi_fasta"),
                                 uiOutput("select_ncbi_gtf")),

                div(style = "margin-top: +20px"),

                radioButtons("gtf_format", "Select Annotation File Format", choices = c("gtf", "gff", "gff3"), inline = T),
                actionBttn("prepare_ann", "Prepare Annotation Files", icon = icon("file-import"), size = "sm")

              ),

              tabPanel("Settings",
                div(style = "margin-top: +20px"),

                uiOutput("select_sample_al"),
                radioButtons("trimmed_untrimmed", "Use Trimmed/Untrimmed Files?", choices = c("Trimmed", "Untrimmed"), inline = T),
                bsPopover("trimmed_untrimmed",
                          title = "Info",
                          content = paste("Trimmed: Used trimmed fastq files from previous trimming function (Rfastp)",
                                          "Untrimmed: Used original uploaded fastq files (that could have been trimmed by external application)",
                                          sep = "</br>"),
                          placement = "right",
                          options=list(container="body")),

                div(style = "margin-top: +20px"),

                conditionalPanel("input.read_type == 'paired'",
                                  radioButtons("paired", "Read Type", choices = c("fr", "ff", "rf"), inline = TRUE),
                                  bsPopover("paired",
                                             title = "Info",
                                             content = paste("fr: (Paired-End) Forward/Reverse",
                                                             "ff: (Paired-End) Foward/Foward",
                                                             "rf: (Paired-End) Reverse/Foward",
                                                             sep = "</br>"),
                                             placement = "right",
                                             options=list(container="body")),

                                 div(style = "margin-top: +20px"),
                ),

                materialSwitch("splicedAlignment", strong("Spiced Alignment?"), status = "success"),
                radioButtons("aligner", "Select Aligner", choices = c("Rbowtie", "Rhisat2"), inline = TRUE),
                bsPopover("aligner",
                          title = "Info",
                          content = paste("Rhisat2 is the recommended setting for spliced alignments.",
                                          "If Rbowtie is selected for spliced alignments, SpliceMap will be used.",
                                          "SpliceMap is about ten-fold slower, less sensitive than Hisat2 and can only be used for",
                                          "reads with a minimal length of 50nt (shorter reads are neither mapped or unmapped).",
                                          sep = "</br>"),
                          placement = "right",
                          options=list(container="body")),

                div(style = "margin-top: +20px"),

                numericInput("maxHits", "Max Mapping Positions Per Read", value = 1, min = 1),
                bsPopover("maxHits",
                          title = "Info",
                          content = paste("Maximal number of allowed mapping positions per read",
                                          "If a read produces more than maxHits alignments, no alignments will be reported for it.",
                                          "In case of a multi-mapping read, a single alignment is randomly selected.",
                                          sep = "</br>"),
                          placement = "right",
                          options=list(container="body")),

                div(style = "margin-top: +20px"),

                radioButtons("selectReadPosition", "Read Position", choices = c("start", "end"), inline = TRUE),
                bsPopover(id = "selectReadPosition",
                          title = "Info",
                          content = paste("The part of the alignment that has to be contained within a query region to produce an overlap.",
                                          "Start: Start of the alignment",
                                          "End: End of the alignment",
                                          sep = "</br>"),
                          placement = "right",
                          options=list(container="body")),

                div(style = "margin-top: +20px"),

                radioButtons("orientation", "Orientation of Alignments", choices = c("any", "same", "opposite"), inline = TRUE),
                bsPopover(id = "orientation",
                          title = "Info",
                          content = paste("The required orientation of the alignments relative to the query region in order to be counted.",
                                          "any: Count alignment on same & opp strand",
                                          "same: Count alignment only on same strand",
                                          "opposite: Count alignment only on opp strand",
                                          sep = "</br>"),
                          placement = "right",
                          options=list(container="body")),

                div(style = "margin-top: +20px"),

                radioButtons("useRead", "Read Mate to Count Alignments", choices = c("any", "first", "last"), inline = TRUE),
                bsPopover(id = "useRead",
                          title = "Info",
                          content = paste("The read mate whose alignments should be counted for paired-end reads.",
                                          "any: count all alignments",
                                          "first: count only alignments from the first read",
                                          "last: count only alignments from the last read",
                                          sep = "</br>"),
                          placement = "right",
                          options=list(container="body")),

                div(style = "margin-top: +20px"),

                materialSwitch("includeSpliced", strong("Count Spliced Alignments?"), status = "success"),
                bsPopover(id = "includeSpliced",
                          title = "Info",
                          content = paste("A spliced alignment is defined as an alignment with a gap in the read",
                                          "of at least 60 bases.",
                                          sep = "</br>"),
                          placement = "right",
                          options=list(container="body")),

                div(style = "margin-top: +20px"),

                materialSwitch("includeSecondary", strong("Count Secondary Alignments?"), status = "success"),

                div(style = "margin-top: +20px"),

                sliderInput("mapqlty", "Mapping Quality of Alignments When Counting", min = 0, max = 255, value = c(0, 255), step = 1),
                bsPopover(id = "mapqlty",
                          title = "Info",
                          content = paste("Min/Max mapping quality of alignments to be included when counting.",
                                          "Min = 0, Max = 255: Include all alignments",
                                          sep = "</br>"),
                          placement = "right",
                          options=list(container="body")),

                div(style = "margin-top: +20px"),

                actionButton("runalign", "Run Aligning")

              ),
            )
          ),

          mainPanel(width = 7,
            tabsetPanel(
              tabPanel("NCBI Database",
                div(style = "margin-top: +20px"),
                dataTableOutput("refgene_table")
              ),

              tabPanel("Counts",
                div(style = "margin-top: +20px"),

                h3(span("Gene Counts", style = "color:#73C6B6;font-weight: bold;")),
                downloadButton("gene_counts", "Download: Gene Counts", size = "xs"),
                h4("Preview:"),
                withLoader(dataTableOutput("gene_preview"), type = "html", loader = 'dnaspin'),
                div(style = "margin-top: +30px"),

                h3(span("Exon Counts", style = "color:#73C6B6;font-weight: bold;")),
                downloadButton("exon_counts", "Download: Exon Counts", size = "xs"),
                h4("Preview:"),
                withLoader(dataTableOutput("exon_preview"), type = "html", loader = 'dnaspin'),
                div(style = "margin-top: +30px"),

                h3(span("Promoter Counts", style = "color:#73C6B6;font-weight: bold;")),
                downloadButton("prm_counts", "Download: Promoter Counts", size = "xs"),
                h4("Preview:"),
                withLoader(dataTableOutput("prm_preview"), type = "html", loader = 'dnaspin'),
                div(style = "margin-top: +30px"),

                conditionalPanel(condition = "input.splicedAlignment",
                                 h3(span("Junction Counts", style = "color:#73C6B6;font-weight: bold;")),
                                 downloadButton("jun_counts", "Download: Junction Counts", size = "xs"),
                                 h4("Preview:"),
                                 withLoader(dataTableOutput("jun_preview"), type = "html", loader = 'dnaspin'))

              )
            )
          )
        )
      ),

      ############################################### UI: ALIGNMENT STATISTICS ###############################################
      tabPanel("Post-Alignment",
        sidebarLayout(
          sidebarPanel(
            uiOutput("select_sample_as"),
            downloadButton("alignment_report", "Download Alignment Report"),

            div(style = "margin-top: +20px"),

            sliderInput("minscore", "Min Alignment Score (MAPQ)", min = 0, max = 255, value = 0),
            numericInput("minfragmentsize", "Min Fragment Size", value = 50),
            numericInput("maxfragmentsize", "Max Fragment Size", value = 250),
            numericInput("maxrepeats", "Max Repeats", value = 1),
            actionButton("view_statistics", "View Bam Statistics")
          ),

          mainPanel(
            fluidRow(
              column(width = 4,
                     h6("No. of reads aligned to the reference genome"),
                     div(style = "margin-top: -10px"),
                     withLoader(plotlyOutput("fig1", width = 250, height = 250), type = "html", loader = "loader6")),

              column(width = 4,
                     h6("No. of recorded reads (passed min score)"),
                     div(style = "margin-top: -10px"),
                     withLoader(plotlyOutput("fig2", width = 250, height = 250), type = "html", loader = "loader6")),

              column(width = 4,
                     h6("No. of reads after removal of duplicate reads"),
                     div(style = "margin-top: -10px"),
                     withLoader(plotlyOutput("fig3", width = 250, height = 250), type = "html", loader = "loader6"))
            ),

            fluidRow(
              column(width = 4, offset = 2,
                     h6("No. of recorded reads aligned to each strand"),
                     div(style = "margin-top: -10px"),
                     withLoader(plotlyOutput("fig4", width = 250, height = 250), type = "html", loader = "loader6")),

              column(width = 4,
                     h6("No. of recorded reads aligned to each strand (exclude repeated reads)"),
                     div(style = "margin-top: -10px"),
                     withLoader(plotlyOutput("fig5", width = 250, height = 250), type = "html", loader = "loader6"))
            ),

            box(width = NULL,
                title = "Alignment Scores",
                verbatimTextOutput("fig6_summary"),
                withLoader(plotOutput("fig6"), type = "html", loader = "dnaspin")),

            box(width = NULL,
                title = "Alignment Length",
                verbatimTextOutput("fig7_summary"),
                withLoader(plotOutput("fig7"), type = "html", loader = "dnaspin")),

            box(width = NULL,
                title = "Edit Distance",
                verbatimTextOutput("fig8_summary"),
                withLoader(plotOutput("fig8"), type = "html", loader = "dnaspin"))
          )
        )
      ),
      
      ############################################### UI: GENOME BROWSER ###############################################
      tabPanel("Genome Browser",
        sidebarLayout(
          sidebarPanel(
            uiOutput("select_sample_gb"),
            uiOutput("select_chrom"),
            uiOutput("select_chrom_range"),

            div(style = "margin-top: +30px"),

            numericInput("extend_right", "Extend plot from the right: ", value = 0),
            numericInput("extend_left", "Extend plot from the left:", value = 0),
            actionButton("view_cov", "View Coverage")
          ),

          mainPanel(
            tabsetPanel(
              tabPanel("View Ranges",
                div(style = "margin-top: +30px"),
                radioButtons("sel_ranges_type", "Genomic Feature Type", choices = c("gene", "transcript", "exon"), inline = TRUE),

                div(style = "margin-top: +20px"),
                
                withLoader(DTOutput("selected_ranges"), type = "html", loader = "dnaspin")
              ),
              
              tabPanel('View Coverage',
                withLoader(plotOutput("cov_plot", width = "100%", height = "1200px"), type = "html", loader = "dnaspin")
              )
            )
          )
        )
      )

    )
  )
)



############################################################# SERVER ############################################################# 
server <- function(input, output, session) {
  #dir.create(session$token)
  ############################################### SERVER: FASTQ Upload ###############################################
  observeEvent(input$addSample, {
    id <- sprintf('fastqSample%s', input$addSample)
    insertUI(
      selector = "#addSample",
      where = "beforeBegin",
      ui = fastqUI(id)
    )
    
    fastqServer(id, input$read_type, session$token)
  })
  
  allsamples <- reactive({
    unlist(lapply(seq_len(input$addSample), function(i) input[[paste0("fastqSample", i, "-fastq_sample")]]))
  })
  
  raw_se_dir <- reactive({
    unlist(lapply(seq_len(input$addSample), function(i) tools::file_path_sans_ext(input[[paste0("fastqSample", i, "-fastq_se")]][["name"]])))
  })
  
  raw_per1_dir <- reactive({
    unlist(lapply(seq_len(input$addSample), function(i) tools::file_path_sans_ext(input[[paste0("fastqSample", i, "-fastq_pe_r1")]][["name"]])))
  })
  
  raw_per2_dir <- reactive({
    unlist(lapply(seq_len(input$addSample), function(i) tools::file_path_sans_ext(input[[paste0("fastqSample", i, "-fastq_pe_r2")]][["name"]])))
  })
  
  fastq_suffix <- reactive({
    unlist(lapply(seq_len(input$addSample), function(i) input[[paste0("fastqSample", i, "-fastq_suffix")]]))
  })
  
  
  
  ############################################### SERVER: FAST QC ###############################################
  # Select FASTQ sample to run FASTQC
  output$select_sample_qc <- renderUI({
    selectInput("fastqSample_qc", "Select FASTQ Sample", choices = allsamples())
  })
  
  # Download FASTQC & save results in fastqc folder inside FASTQ sample folder
  observeEvent(input$runfastqc, {
    showModal(modalDialog("Loading FASTQC ... ", footer=NULL))
    
    index <- match(input$fastqSample_qc, allsamples())
    sample_dir <- paste0(session$token, "/", input$fastqSample_qc)
    
    if (!dir.exists('FastQC')) {
      fastqcr::fastqc_install(dest.dir = getwd())
    }
    
    if (input$read_type == "single") {
      fastqc(fq.dir = paste0(sample_dir, "/", raw_se_dir()[index]),
             qc.dir = paste0(sample_dir, "/", "fastqc"),
             fastqc.path = "FastQC/fastqc")
    } else {
      fastqc(fq.dir = paste0(sample_dir, "/", raw_per1_dir()[index]),
             qc.dir = paste0(sample_dir, "/", "fastqc"),
             fastqc.path = "FastQC/fastqc")
      
      fastqc(fq.dir = paste0(sample_dir, "/", raw_per2_dir()[index]),
             qc.dir = paste0(sample_dir, "/", "fastqc"),
             fastqc.path = "FastQC/fastqc")
    }
    
    removeModal()
  })
  
  # List fastqc results for multi-plot & single-plot
  fastqc_files <- reactive({
    req(input$runfastqc)
    fastqc_dir <- paste0(session$token, "/", input$fastqSample_qc, "/", "fastqc")
    list.files(fastqc_dir, pattern = "fastqc.zip$", full.names = TRUE)
  })
  
  fdl <- reactive({FastqcDataList(fastqc_files())})
  
  output$multiqc_plot <- renderPlot({
    req(input$multiqc_plot_type)
    fastqc_plot(fdl(), input$multiqc_plot_type)
  })
  
  # Select fastqc result from sample for single-plot
  output$select_qc_file <- renderUI({
    selectInput("raw_fastq_file", "Select FASTQC File", choices = basename(fastqc_files()))
  })
  
  output$singleqc_plot <- renderPlot({
    req(input$raw_fastq_file, input$singleqc_plot_type)
    index <- match(input$raw_fastq_file, basename(fastqc_files()))
    which_fdl <- fdl()[[fastqc_files()[index]]]
    fastqc_plot(which_fdl, input$singleqc_plot_type)
  })
  
  
  ############################################### SERVER: TRIM SETTINGS ###############################################
  # Gather all trim settings
  read_type <- reactive({input$read_type})
  
  adapterTrimming <- reactive({input$adapterTrimming})
  adapterSequenceRead1 <- reactive({switch(input$adapter1, "auto"="auto", "manual"=toupper(input$adapter1_manual))})
  adapterSequenceRead2 <- reactive({
    if (input$adapter2 != "auto") {
      if (input$adapter2_manual == "") {
        adapterSequenceRead1()
      } else {
        toupper(input$adapter2_manual)
      }
    } else {
      "auto"
    }
  })
  
  qualityFiltering <- reactive({input$qualityFiltering})
  qualityFilterPhred <- reactive({input$qualityFilterPhred})
  qualityFilterPercent <- reactive({input$qualityFilterPercent})
  averageQualFilter <- reactive({input$averageQualFilter})
  maxNfilter <- reactive({input$maxNfilter})
  
  cutLowQualTail <- reactive({input$cutLowQualTail})
  cutTailWindowSize <- reactive({input$cutTailWindowSize})
  cutTailMeanQual <- reactive({input$cutTailMeanQual})
  
  cutSlideWindowRight <- reactive({input$cutSlideWindowRight})
  cutSlideWindowSize <- reactive({input$cutSlideWindowSize})
  cutSlideWindowQual <- reactive({input$cutSlideWindowQual})
  
  cutLowQualFront <- reactive({input$cutLowQualFront})
  cutFrontWindowSize <- reactive({input$cutFrontWindowSize})
  cutFrontMeanQual <- reactive({input$cutFrontMeanQual})
  
  lengthFiltering <- reactive({input$lengthFiltering})
  minReadLength <- reactive({input$minReadLength})
  maxReadLength <- reactive({input$maxReadLength})
  
  trimFrontRead1 <- reactive({input$trimFrontRead1})
  trimTailRead1 <- reactive({input$trimTailRead1})
  trimFrontRead2 <- reactive({input$trimFrontRead2})
  trimTailRead2 <- reactive({input$trimTailRead2})
  maxLengthRead1 <- reactive({input$maxLengthRead1})
  maxLengthRead2 <- reactive({input$maxLengthRead2})
  
  forceTrimPolyG <- reactive({input$forceTrimPolyG})
  disableTrimPolyG <- reactive({input$disableTrimPolyG})
  minLengthPolyG <- reactive({input$minLengthPolyG})
  trimPolyX <- reactive({input$trimPolyX})
  minLengthPolyX <- reactive({input$minLengthPolyX})
  
  overrepresentationAnalysis <- reactive({input$overrepresentationAnalysis})
  overrepresentationSampling <- reactive({input$overrepresentationSampling})
  
  correctionOverlap <- reactive({input$correctionOverlap})
  minOverlapLength <- reactive({input$minOverlapLength})
  maxOverlapMismatch <- reactive({input$maxOverlapMismatch})
  maxOverlapMismatchPercentage <- reactive({input$maxOverlapMismatchPercentage})
  
  umi <- reactive({input$umi})
  umiLoc <- reactive({input$umiLoc})
  umiLength <- reactive({input$umiLength})
  umiSkipBaseLength <- reactive({input$umiSkipBaseLength})
  umiPrefix <- reactive({input$umiPrefix})
  umiNoConnection <- reactive({input$umiNoConnection})
  
  output$trim_settings <- DT::renderDataTable({
    tbl <- tibble("setting" = character(), "value" = character())
    
    tbl <- tbl %>% add_row(setting = "Adapter Trimmming", value = as.character(adapterTrimming()))
    if (adapterTrimming()) {
      tbl <- tbl %>% add_row(setting = "Adapter(R1)", value = adapterSequenceRead1())
      if (read_type() == "paired") {
        tbl <- tbl %>% add_row(setting = "Adapter(R2)", value = adapterSequenceRead2())
      }
    }
    
    tbl <- tbl %>% add_row(setting = "Quality Filtering", value = as.character(qualityFiltering()))
    if (qualityFiltering()) {
      tbl <- tbl %>%
        add_row(setting = "Min base quality", value = as.character(qualityFilterPhred())) %>%
        add_row(setting = "Max % unqualified bases", value = as.character(qualityFilterPercent())) %>%
        add_row(setting = "Avg. base quality", value = as.character(averageQualFilter()))
    }
    
    tbl <- tbl %>% add_row(setting = "Max N bases", value = as.character(maxNfilter()))
    
    tbl <- tbl %>% add_row(setting = "Cut low quality tail", value = as.character(cutLowQualTail()))
    if (cutLowQualTail()) {
      tbl <- tbl %>%
        add_row(setting = "Window size: Cut tail", value = as.character(cutTailWindowSize())) %>%
        add_row(setting = "Mean quality: Cut tail", value = as.character(cutTailMeanQual()))
    }
    
    tbl <- tbl %>% add_row(setting = "Cut low quality right", value = as.character(cutSlideWindowRight()))
    if (cutSlideWindowRight()) {
      tbl <- tbl %>%
        add_row(setting = "Window size: Cut right", value = as.character(cutSlideWindowSize())) %>%
        add_row(setting = "Mean quality: Cut right", value = as.character(cutSlideWindowQual()))
    }
    
    tbl <- tbl %>% add_row(setting = "Cut low quality front", value = as.character(cutLowQualFront()))
    if (cutLowQualFront()) {
      tbl <- tbl %>%
        add_row(setting = "Window size: Cut front", value = as.character(cutFrontWindowSize())) %>%
        add_row(setting = "Mean quality: Cut front", value = as.character(cutFrontMeanQual()))
    }
    
    tbl <- tbl %>% add_row(setting = "Length Filtering", value = as.character(lengthFiltering()))
    if (lengthFiltering()) {
      tbl <- tbl %>%
        add_row(setting = "Min length", value = as.character(minReadLength())) %>%
        add_row(setting = "Max length", value = as.character(maxReadLength()))
    }
    
    tbl <- tbl %>%
      add_row(setting = "Trim front(R1)", value = as.character(trimFrontRead1())) %>%
      add_row(setting = "Trim tail(R1)", value = as.character(trimTailRead1())) %>%
      add_row(setting = "Trim to length(R1)", value = as.character(maxLengthRead1()))
    if (read_type() == "paired") {
      tbl <- tbl %>%
        add_row(setting = "Trim front(R2)", value = as.character(trimFrontRead2())) %>%
        add_row(setting = "Trim tail(R2)", value = as.character(trimTailRead2())) %>%
        add_row(setting = "Trim to length(R2)", value = as.character(maxLengthRead2()))
    }
    
    tbl <- tbl %>%
      add_row(setting = "Trim polyG tail", value = as.character(forceTrimPolyG())) %>%
      add_row(setting = "Disable polyG trimming", value = as.character(disableTrimPolyG()))
    if (!disableTrimPolyG()) {
      tbl <- tbl %>% add_row(setting = "Min polyG length", value = as.character(minLengthPolyG()))
    }
    
    tbl <- tbl %>% add_row(setting = "Trim polyX tail", value = as.character(trimPolyX()))
    if (trimPolyX()) {
      tbl <- tbl %>% add_row(setting = "Min polyX length", value = as.character(minLengthPolyX()))
    }
    
    tbl <- tbl %>% add_row(setting = "Overrepresentation Analysis", value = as.character(overrepresentationAnalysis()))
    if (overrepresentationAnalysis()) {
      tbl <- tbl %>% add_row(setting = "Overrepresentation Size", value = as.character(overrepresentationSampling()))
    }
    
    tbl <- tbl %>% add_row(setting = "Base Correction", value = as.character(correctionOverlap()))
    if (correctionOverlap()) {
      tbl <- tbl %>%
        add_row(setting = "Min overlap length", value = as.character(minOverlapLength())) %>%
        add_row(setting = "Max no. of mismatch in overlap", value = as.character(maxOverlapMismatch())) %>%
        add_row(setting = "Max % of mismatch in overlap", value = as.character(maxOverlapMismatchPercentage()))
      
    }
    
    tbl <- tbl %>% add_row(setting = "UMI", value = as.character(umi()))
    if (umi()) {
      tbl <- tbl %>%
        add_row(setting = "Location of UMI", value = umiLoc()) %>%
        add_row(setting = "Length of UMI in read1/read2", value = as.character(umiLength())) %>%
        add_row(setting = "Skip bases after UMI to trim", value = as.character(umiSkipBaseLength())) %>%
        add_row(setting = "UMI Prefix", value = umiPrefix()) %>%
        add_row(setting = "Remove _ between prefix & UMI", value = as.character(umiNoConnection()))
    }
    
    return(tbl)
  }, options = list(dom = 'ltipr', paging = FALSE))
  
  
  
  ############################################### SERVER: TRIMMING ###############################################
  output$select_sample_tr <- renderUI({
    selectInput("fastqSample_tr", "Select FASTQ Sample", choices = allsamples())
  })
  
  observeEvent(input$runtrimming, {
    showModal(modalDialog("Trimming with Rfastp ... ", footer=NULL))
    
    i <- match(input$fastqSample_tr, allsamples())
    trimmed_dir <- paste0(session$token, "/", input$fastqSample_tr)
    
    if (read_type() == 'single') {
      se_dir <- paste0(trimmed_dir, "/", raw_se_dir()[i])
      se_fastq <- list.files(path = se_dir, full.names = T)
      SampleName <- gsub(pattern = fastq_suffix()[i], replacement = "", x = basename(se_fastq))
      
      mapply(rfastp,
             read1 = se_fastq,
             outputFastq = paste0(trimmed_dir, "/", SampleName, "_trimmed"),
             MoreArgs = list(adapterTrimming = adapterTrimming(),
                             adapterSequenceRead1 = adapterSequenceRead1(),
                             trimFrontRead1 = trimFrontRead1(),
                             trimTailRead1 = trimTailRead1(),
                             maxLengthRead1 = maxLengthRead1(),
                             forceTrimPolyG = forceTrimPolyG(),
                             disableTrimPolyG = disableTrimPolyG(),
                             minLengthPolyG = minLengthPolyG(),
                             trimPolyX = trimPolyX(),
                             minLengthPolyX = minLengthPolyX(),
                             cutLowQualTail = cutLowQualTail(),
                             cutSlideWindowRight = cutSlideWindowRight(),
                             cutLowQualFront = cutLowQualFront(),
                             cutFrontWindowSize = cutFrontWindowSize(),
                             cutFrontMeanQual = cutFrontMeanQual(),
                             cutTailWindowSize = cutTailWindowSize(),
                             cutTailMeanQual = cutTailMeanQual(),
                             cutSlideWindowSize = cutSlideWindowSize(),
                             cutSlideWindowQual = cutSlideWindowQual(),
                             qualityFiltering = qualityFiltering(),
                             qualityFilterPhred = qualityFilterPhred(),
                             qualityFilterPercent = qualityFilterPercent(),
                             maxNfilter = maxNfilter(),
                             lengthFiltering = lengthFiltering(),
                             minReadLength = minReadLength(),
                             maxReadLength = maxReadLength(),
                             overrepresentationAnalysis = overrepresentationAnalysis(),
                             overrepresentationSampling = overrepresentationSampling(),
                             umi = umi(),
                             umiLoc = umiLoc(),
                             umiLength = umiLength(),
                             umiPrefix = umiPrefix(),
                             umiSkipBaseLength = umiSkipBaseLength(),
                             umiNoConnection = umiNoConnection()))
      
      FileName <- list.files(path = trimmed_dir, pattern = "_trimmed_R1.fastq.gz$", full.names = F)
      df <- data.frame(FileName, SampleName)
      names(df) <- c("FileName", "SampleName")
      
    } else {
      r1_dir <- paste0(trimmed_dir, "/", raw_per1_dir()[i])
      r2_dir <- paste0(trimmed_dir, "/", raw_per2_dir()[i])
      r1_fastq <- list.files(path = r1_dir, full.names = T)
      r2_fastq <- list.files(path = r2_dir, full.names = T)
      SampleName <- gsub(pattern = fastq_suffix()[i], replacement = "", x = basename(r1_fastq))
      
      mapply(rfastp,
             read1 = r1_fastq,
             read2 = r2_fastq,
             outputFastq = paste0(trimmed_dir, "/", SampleName, "_trimmed"),
             MoreArgs = list(
               adapterTrimming = adapterTrimming(),
               adapterSequenceRead1 = adapterSequenceRead1(),
               adapterSequenceRead2 = adapterSequenceRead2(),
               trimFrontRead1 = trimFrontRead1(),
               trimTailRead1 = trimTailRead1(),
               trimFrontRead2 = trimFrontRead2(),
               trimTailRead2 = trimTailRead2(),
               maxLengthRead1 = maxLengthRead1(),
               maxLengthRead2 = maxLengthRead2(),
               forceTrimPolyG = forceTrimPolyG(),
               disableTrimPolyG = disableTrimPolyG(),
               minLengthPolyG = minLengthPolyG(),
               trimPolyX = trimPolyX(),
               minLengthPolyX = minLengthPolyX(),
               cutLowQualTail = cutLowQualTail(),
               cutSlideWindowRight = cutSlideWindowRight(),
               cutLowQualFront = cutLowQualFront(),
               cutFrontWindowSize = cutFrontWindowSize(),
               cutFrontMeanQual = cutFrontMeanQual(),
               cutTailWindowSize = cutTailWindowSize(),
               cutTailMeanQual = cutTailMeanQual(),
               cutSlideWindowSize = cutSlideWindowSize(),
               cutSlideWindowQual = cutSlideWindowQual(),
               qualityFiltering = qualityFiltering(),
               qualityFilterPhred = qualityFilterPhred(),
               qualityFilterPercent = qualityFilterPercent(),
               maxNfilter = maxNfilter(),
               lengthFiltering = lengthFiltering(),
               minReadLength = minReadLength(),
               maxReadLength = maxReadLength(),
               correctionOverlap = correctionOverlap(),
               minOverlapLength = minOverlapLength(),
               maxOverlapMismatch = maxOverlapMismatch(),
               maxOverlapMismatchPercentage = maxOverlapMismatchPercentage(),
               overrepresentationAnalysis = overrepresentationAnalysis(),
               overrepresentationSampling = overrepresentationSampling(),
               umi = umi(),
               umiLoc = umiLoc(),
               umiLength = umiLength(),
               umiPrefix = umiPrefix(),
               umiSkipBaseLength = umiSkipBaseLength(),
               umiNoConnection = umiNoConnection()))
      
      FileName1 <- list.files(path = trimmed_dir, pattern = "_trimmed_R1.fastq.gz$", full.names = F)
      FileName2 <- list.files(path = trimmed_dir, pattern = "_trimmed_R2.fastq.gz$", full.names = F)
      df <- data.frame(FileName1, FileName2, SampleName)
      names(df) <- c("FileName1", "FileName2", "SampleName")
      
    }
    write.table(df, file = paste0(trimmed_dir, "/", "trimmed.txt"), sep="\t", row.names=FALSE, quote = FALSE)
    
    removeModal()
  })
  
  # Results of trimming
  output$select_trim_fastq <- renderUI({
    json_files <- list.files(path = paste0(session$token, "/", input$fastqSample_tr), pattern = "_trimmed.json$", full.names = F)
    selectInput("trimmed_fastq_file", "Select Trimmed Fastq JSON", choices = json_files)
  })
  
  output$table_trim_summary <- renderTable({
    json <- rjson::fromJSON(file = paste0(session$token, "/", input$fastqSample_tr, "/", input$trimmed_fastq_file))
    qcSummary(json)
  }, rownames = TRUE)
  
  output$trim_base_quality <- renderPlot({
    json <- rjson::fromJSON(file = paste0(session$token, "/", input$fastqSample_tr, "/", input$trimmed_fastq_file))
    curvePlot(json)
  })
  
  output$trim_gc_content <- renderPlot({
    json <- rjson::fromJSON(file = paste0(session$token, "/", input$fastqSample_tr, "/", input$trimmed_fastq_file))
    curvePlot(json, curves = "content_curves")
  })
  
  output$trim_overrepresented <- renderDataTable({
    json <- rjson::fromJSON(file = paste0(session$token, "/", input$fastqSample_tr, "/", input$trimmed_fastq_file))
    
    if (read_type() == "single") {
      df <- data.frame(unlist(json$read1_after_filtering$overrepresented_sequences))
      names(df) <- c("Count")
    } else {
      df1 <- data.frame(unlist(json$read1_after_filtering$overrepresented_sequences))
      names(df1) <- c("Count")
      df1$Read <- "R1"
      
      df2 <- data.frame(unlist(json$read2_after_filtering$overrepresented_sequences))
      names(df2) <- c("Count")
      df2$Read <- "R2"
      
      df <- rbind(df1, df2)
    }
    return(df)
  }, options = list(scrollX = TRUE, scrollY = TRUE))
  
  output$trim_report <- downloadHandler(
    filename = "trim_report.html",
    content = function(file) {
      file.copy(from = paste0(session$token, "/", input$fastqSample_tr, "/", tools::file_path_sans_ext(input$trimmed_fastq_file), ".html"),
                to = file)
    }
  )
  
  ############################################### SERVER: ANNOTATION FILES ###############################################
  clickValues <- reactiveValues(fasta = NULL, txdb = NULL, dxTrack = NULL, trTrack = NULL, geTrack = NULL)
  
  output$refgene_table <- renderDataTable(refgene_info, server = T, filter = 'top', options = list(scrollX = TRUE))
  
  # Downlaod NCBI Files
  output$select_ncbi_fasta <- renderUI({
    validate(need(input$ncbi_accession %in% refgene_info$assembly_accession, message = "Not Valid"))
    
    ftp_path <- paste0(refgene_info[refgene_info$assembly_accession == input$ncbi_accession, "ftp_path"], "/")
    refgene_url <- RCurl::getURL(ftp_path)
    ncbi_fa_options <- rvest::read_html(refgene_url) %>%
      rvest::html_elements(., xpath = ".//a[contains(@href, 'a.gz')]") %>%
      rvest::html_text()
    selectInput("ncbi_fa", "Select NCBI Fasta File", choices = ncbi_fa_options)
  })
  
  output$select_ncbi_gtf <- renderUI({
    validate(need(input$ncbi_accession %in% refgene_info$assembly_accession, message = "Not Valid"))
    
    ftp_path <- paste0(refgene_info[refgene_info$assembly_accession == input$ncbi_accession, "ftp_path"], "/")
    refgene_url <- RCurl::getURL(ftp_path)
    ncbi_gtf_options <- rvest::read_html(refgene_url) %>%
      rvest::html_elements(., xpath = ".//a[contains(@href, 'f.gz')]") %>%
      rvest::html_text()
    selectInput("ncbi_gtf", "Select NCBI Annotation File", choices = ncbi_gtf_options)
  })
  
  # Check whether gene ids of fasta & gtf match
  observeEvent(input$prepare_ann, {
    showModal(modalDialog("Prepare for Aligning ... ", footer=NULL))
    
    if (input$upload_how) {
      full_fasta <- Biostrings::readDNAStringSet(input$upload_fasta$datapath)
      
      txdb <- GenomicFeatures::makeTxDbFromGFF(input$upload_gtf$datapath, format = input$gtf_format)
      txdb_names <- GenomeInfoDb::seqlevels(txdb)
      
    } else {
      options(timeout = 600)
      ftp_path <- paste0(refgene_info[refgene_info$assembly_accession == input$ncbi_accession, "ftp_path"], "/")
      
      fasta_url <- paste0(ftp_path, input$ncbi_fa)
      threadr::download_ftp_file(fasta_url, input$ncbi_fa)
      
      gtf_url <- paste0(ftp_path, input$ncbi_gtf)
      threadr::download_ftp_file(gtf_url, input$ncbi_gtf)
      
      full_fasta <- Biostrings::readDNAStringSet(input$ncbi_fa)
      
      txdb <- GenomicFeatures::makeTxDbFromGFF(input$ncbi_gtf, format = input$gtf_format)
      txdb_names <- GenomeInfoDb::seqlevels(txdb)
      
    }
    
    names(full_fasta) <- sapply(strsplit(names(full_fasta)," "), `[`, 1)
    fasta_names <- names(full_fasta)
    fasta_names <- fasta_names[order(fasta_names, decreasing = F)]
    txdb_names <- txdb_names[order(txdb_names, decreasing = F)]
    
    if (input$upload_how & identical(fasta_names, txdb_names)) {
      clickValues$fasta <- input$upload_fasta$datapath
      clickValues$txdb <- txdb
      
    } else if (!input$upload_how & identical(fasta_names, txdb_names)) {
      Biostrings::writeXStringSet(full_fasta, filepath = paste0(input$ncbi_accession, ".fasta"))
      clickValues$fasta <- paste0(input$ncbi_accession, ".fasta")
      clickValues$txdb <- input$ncbi_gtf
      
    } else if (length(fasta_names[!fasta_names %in% txdb_names]) != 0) {
      if (input$upload_how) {
        filename <- tools::file_path_sans_ext(input$upload_fasta$name)
      } else {
        filename <- input$ncbi_accession
      }
      
      Biostrings::writeXStringSet(full_fasta[txdb_names], filepath = paste0(filename, ".fasta"))
      clickValues$fasta <- paste0(filename, ".fasta")
      clickValues$txdb <- txdb
      
    } else {
      GenomeInfoDb::seqlevels(txdb) <- fasta_names
      
      if (input$upload_how) {
        clickValues$fasta <-input$upload_fasta$datapath
      } else {
        Biostrings::writeXStringSet(full_fasta, filepath = paste0(input$ncbi_accession, ".fasta"))
        clickValues$fasta <- paste0(input$ncbi_accession, ".fasta")
      }
      clickValues$txdb <- txdb
    }
    
    output$select_chrom <- renderUI({
      chrom <- seqlevels(clickValues$txdb)
      pickerInput("plot_chrom", "Select Chromosome to View Coverage", choices = chrom)
    })
    
    output$select_chrom_range <- renderUI({
      genes <- genes(clickValues$txdb,filter=list(tx_chrom = input$plot_chrom))
      df <- data.frame(genes)
      
      min_range <- min(df$start)
      max_range <- max(df$end)
      
      q1 <- floor((max_range - min_range)/4 + min_range)
      q2 <- floor((max_range - min_range)/2 + min_range)
      sliderInput("plot_chrom_range", "Range to plot", min = min_range, max = max_range, value = c(q1, q2))
    })
    
    removeModal()
    updateActionButton(session, "prepare_ann", label = "Ann. Files Ready!", icon = icon("check"))
  })
  
  
  
  ############################################### SERVER: ALIGNING ###############################################
  output$select_sample_al <- renderUI({
    selectInput("fastqSample_al", "Select FASTQ Sample", choices = allsamples())
  })
  
  observeEvent(input$runalign, {
    showModal(modalDialog("Aligning ... ", footer=NULL))
    sample_name <- input$fastqSample_al
    fastq_dir <- paste0(session$token, "/", sample_name)
    sampleFile <- switch (input$trimmed_untrimmed,
                          "Trimmed" = paste0(fastq_dir, "/", "trimmed.txt"),
                          "Untrimmed" = paste0(fastq_dir, "/", "untrimmed.txt"))
    
    # Run aligning
    cl <- parallel::makeCluster(parallel::detectCores()/2)
    
    proj <- qAlign(sampleFile = sampleFile,
                   genome = clickValues$fasta,
                   aligner = input$aligner,
                   maxHits = input$maxHits,
                   paired = switch(read_type(), 'single' = 'no', 'paired' = input$paired),
                   splicedAlignment = input$splicedAlignment,
                   alignmentsDir = fastq_dir,
                   clObj = cl)
    
    # Generate counts
    showModal(modalDialog("Generating Counts ... ", footer=NULL))
    # GENE
    gc <- qCount(
      proj = proj,
      query = clickValues$txdb,
      reportLevel = "gene",
      selectReadPosition = input$selectReadPosition,
      orientation = input$orientation,
      useRead = input$useRead,
      includeSpliced = input$includeSpliced,
      includeSecondary = input$includeSecondary,
      mapqMin = input$mapqlty[1],
      mapqMax = input$mapqlty[2],
      clObj = cl
    )
    
    output$gene_counts <- downloadHandler(
      filename = paste0(sample_name, "_gene.csv"),
      content = function(file) {
        write.csv(gc, file, row.names=T, quote = F)
      }
    )
    
    output$gene_preview <- renderDataTable({
      datatable(head(gc), options = list(dom = 't'))
    })
    
    
    # EXON
    ec <- qCount(
      proj = proj,
      query = clickValues$txdb,
      reportLevel = "exon",
      selectReadPosition = input$selectReadPosition,
      orientation = input$orientation,
      useRead = input$useRead,
      includeSpliced = input$includeSpliced,
      includeSecondary = input$includeSecondary,
      mapqMin = input$mapqlty[1],
      mapqMax = input$mapqlty[2],
      clObj = cl
    )
    
    output$exon_counts <- downloadHandler(
      filename = paste0(sample_name, "_exon.csv"),
      content = function(file) {
        write.csv(ec, file, row.names=T, quote = F)
      }
    )
    
    output$exon_preview <- renderDataTable({
      datatable(head(ec), options = list(dom = 't'))
    })
    
    
    # PROMOTER
    pc <- qCount(
      proj = proj,
      query = clickValues$txdb,
      reportLevel = "promoter",
      selectReadPosition = input$selectReadPosition,
      orientation = input$orientation,
      useRead = input$useRead,
      includeSpliced = input$includeSpliced,
      includeSecondary = input$includeSecondary,
      mapqMin = input$mapqlty[1],
      mapqMax = input$mapqlty[2],
      clObj = cl
    )
    
    output$prm_counts <- downloadHandler(
      filename = paste0(sample_name, "_promoter.csv"),
      content = function(file) {
        write.csv(pc, file, row.names=T, quote = F)
      }
    )
    
    output$prm_preview <- renderDataTable({
      datatable(head(pc), options = list(dom = 't'))
    })
    
    # JUNCTION
    if (input$splicedAlignment) {
      jc <- qCount(
        proj = proj,
        query = NULL,
        reportLevel = "junction",
        selectReadPosition = input$selectReadPosition,
        orientation = input$orientation,
        useRead = input$useRead,
        includeSpliced = input$includeSpliced,
        includeSecondary = input$includeSecondary,
        mapqMin = input$mapqlty[1],
        mapqMax = input$mapqlty[2],
        clObj = cl
      )
      
      output$jun_counts <- downloadHandler(
        filename = paste0(sample_name, "_junction.csv"),
        content = function(file) {
          write.csv(jc, file, row.names=T, quote = F)
        }
      )
      
      output$jun_preview <- renderDataTable({
        datatable(head(data.frame(jc)), options = list(dom = 't'))
      })
    }
    
    # Generate reports
    showModal(modalDialog("Quality Reporting ... ", footer=NULL))
    qQCReport(proj, pdfFilename = paste0(session$token, "/", sample_name, "/", sample_name, ".pdf"))
    
    # Export wig files
    fastq_names <- read.table(sampleFile, header = T)$SampleName
    qExportWig(proj, file = paste0(fastq_dir, "/", fastq_names, ".wig.gz"))
    
    parallel::stopCluster(cl)
    removeModal()
  })
  
  ############################################### SERVER: ALIGNMENT STATISTICS ###############################################
  output$select_sample_as <- renderUI({
    selectInput("fastqSample_as", "Select FASTQ Sample", choices = allsamples())
  })
  
  output$alignment_report <- downloadHandler(
    filename = "alignment_report.pdf",
    content = function(file) {
      file.copy(from = paste0(session$token, "/", input$fastqSample_as, "/", input$fastqSample_as, ".pdf"),
                to = file)
    }
  )
  
  observeEvent(input$view_statistics, {
    showModal(modalDialog("Generating Alignment Statistics ... ", footer=NULL))
    # get bam files
    bam_files <- list.files(paste0(session$token, "/", input$fastqSample_as), pattern = ".bam$", full.names = T)
    
    # generate rds files
    params <- ramwas::ramwasParameters(dirbam = getwd(),
                                       bamnames = bam_files,
                                       dirrqc = paste0(session$token, "/", input$fastqSample_as),
                                       minscore = input$minscore,
                                       minfragmentsize = input$minfragmentsize,
                                       maxfragmentsize = input$maxfragmentsize,
                                       maxrepeats = input$maxrepeats,
                                       cputhreads = 4)
    ramwas1scanBams(params)
    
    # compile qc statistics
    bamname <- tools::file_path_sans_ext(basename(bam_files))
    rbamlist = vector("list")
    for (i in bamname) {
      rdsqcfile = paste0(session$token, "/", input$fastqSample_as, "/", i, ".qc.rds")
      rbamlist[[i]] = readRDS(rdsqcfile)
    }
    qclist <- lapply(rbamlist, `[[`, "qc")
    
    ## The number of BAM files
    nbams <- Reduce('+', lapply(qclist, `[[`, 'nbams'))
    
    ## Total number of reads in the BAM file(s)
    reads.total <- Reduce('+', lapply(qclist, `[[`, 'reads.total'))
    
    ## Number of reads aligned to the reference genome
    reads.aligned <- Reduce('+', lapply(qclist, `[[`, 'reads.aligned'))
    
    ## Number of reads that passed minimum score filter and are recorded
    reads.recorded <- Reduce('+', lapply(qclist, `[[`, 'reads.recorded'))
    
    ## Number of reads after removal of duplicate reads
    reads.recorded.no.repeats <- Reduce('+', lapply(qclist, `[[`, 'reads.recorded.no.repeats'))
    
    ## Number of recorded reads aligned to the forward and reverse strands respectively (with duplicates)
    frwrev <- Reduce('+', lapply(qclist, `[[`, 'frwrev'))
    
    ## Number of recorded reads aligned to the forward and reverse strands respectively (without duplicates)
    frwrev.no.repeats <- Reduce('+', lapply(qclist, `[[`, 'frwrev.no.repeats'))
    
    ## Distribution of the read scores
    hist.score1 <- floor(Reduce('+', lapply(qclist, `[[`, 'hist.score1')) / nbams)
    bf.hist.score1 <- floor(Reduce('+', lapply(qclist, `[[`, 'bf.hist.score1')) / nbams)
    
    ## Distribution of the length of the aligned part of the reads
    hist.length.matched <- floor(Reduce('+', lapply(qclist, `[[`, 'hist.length.matched')) / nbams)
    bf.hist.length.matched <- floor(Reduce('+', lapply(qclist, `[[`, 'bf.hist.length.matched')) / nbams)
    
    ## Distribution of edit distance between the aligned part of the read and the reference genome
    hist.edit.dist1 <- floor(Reduce('+', lapply(qclist, `[[`, 'hist.edit.dist1')) / nbams)
    bf.hist.edit.dist1 <- floor(Reduce('+', lapply(qclist, `[[`, 'bf.hist.edit.dist1')) / nbams)
    
    output$fig1 <- renderPlotly({
      plot_ly(data = data.frame(group = c("Aligned", "Not Aligned"), value = c(reads.aligned, reads.total - reads.aligned)),
              labels = ~group,
              values = ~value,
              insidetextorientation='horizontal',
              textposition = 'inside',
              textinfo = 'label+percent',
              hoverinfo = 'text',
              text = ~paste0(group, '\n', 'Total: ', value, "\n", "Average: ", floor(value / nbams)),
              marker =list(colors=c("midnightblue", "powderblue"))) %>%
        add_pie(hole = 0.6) %>%
        layout(showlegend = F,
               margin = list(l = 20, r = 20),
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = F),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = F),
               annotations = list(text = paste0(round(reads.aligned/reads.total * 100), "%", "\nAligned"),
                                  "showarrow"=F,
                                  font=list(size = 15, color = "grey"))
        )
    })
    
    output$fig2 <- renderPlotly({
      plot_ly(data = data.frame(group = c("Pass", "Failed"), value = c(reads.recorded, reads.aligned - reads.recorded)),
              labels = ~group,
              values = ~value,
              insidetextorientation='horizontal',
              textposition = 'inside',
              textinfo = 'label+percent',
              hoverinfo = 'text',
              text = ~paste0(group, '\n', 'Total: ', value, "\n", "Average: ", floor(value / nbams)),
              marker =list(colors=c("midnightblue", "powderblue"))) %>%
        add_pie(hole = 0.6) %>%
        layout(showlegend = F,
               margin = list(l = 20, r = 20),
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = F),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = F),
               annotations = list(text = paste0(round(reads.recorded/reads.aligned * 100), "%", "\nPassed"),
                                  "showarrow"=F,
                                  font=list(size = 15, color = "grey"))
        )
    })
    
    output$fig3 <- renderPlotly({
      plot_ly(data = data.frame(group = c("w/o Rpts", "w Rpts"), value = c(reads.recorded.no.repeats, reads.recorded - reads.recorded.no.repeats)),
              labels = ~group,
              values = ~value,
              insidetextorientation='horizontal',
              textposition = 'inside',
              textinfo = 'label+percent',
              hoverinfo = 'text',
              text = ~paste0(group, '\n', 'Total: ', value, "\n", "Average: ", floor(value / nbams)),
              marker =list(colors=c("midnightblue", "powderblue"))) %>%
        add_pie(hole = 0.6) %>%
        layout(showlegend = F,
               margin = list(l = 20, r = 20),
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = F),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = F),
               annotations = list(text = paste0(round(reads.recorded.no.repeats/reads.recorded * 100), "%", "\nNot Repeated"),
                                  "showarrow"=F,
                                  font=list(size = 15, color = "grey"))
        )
    })
    
    output$fig4 <- renderPlotly({
      plot_ly(data = data.frame(group = c("Forward", "Reverse"), value = c(frwrev[1], frwrev[2])),
              labels = ~group,
              values = ~value,
              insidetextorientation='horizontal',
              textposition = 'inside',
              textinfo = 'label+percent',
              hoverinfo = 'text',
              text = ~paste0(group, '\n', 'Total: ', value, "\n", "Average: ", floor(value / nbams)),
              marker =list(colors=c("midnightblue", "powderblue"))) %>%
        add_pie(hole = 0.6) %>%
        layout(showlegend = F,
               margin = list(l = 20, r = 20),
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = F),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = F),
               annotations = list(text = paste0(round(frwrev[1]/(frwrev[1] + frwrev[2]) * 100), "%", "\nForward"),
                                  "showarrow"=F,
                                  font=list(size = 15, color = "grey"))
        )
    })
    
    output$fig5 <- renderPlotly({
      plot_ly(data = data.frame(group = c("Forward", "Reverse"), value = c(frwrev.no.repeats[1], frwrev.no.repeats[2])),
              labels = ~group,
              values = ~value,
              insidetextorientation='horizontal',
              textposition = 'inside',
              textinfo = 'label+percent',
              hoverinfo = 'text',
              text = ~paste0(group, '\n', 'Total: ', value, "\n", "Average: ", floor(value / nbams)),
              marker =list(colors=c("midnightblue", "powderblue"))) %>%
        add_pie(hole = 0.6) %>%
        layout(showlegend = F,
               margin = list(l = 20, r = 20),
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = F),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = F),
               annotations = list(text = paste0(round(frwrev.no.repeats[1]/(frwrev.no.repeats[1] + frwrev.no.repeats[2]) * 100), "%", "\nForward"),
                                  "showarrow"=F,
                                  font=list(size = 15, color = "grey"))
        )
    })
    
    output$fig6 <- renderPlot({
      par(mfrow=c(1,2))
      plot(hist.score1)
      plot(bf.hist.score1)
    })
    
    output$fig6_summary <- renderText({
      return(paste0("Average alignment score, after filter: ", qcmean(hist.score1), "\n",
                    "Average alignment score, no filter: ", qcmean(bf.hist.score1))
      )
    })
    
    output$fig7 <- renderPlot({
      par(mfrow = c(1,2))
      plot(hist.length.matched)
      plot(bf.hist.length.matched)
    })
    
    output$fig7_summary <- renderText({
      return(paste0("Average aligned length: ", round(qcmean(hist.length.matched)), "\n",
                    "Average aligned length, no filter: ", round(qcmean(bf.hist.length.matched)))
      )
    })
    
    output$fig8 <- renderPlot({
      par(mfrow = c(1,2))
      plot(hist.edit.dist1)
      plot(bf.hist.edit.dist1)
    })
    
    output$fig8_summary <- renderText({
      return(paste0("Average edit distance, after filter: ", qcmean(hist.edit.dist1), "\n",
                    "Average edit distance, no filter: ", qcmean(bf.hist.edit.dist1))
      )
    })
    
    removeModal()
  })
  
  
  
  ############################################### SERVER: GENOME BROWSER ###############################################
  output$select_sample_gb <- renderUI({
    selectInput("fastqSample_gb", "Select FASTQ Sample", choices = allsamples())
  })
  
  
  output$selected_ranges <- renderDT({
    gr <- GRanges(seqnames = input$plot_chrom,
                  ranges = IRanges(start = input$plot_chrom_range[1], input$plot_chrom_range[2]))
    
    x <- switch (input$sel_ranges_type,
                 "gene" = genes(clickValues$txdb),
                 "transcript" = transcripts(clickValues$txdb),
                 "exon" = exons(clickValues$txdb)
    )
    
    return(data.frame(subsetByOverlaps(x, gr)))
  }, filter = "top", selection = "none", options = list(pageLength = 15))
  
  
  observeEvent(input$view_cov, {
    wig_file <- list.files(paste0(session$token, "/", input$fastqSample_gb), pattern = ".wig.gz$", full.names = T)
    hist_col <- colorspace::qualitative_hcl(n = length(wig_file), palette = "Dark 3")
    
    axTrack <- GenomeAxisTrack()
    dxTrack <- c()
    for (i in 1:length(wig_file)) {
      x <- DataTrack(range = wig_file[i],
                     name = gsub(pattern = ".wig.gz", replacement = "", x = basename(wig_file[i])),
                     type = "h",
                     col = hist_col[i])
      dxTrack <- c(dxTrack, x)
    }
    
    ge <- genes(clickValues$txdb)
    geTrack <- GeneRegionTrack(ge, showId = T,  name = "Gene",
                               background.title = "darkblue", background.panel = "#FFFEDB")
    symbol(geTrack) <- ge$gene_id
    
    tr <- transcripts(clickValues$txdb)
    if (length(tr[is.na(tr$tx_name)]) == 0) {
      trTrack <- GeneRegionTrack(clickValues$txdb, showId = T,  name = "Transcript",
                                 background.title = "red", background.panel = "#FFFEDB")
      
      output$cov_plot <- renderPlot({
        plotTracks(as.list(c(axTrack, dxTrack, trTrack, geTrack)),
                   chromosome = input$plot_chrom,
                   from = input$plot_chrom_range[1],
                   to = input$plot_chrom_range[2],
                   extend.right = input$extend_right,
                   extend.left = input$extend_left)
      })
      
    } else {
      output$cov_plot <- renderPlot({
        plotTracks(as.list(c(axTrack, dxTrack,geTrack)),
                   chromosome = input$plot_chrom,
                   from = input$plot_chrom_range[1],
                   to = input$plot_chrom_range[2],
                   extend.right = input$extend_right,
                   extend.left = input$extend_left)
      })
    }
  })
  
}


# Run the application 
shinyApp(ui = ui, server = server)
