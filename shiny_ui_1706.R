ui <- tagList(useShinydashboard(),
  navbarPage(id = "navbar", title = "",
    navbarMenu("Quality Control",
      tabPanel("FastQC",
        sidebarLayout(
          sidebarPanel(
            radioButtons("read_type", "Choose Read Type", choices = c("single", "paired"), inline = T),
            actionButton('addSample', '', icon = icon('plus'))
          ),
          
          mainPanel(
            tabsetPanel(
              tabPanel("Multi-QC",
                       fluidRow(
                         div(style = "margin-top: +20px"),
                         column(4, uiOutput("select_sample_qc")),
                         column(4, selectInput("multiqc_plot_type", "Type of Multi-QC Plot", choices = multiqc.plot.type)),
                         column(4,
                                div(style = "margin-top: +25px"),
                                actionButton("runfastqc", "Run FASTQC")
                         )
                       ),
                       
                       plotOutput("multiqc_plot")
              ),
              
              tabPanel("Single-QC",
                       fluidRow(
                         div(style = "margin-top: +20px"),
                         column(4, uiOutput("select_qc_file")),
                         column(4, selectInput("singleqc_plot_type", "Type of Single-QC Plot", choices = singleqc.plot.type))
                       ),
                       
                       plotOutput("singleqc_plot")
              ),
              
              tabPanel("More Info",
                verbatimTextOutput("test")
              )
            )
          )
        )
      ),
      
      tabPanel("Trimming",
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
                      radioButtons("adapterTrimming", "Trim adaptor", choices = c("no", "yes"), inline = TRUE),
  
                      conditionalPanel(condition = "input.adapterTrimming == 'yes'",
  
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
  
                      radioButtons("qualityFiltering", "Filter by quality", choices = c("no", "yes"), inline = TRUE),
  
                      conditionalPanel(condition = "input.qualityFiltering == 'yes'",
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
                               conditionalPanel(condition = "input.cutLowQualTail == 'yes'",
                                                sliderInput("cutTailWindowSize", "Window size: Cut tail", min = 1, max=1000, value = 4),
                                                numericInput("cutTailMeanQual", "Mean quality: Cut tail", value = 20)),
  
                               conditionalPanel(condition = "input.cutSlideWindowRight == 'yes'",
                                                sliderInput("cutSlideWindowSize", "Window size: Cut right", min = 1, max=1000, value = 4),
                                                numericInput("cutSlideWindowQual", "Mean quality: Cut right", value = 20)),
  
                               conditionalPanel(condition = "input.cutLowQualFront == 'yes'",
                                                sliderInput("cutFrontWindowSize", "Window size: Cut front", min = 1, max=1000, value = 4),
                                                numericInput("cutFrontMeanQual", "Mean quality: Cut front", value = 20))
                      ),
  
                      tabPanel("Enable", value = "Enable",
                        radioButtons("cutLowQualTail", "Cut tail", choices = c("no", "yes"), inline = TRUE),
                        div(style = "margin-top: -15px"),
                        helpText(HTML("<font size=-1>
                                        Move a sliding window from tail (3') to the front, drop the bases in the window if its
                                        mean quality < threshold, stop otherwise. This is similar to Trimmomatic TRAILING method.
                                      </font>")),
  
                        div(style = "margin-top: +20px"),
  
                        radioButtons("cutSlideWindowRight", "Cut right", choices = c("no", "yes"), inline = TRUE),
                        div(style = "margin-top: -15px"),
                        helpText(HTML("<font size=-1>
                                        Move a sliding window from front to tail, if meet one window with mean quality < threshold,
                                        drop the bases in the window and the right part, and then stop. This is simliar to Trimmomatic
                                        SLIDINGWINDOW method.
                                      </font>")),
  
                        div(style = "margin-top: +20px"),
  
                        radioButtons("cutLowQualFront", "Cut front", choices = c("no", "yes"), inline = TRUE),
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
  
                        radioButtons("lengthFiltering", "Filter by length", choices = c("no", "yes"), inline = TRUE),
                        div(style = "margin-top: -15px"),
                        helpText(HTML("<font size=-1>Note that length filtering is applied as a last step.</font>")),
                        
                        div(style = "margin-top: +20px"),
  
                        conditionalPanel(condition = "input.lengthFiltering == 'yes'",
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
  
                      radioButtons("forceTrimPolyG", "Force trim polyG tail", choices = c("no", "yes"), inline = TRUE),
                      div(style = "margin-top: -15px"),
                      helpText(HTML("<font size=-1>
                                      <li>For Illumina NextSeq/NovaSeq data, polyG can happen in read tails since G means no signal in the Illumina two-color systems.</li>
                                      <li>NextSeq/NovaSeq data is detected by the machine ID in the FASTQ records and this feature is automatically enabled for such data.</li>
                                    </font>")),
  
                      div(style = "margin-top: +20px"),
  
                      radioButtons("disableTrimPolyG", "Disable polyG tail trimming", choices = c("no", "yes"), inline = TRUE),
                      div(style = "margin-top: -15px"),
                      helpText(HTML("<font size=-1>Disable automatic polyG trimming for Illumina NextSeq/NovaSeq data.</font>")),
  
                      div(style = "margin-top: +20px"),
  
                      numericInput("minLengthPolyG", "Min polyG tail length", value = 10),
                      div(style = "margin-top: -15px"),
                      helpText(HTML("<font size=-1>Min length to detect polyG</font>")),
  
                      div(style = "margin-top: +20px"),
  
                      radioButtons("trimPolyX",  "Force trim polyX tail", choices = c("no", "yes"), inline = TRUE),
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
  
                      radioButtons("overrepresentationAnalysis", "Enable overrepresentation analysis", choices = c("no", "yes"), inline = TRUE),
  
                      conditionalPanel(condition = "input.overrepresentationAnalysis == 'yes'",
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
                                       radioButtons("correctionOverlap", "Run correction overlap", choices = c("no", "yes"), inline = TRUE),
                                       div(style = "margin-top: -15px"),
                                       helpText(HTML('<font size=-1>
                                                        If a proper overlap between a pair of reads is found, it can correct mismatched base pairs or if 1 base is with
                                                        high quality while the other is with ultra low quality. If a base is corrected, the quality of its paired base
                                                        will be assigned to it so that they will share the same quality.
                                                     </font>')),
                                       
                                       div(style = "margin-top: +20px"),
  
                                       conditionalPanel(condition = "input.correctionOverlap == 'yes'",
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
                      radioButtons("umi", "Preprocess UMI", choices = c("no", "yes"), inline = TRUE),
  
                      conditionalPanel(condition = "input.umi == 'yes'",
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
                                       
                                       radioButtons("umiNoConnection", "Remove '_' between UMI prefix & UMI string", choices = c("no", "yes"), inline = T))
                    )
                  )
                )
              ),
  
              tabPanel("Results",
                h3(strong("Trim Result")),
                uiOutput("select_trim_file"),
                
                fluidRow(
                  column(6, 
                    h4(span("Summary:", style = "color:#73C6B6;font-weight: bold;")),
                    tableOutput("table_trim_summary")
                  ),
                  
                  column(6,
                    h4(span("Over-Represented Sequences:", style = "color:#73C6B6;font-weight: bold;")),
                    div(dataTableOutput("trim_overrepresented"), style = "font-size:80%; height:330px; overflow-y: scroll; overflow-x: scroll")
                  )
                ),
                
                h4(span("Base Quality Comparison:", style = "color:#73C6B6;font-weight: bold;")),
                withLoader(plotOutput("trim_base_quality"), type = "html", loader = "dnaspin"),
  
                div(style = "margin-top: +20px"),
  
                h4(span("Base Content Comparison:", style = "color:#73C6B6;font-weight: bold;")),
                withLoader(plotOutput("trim_gc_content"), type = "html", loader = "dnaspin")
              )
            )
          )
        )
      )
    ),
    
    navbarMenu("Aligning",
      tabPanel("Generate Counts",
        sidebarLayout(
          sidebarPanel(
            tabsetPanel(
              tabPanel("Upload",
                div(style = "margin-top: +20px"),
                radioButtons("kingdom", "Select Genome Kingdom", choices = kingdoms, inline = T),
                uiOutput("select_organism"),
                uiOutput("select_accession"),
               
                radioButtons("fasta_how", "Upload FASTA File?", choices = c("no", "yes"), inline = T),
                conditionalPanel(condition = "input.fasta_how == 'no'", 
                                 actionButton("get_fasta", "Get FASTA"),
                                 div(style = "margin-top: +20px"),
                                 uiOutput("select_fasta_chrom"),
                                 actionButton("load_fasta", "Load FASTA"),
                                 div(style = "margin-top: +20px")),
                conditionalPanel(condition = "input.fasta_how == 'yes'",
                                 fileInput("fasta_file", "Choose FASTA File")),
                
                radioButtons("gtf_how", "Upload GTF File?", choices = c("no", "yes"), inline = T),
                conditionalPanel(condition = "input.gtf_how == 'no'", 
                                 actionButton("get_txdb", "Get TXDB"),
                                 div(style = "margin-top: +20px"),
                                 uiOutput("select_txdb_chrom"),
                                 actionButton("load_txdb", "Load TXDB")),
                conditionalPanel(condition = "input.gtf_how == 'yes'",
                                 fileInput("gtf_file", "Choose GTF File"))
              ),
              
              tabPanel("Settings",
                div(style = "margin-top: +20px"),
                uiOutput("select_sample_al"),
                radioButtons("trim_untrimmed", "Select Trimmed for Aligning?", choices = c("no", "yes"), inline = T),
                
                radioButtons("paired", "Type of (PE) Library", choices = c("no", "fr", "ff", "rf"), inline = TRUE),
                div(style = "margin-top: -15px"),
                helpText(HTML("<font size=-1>
                                no: Single reads, fr: Forward/Reverse, ff: Foward/Foward, rf: Reverse/Foward
                              </font>")),
                
                div(style = "margin-top: +20px"),
                
                prettySwitch("splicedAlignment", strong("Spiced Alignment?"), bigger = TRUE, fill = TRUE, status = "success"),
                radioButtons("aligner", "Select Aligner", choices = c("Rbowtie", "Rhisat2"), inline = TRUE),
                div(style = "margin-top: -15px"),
                helpText(HTML("<font size=-1>Hisat2 is the recommended setting for spliced alignment.</font>")),
                
                div(style = "margin-top: +20px"),
                
                numericInput("maxHits", "Max Mapping Positions Per Read", value = 1),
                div(style = "margin-top: -15px"),
                helpText(HTML("<font size=-1>
                                Set the maximal no. of allowed positions per read.</br>
                                If a read produces more that the specified no. of alignments, no alignments will be reported.</br>
                                In case of a multi-mapping read, a single alignment is randomly selected.
                              </font>")),
                
                div(style = "margin-top: +20px"),
                
                radioButtons("selectReadPosition", "Read Position", choices = c("start", "end"), inline = TRUE),
                div(style = "margin-top: -15px"),
                helpText(HTML("<font size=-1>
                                The part of the alignment that has to be contained within a query region to produce an overlap.<br>
                                Start: Start of the alignment<br>
                                End: End of the alignment
                              </font>")),
                
                div(style = "margin-top: +20px"),
                
                radioButtons("orientation", "Orientation of Alignments", choices = c("any", "same", "opposite"), inline = TRUE),
                div(style = "margin-top: -15px"),
                helpText(HTML("<font size=-1>
                                The required orientation of the alignments relative to the query region in order to be counted.<br>
                                any: Count alignment on same & opp strand<br>
                                same: Count alignment only on same strand<br>
                                opposite: Count alignment only on opp strand
                              </font>")),
                
                div(style = "margin-top: +20px"),
                
                radioButtons("useRead", "Read Mate to Count Alignments", choices = c("any", "first", "last"), inline = TRUE),
                div(style = "margin-top: -15px"),
                helpText(HTML("<font size=-1>
                                The read mate whose alignments should be counted.<br>
                                any: Count all alignments<br>
                                first: Count alignment only on the 1st strand<br>
                                last: Count alignment only on the last strand
                              </font>")),
                
                div(style = "margin-top: +20px"),
                
                prettySwitch("includeSpliced", strong("Include Spliced Alignments When Counting?"), bigger = TRUE, fill = TRUE, status = "success", value = TRUE),
                prettySwitch("includeSecondary", strong("Include Flagged Alignments When Counting?"), bigger = TRUE, fill = TRUE, status = "success", value = TRUE),
                
                div(style = "margin-top: +30px"),
                
                sliderInput("mapqlty", "Mapping Quality of Alignments When Counting", min = 0, max = 255, value = c(0, 255), step = 1),
                
                actionButton("runalign", "Align")
              )
            )
          ),
          
          mainPanel(
            dataTableOutput("refseq_organisms")
          )
        )
      ),
      
      tabPanel("Genome Browser")
    )
  )
)
