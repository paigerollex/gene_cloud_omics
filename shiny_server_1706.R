server <- function(input, output, session) {
  dir.create(session$token)
  clickValues <- reactiveValues(full_fasta = NULL, full_txdb = NULL)
  
  observeEvent(input$addSample, {
    id <- sprintf('fastqSample%s', input$addSample)
    insertUI(
      selector = "#addSample",
      where = "beforeBegin",
      ui = fastqUI(id)
    )
    
    fastqServer(id, input$read_type, session$token)
  })
  
  # Get all file inputs
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

  output$select_sample_qc <- renderUI({
    selectInput("fastqSample_qc", "Select FASTQ Sample", choices = allsamples())
  })
  
  output$select_sample_tr <- renderUI({
    selectInput("fastqSample_tr", "Select FASTQ Sample", choices = allsamples())
  })
  
  output$select_sample_al <- renderUI({
    selectInput("fastqSample_al", "Select FASTQ Sample", choices = allsamples())
  })
  
  
  ################################################################ FASTQC ################################################################ 
  observeEvent(input$runfastqc, {
    index <- match(input$fastqSample_qc, allsamples())
    sample_dir <- paste0(session$token, "/", input$fastqSample_qc)
    
    if (input$read_type == "single") {
      fastqc(fq.dir = paste0(sample_dir, "/", raw_se_dir()[index]),
             qc.dir = paste0(sample_dir, "/", "fastqc"),
             fastqc.path = paste0(getwd(), "/", "FastQC/fastqc"))
    } else {
      fastqc(fq.dir = paste0(sample_dir, "/", raw_per1_dir()[index]),
             qc.dir = paste0(sample_dir, "/", "fastqc"),
             fastqc.path = paste0(getwd(), "/", "FastQC/fastqc"))
      
      fastqc(fq.dir = paste0(sample_dir, "/", raw_per2_dir()[index]),
             qc.dir = paste0(sample_dir, "/", "fastqc"),
             fastqc.path = paste0(getwd(), "/", "FastQC/fastqc"))
    }
  })
  
  fastqc_dir <- reactive({paste0(session$token, "/", input$fastqSample_qc, "/", "fastqc")})
  
  fastqc_files <- reactive({
    req(input$runfastqc)
    list.files(fastqc_dir(), pattern = "fastqc.zip$", full.names = TRUE)
  })
  
  fdl <- reactive({FastqcDataList(fastqc_files())})
  
  output$multiqc_plot <- renderPlot({
    req(input$multiqc_plot_type)
    fastqc_plot(fdl(), input$multiqc_plot_type)
  })
  
  output$select_qc_file <- renderUI({
    selectInput("raw_fastq_file", "Select FASTQC File", choices = basename(fastqc_files()))
  })
  
  output$singleqc_plot <- renderPlot({
    req(input$raw_fastq_file, input$singleqc_plot_type)
    index <- match(input$raw_fastq_file, basename(fastqc_files()))
    which_fdl <- fdl()[[fastqc_files()[index]]]
    fastqc_plot(which_fdl, input$singleqc_plot_type)
  })
  
  ######################################################### Trim Settings ######################################################### 
  read_type <- reactive({input$read_type})
  
  adapterTrimming <- reactive({switch(input$adapterTrimming, "no" = FALSE, "yes" = TRUE)})
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
  
  qualityFiltering <- reactive({switch(input$qualityFiltering, "no" = FALSE, "yes" = TRUE)})
  qualityFilterPhred <- reactive({input$qualityFilterPhred})
  qualityFilterPercent <- reactive({input$qualityFilterPercent})
  averageQualFilter <- reactive({input$averageQualFilter})
  maxNfilter <- reactive({input$maxNfilter})
  
  cutLowQualTail <- reactive({switch(input$cutLowQualTail, "no" = FALSE, "yes" = TRUE)})
  cutTailWindowSize <- reactive({input$cutTailWindowSize})
  cutTailMeanQual <- reactive({input$cutTailMeanQual})
  
  cutSlideWindowRight <- reactive({switch(input$cutSlideWindowRight, "no" = FALSE, "yes" = TRUE)})
  cutSlideWindowSize <- reactive({input$cutSlideWindowSize})
  cutSlideWindowQual <- reactive({input$cutSlideWindowQual})
  
  cutLowQualFront <- reactive({switch(input$cutLowQualFront, "no" = FALSE, "yes" = TRUE)})
  cutFrontWindowSize <- reactive({input$cutFrontWindowSize})
  cutFrontMeanQual <- reactive({input$cutFrontMeanQual})
  
  lengthFiltering <- reactive({switch(input$lengthFiltering, "no" = FALSE, "yes" = TRUE)})
  minReadLength <- reactive({input$minReadLength})
  maxReadLength <- reactive({input$maxReadLength})
  
  trimFrontRead1 <- reactive({input$trimFrontRead1})
  trimTailRead1 <- reactive({input$trimTailRead1})  
  trimFrontRead2 <- reactive({input$trimFrontRead2})
  trimTailRead2 <- reactive({input$trimTailRead2})
  maxLengthRead1 <- reactive({input$maxLengthRead1})
  maxLengthRead2 <- reactive({input$maxLengthRead2})
  
  forceTrimPolyG <- reactive({switch(input$forceTrimPolyG, "no" = FALSE, "yes" = TRUE)})
  disableTrimPolyG <- reactive({switch(input$disableTrimPolyG, "no" = FALSE, "yes" = TRUE)})
  minLengthPolyG <- reactive({input$minLengthPolyG})
  trimPolyX <- reactive({switch(input$trimPolyX, "no" = FALSE, "yes" = TRUE)})
  minLengthPolyX <- reactive({input$minLengthPolyX})
  
  overrepresentationAnalysis <- reactive({switch(input$overrepresentationAnalysis, "no" = FALSE, "yes" = TRUE)})
  overrepresentationSampling <- reactive({input$overrepresentationSampling})
  
  correctionOverlap <- reactive({switch(input$correctionOverlap, "no" = FALSE, "yes" = TRUE)})
  minOverlapLength <- reactive({input$minOverlapLength})
  maxOverlapMismatch <- reactive({input$maxOverlapMismatch})
  maxOverlapMismatchPercentage <- reactive({input$maxOverlapMismatchPercentage})
  
  umi <- reactive({switch(input$umi, "no" = FALSE, "yes" = TRUE)})
  umiLoc <- reactive({input$umiLoc})
  umiLength <- reactive({input$umiLength})
  umiSkipBaseLength <- reactive({input$umiSkipBaseLength})
  umiPrefix <- reactive({input$umiPrefix})
  umiNoConnection <- reactive({switch(input$umiNoConnection, "no" = FALSE, "yes" = TRUE)})
  
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
        add_row(setting = "Add _ between prefix & UMI", value = as.character(umiNoConnection()))
    }
    
    return(tbl)
  }, options = list(dom = 'ltipr', paging = FALSE))
  
  ########################################################### Trimming ########################################################### 
  observeEvent(input$runtrimming, {
    i <- match(input$fastqSample_tr, allsamples())

    if (read_type() == 'single') {
      se_dir <- paste0(session$token, "/", input$fastqSample_tr, "/", raw_se_dir()[i])
      se_fastq <- list.files(path = se_dir, full.names = T)
      se_trimmed_dir <- paste0(session$token, "/", input$fastqSample_tr)
      se_names <- gsub(pattern = fastq_suffix()[i], replacement = "", x = list.files(path = se_dir, full.names = F))
      
      for (j in 1:length(se_fastq)) {
        rfastp(read1 = se_fastq[j],
               outputFastq = paste0(se_trimmed_dir, "/", se_names[j], "_trimmed"),
               adapterTrimming = adapterTrimming(),
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
               umiSkipBaseLength = umiSkipBaseLength())
      }
      
      # Create summary file for QuasR
      trimmed <- list.files(path = se_trimmed_dir, pattern = "_trimmed_R1.fastq.gz$", full.names = F)
      df <- data.frame(trimmed, se_names)
      names(df) <- c("FileName", "SampleName")
      write.table(x = df, file = paste0(se_trimmed_dir, "/", "trimmed.txt"), sep="\t", row.names=FALSE, quote = FALSE)
      
    } else {
      r1_dir <- paste0(session$token, "/", input$fastqSample_tr, "/", raw_per1_dir()[i])
      r2_dir <- paste0(session$token, "/", input$fastqSample_tr, "/", raw_per2_dir()[i])
      r1_fastq <- list.files(path = r1_dir, full.names = T)
      r2_fastq <- list.files(path = r2_dir, full.names = T)
      pe_trimmed_dir <- paste0(session$token, "/", input$fastqSample_tr)
      pe_names <- gsub(pattern = fastq_suffix()[i], replacement = "", x = list.files(path = r1_dir, full.names = F))
      
      for (j in 1:length(r1_fastq)) {
        rfastp(read1 = r1_fastq[j],
               read2 = r2_fastq[j],
               outputFastq = paste0(pe_trimmed_dir, "/", pe_names[j], "_trimmed"),
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
               umiSkipBaseLength = umiSkipBaseLength())
      }
      
      # Create summary file for QuasR
      trimmed_r1 <- list.files(path = pe_trimmed_dir, pattern = "_trimmed_R1.fastq.gz$", full.names = T)
      trimmed_r2 <- list.files(path = pe_trimmed_dir, pattern = "_trimmed_R2.fastq.gz$", full.names = T)
      df <- data.frame(trimmed_r1, trimmed_r2, pe_names)
      names(df) <- c("FileName1", "FileName2", "SampleName")
      write.table(x = df, file = paste0(pe_trimmed_dir, "/", "trimmed.txt"), sep="\t", row.names=FALSE, quote = FALSE)
    }
  })
  
  # Plot trimming statistics
  output$select_trim_file <- renderUI({
    json_files <- list.files(path = paste0(session$token, "/", input$fastqSample_tr), pattern = "_trimmed.json$", full.names = F)
    selectInput("trimmed_fastq_file", "Select Trimmed Fastq JSON", choices = json_files)
  })
  
  output$trim_overrepresented <- renderDataTable({
    json <- fromJSON(file = paste0(session$token, "/", input$fastqSample_tr, "/", input$trimmed_fastq_file))
    
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
  
  output$table_trim_summary <- renderTable({
    json <- fromJSON(file = paste0(session$token, "/", input$fastqSample_tr, "/", input$trimmed_fastq_file))
    qcSummary(json)
  }, rownames = TRUE)
  
  output$trim_base_quality <- renderPlot({
    json <- fromJSON(file = paste0(session$token, "/", input$fastqSample_tr, "/", input$trimmed_fastq_file))
    curvePlot(json)
  })
  
  output$trim_gc_content <- renderPlot({
    json <- fromJSON(file = paste0(session$token, "/", input$fastqSample_tr, "/", input$trimmed_fastq_file))
    curvePlot(json, curves = "content_curves")
  })
  
  ######################################################### FASTA ######################################################### 
  output$refseq_organisms <- renderDataTable({
    refseq_table[refseq_table$kingdoms == input$kingdom, ]
  }, filter = "top", selection = "none",options = list(pageLength = 10, scrollX = TRUE))
  
  output$select_organism <- renderUI({
    x <- refseq_table[refseq_table$kingdoms == input$kingdom, ]
    y <- unique(x$organism_name)
    pickerInput("organism", "Select Organism From NCBI", choices = y, options = pickerOptions(liveSearch = T))
  })
  
  output$select_accession <- renderUI({
    x <- refseq_table[refseq_table$kingdoms == input$kingdom & refseq_table$organism_name == input$organism, ]
    pickerInput("accession", "Select Accession No. From NCBI", choices = x$assembly_accession, options = pickerOptions(liveSearch = T))
  })
  
  observeEvent(input$get_fasta, {
    fasta_file <- getGenome(db = "refseq", organism = input$accession, reference = TRUE, path = session$token)
    clickValues$full_fasta <- Biostrings::readDNAStringSet(fasta_file)
    
    output$select_fasta_chrom <- renderUI({
      pickerInput("fasta_chrom", "Select Chromosomes from FASTA File", 
                  choices = names(clickValues$full_fasta), multiple = TRUE, 
                  options = pickerOptions(actionsBox = TRUE, liveSearch = TRUE))
    })
  })
  
  observeEvent(input$load_fasta, {
    full_fasta <- clickValues$full_fasta
    fasta_chrom <- input$fasta_chrom
    Biostrings::writeXStringSet(full_fasta[fasta_chrom],append=FALSE,filepath = paste0(session$token, ".fasta"))
  })
  
  fasta_file <- reactive({
    if (input$fasta_how == "yes") {
      return(input$fasta_file$datapath)
    } else {
      return(paste0(session$token, ".fasta"))
    }
  })
  
  observeEvent(input$get_txdb, {
    gtf_file <- getGFF(db = "refseq", organism = input$accession, reference = TRUE, path = session$token)
    clickValues$full_txdb <- makeTxDbFromGFF(gtf_file)
    
    output$select_txdb_chrom <- renderUI({
      pickerInput("gtf_chrom", "Select Chromosomes from Txdb", 
                  choices = seqlevels(clickValues$full_txdb), multiple = TRUE, 
                  options = pickerOptions(actionsBox = TRUE, liveSearch = TRUE))
    })
  })
  
  observeEvent(input$load_txdb, {
    seqlevels(clickValues$full_txdb) <- input$gtf_chrom
  })
  
  txdb <- reactive({
    if (input$gtf_how == "yes") {
      x <- makeTxDbFromGFF(input$gtf_file$datapath)
      return(x)
    } else {
      return(clickValues$full_txdb)
    }
  })
  
  #################################################### Aligning #################################################### 

  observeEvent(input$runalign, {
    dir.create(paste0(session$token, "/", input$fastqSample_al, "/", "alignments"), recursive = T)
    proj <- qAlign(sampleFile = paste0(session$token, "/", input$fastqSample_al, "/", "trimmed.txt"),
                   genome = fasta_file(),
                   aligner = "Rhisat2",
                   alignmentsDir = paste0(session$token, "/", input$fastqSample_al, "/", "alignments"))
    
    qCount(proj, txdb())
  })
}

shinyApp(ui, server)


