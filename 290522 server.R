.check_if_unix <<- function() {
  return(NULL)
}
assignInNamespace(".check_if_unix", .check_if_unix, ns = "fastqcr")

options(shiny.maxRequestSize = 5000000*1024^2)

server <- function(input, output, session) {
  dir.create(paste0(session$token, "/alignments"), recursive = TRUE)
  
  clickValues <- reactiveValues(sample_name = NULL, raw_fastq = NULL, trimmed_fastq = NULL,
                                full_fasta = NULL, full_txdb = NULL, full_qc = NULL, dxTrack = NULL)
  
  ########################################### FASTQC #############################################
  fastqc_dir <- reactive({paste0(session$token, "/fastqc")})
  
  # Paired-end files
  raw_fastq_pe_r1_path <- reactive({
    req(input$raw_fastq_pe_r1)
    input$raw_fastq_pe_r1$datapath
  })
  raw_fastq_pe_r1_files <- reactive({
    req(input$raw_fastq_pe_r1)
    unzip(raw_fastq_pe_r1_path(), list = TRUE)$Name
  })
  raw_per1_dir <- reactive(tools::file_path_sans_ext(input$raw_fastq_pe_r1$name))

  raw_fastq_pe_r2_path <- reactive({
    req(input$raw_fastq_pe_r2)
    input$raw_fastq_pe_r2$datapath
  })
  raw_fastq_pe_r2_files <- reactive({
    req(input$raw_fastq_pe_r2)
    unzip(raw_fastq_pe_r2_path(), list = TRUE)$Name
  })
  raw_per2_dir <- reactive(tools::file_path_sans_ext(input$raw_fastq_pe_r2$name))

  # Single-end files
  raw_fastq_se_path <- reactive({
    req(input$raw_fastq_se)
    input$raw_fastq_se$datapath
  })
  raw_fastq_se_files <- reactive({
    req(input$raw_fastq_se)
    unzip(raw_fastq_se_path(), list = TRUE)$Name
  })
  raw_se_dir <- reactive(tools::file_path_sans_ext(input$raw_fastq_se$name))

  suffix <- reactive({input$suffix})
  
  observeEvent(input$run_fastqc, {
    showModal(modalDialog("Loading FASTQC ... ", footer=NULL))
    if (!dir.exists("FastQC")) {
      fastqcr::fastqc_install(dest.dir = ".")
    }
    
    if (input$read_type == "Single") {
      unzip(raw_fastq_se_path(), exdir = session$token)
      fastqc(fq.dir = paste0(session$token, "/", raw_se_dir()), qc.dir = fastqc_dir(),
             fastqc.path = "wsl ./FastQC/fastqc")
      
      sample_name <- gsub(pattern = suffix(), replacement = "", x=basename(unlist(raw_fastq_se_files())))
      df <- data.frame(raw_fastq_se_files(), sort(sample_name))
      names(df) <- c("FileName", "SampleName")
      write.table(df[c("FileName", "SampleName")], file = paste0(session$token, sep = "/", "untrimmed.txt"), sep="\t",row.names=FALSE, quote = FALSE)
      
    } else {
      unzip(raw_fastq_pe_r1_path(), exdir = session$token)
      unzip(raw_fastq_pe_r2_path(), exdir = session$token)
      fastqc(fq.dir = paste0(session$token, "/", raw_per1_dir()), qc.dir = fastqc_dir(),
             fastqc.path = "wsl ./FastQC/fastqc")
      fastqc(fq.dir = paste0(session$token, "/", raw_per2_dir()), qc.dir = fastqc_dir(),
             fastqc.path = "wsl ./FastQC/fastqc")
      
      sample_name <- gsub(pattern = suffix(), replacement = "", basename(unlist(raw_fastq_pe_r1_files())))
      df <- data.frame(raw_fastq_pe_r1_files(), raw_fastq_pe_r2_files(), sort(sample_name))
      names(df) <- c("FileName1", "FileName2", "SampleName")
      write.table(df[c("FileName1", "FileName2", "SampleName")], file = paste0(session$token, sep = "/", "untrimmed.txt"), sep="\t",row.names=FALSE, quote = FALSE)
    }
    
    clickValues$sample_name <- sort(sample_name)
    clickValues$raw_fastq <- df
    removeModal()
  })
  
  fastqc_files <- reactive({
    req(input$run_fastqc)
    list.files(fastqc_dir(), pattern = "fastqc.zip$", full.names = TRUE)
  })
  
  output$select_raw_fastq <- renderUI({
    req(input$run_fastqc)
    selectInput("raw_fastq_file", "Select FASTQ File", choices = basename(fastqc_files()))
  }) 
  
  output$select_trimmed_fastq <- renderUI({
    req(input$run_fastqc)
    selectInput("trimmed_fastq_file", "Select trimmed FASTQ File", choices = clickValues$sample_name)
  })
  
  fdl <- reactive({FastqcDataList(fastqc_files())})
  
  output$multiqc_plot <- renderPlot({
    req(input$run_fastqc, input$multiqc_plot_type)
    fastqc_plot(fdl(), input$multiqc_plot_type)
  })
  
  output$singleqc_plot <- renderPlot({
    req(input$run_fastqc, input$raw_fastq_file, input$singleqc_plot_type)
    which_fdl <- fdl()[[toString(fastqc_files()[which(grepl(input$raw_fastq_file, fastqc_files()))])]]
    fastqc_plot(which_fdl, input$singleqc_plot_type)
  })
  
  
  ########################################### Configure Trim Settings #############################################
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
  
  output$table_ts <- DT::renderDataTable({
    tbl <- tibble("setting" = character(), "value" = character())
    
    tbl <- tbl %>% add_row(setting = "Adapter Trimmming", value = as.character(adapterTrimming()))
    if (adapterTrimming()) {
      tbl <- tbl %>% add_row(setting = "Adapter(R1)", value = adapterSequenceRead1())
      if (read_type() == "Paired") {
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
    if (read_type() == "Paired") {
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
    
    return(tbl)
  }, options = list(dom = 'ltipr', paging = FALSE))
  
  ################################################# RUN TRIMMING ################################################# 
  observeEvent(input$trimming, {
    show_modal_progress_line(
      text = "Start Trimming...",
      color = "#73C6B6"
    )
    
    df <- clickValues$raw_fastq
    
    if (read_type()=='Single') {
      for (i in 1:nrow(df)) {
        rfastp(read1 = paste0(session$token, "/", df$FileName[i]),
               outputFastq = paste0(session$token, "/",df$SampleName[i], "_trimmed"),
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
               overrepresentationSampling = overrepresentationSampling())
        
        update_modal_progress(i / nrow(df), text = paste0("Trimming Completed For: ", i, " Samples"))
      }
      
      # This is for Quasr txt file
      trimmed <- list.files(path = session$token, pattern = "_trimmed_R1.fastq.gz$", full.names = FALSE)
      json_summary <- list.files(path = session$token, pattern = "_trimmed.json$", full.names = TRUE)
      df <- data.frame(clickValues$sample_name, trimmed, json_summary)
      names(df) <- c("SampleName", "FileName", "Summary")
      write.table(df[c("FileName", "SampleName")], file = paste0(session$token, sep = "/", "trimmed.txt"), sep="\t",row.names=FALSE, quote = FALSE)

    } else {
      for (i in 1:nrow(df)) {
        rfastp(read1 = paste0(session$token, "/", df$FileName1[i]),
               read2 = paste0(session$token, "/", df$FileName2[i]),
               outputFastq = paste0(session$token, "/", df$SampleName[i], "_trimmed"),
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
               overrepresentationSampling = overrepresentationSampling())
        
        update_modal_progress(i / nrow(df), text = paste0("Trimming Completed For: ", i, " Samples"))
      }
      
      # This is for Quasr txt file
      trimmed_1 <- list.files(path = session$token, pattern = "_trimmed_R1.fastq.gz$", full.names = FALSE)
      trimmed_2 <- list.files(path = session$token, pattern = "_trimmed_R2.fastq.gz$", full.names = FALSE)
      json_summary <- list.files(path = session$token, pattern = "_trimmed.json$", full.names = TRUE)
      df <- data.frame(clickValues$sample_name, trimmed_1, trimmed_2, json_summary)
      names(df) <- c("SampleName", "FileName1", "FileName2", "Summary")
      write.table(df[c("FileName1", "FileName2", "SampleName")], file = paste0(session$token, sep = "/", "trimmed.txt"), sep="\t",row.names=FALSE, quote = FALSE)
      
    }
    
    clickValues$trimmed_fastq <- df
    remove_modal_progress()
  })
  
  #PLOT TRIM STATISTICS
  output$table_trim_summary <- renderTable({
    req(input$trimming, input$trimmed_fastq_file)
    df <- clickValues$trimmed_fastq
    json <- fromJSON(file = df[df$SampleName == input$trimmed_fastq_file, ]$Summary)
    qcSummary(json)
  }, rownames = TRUE)
  
  output$raw_overrepresented <- renderDataTable({
    req(input$trimming, input$trimmed_fastq_file)
    df <- clickValues$trimmed_fastq
    json <- fromJSON(file = df[df$SampleName == input$trimmed_fastq_file, ]$Summary)
    
    if (read_type() == "Single") {
      df <- data.frame(unlist(json$read1_before_filtering$overrepresented_sequences))
      names(df) <- c("Count")
    } else {
      df1 <- data.frame(unlist(json$read1_before_filtering$overrepresented_sequences))
      names(df1) <- c("Count")
      df1$Read <- "R1"
      
      df2 <- data.frame(unlist(json$read2_before_filtering$overrepresented_sequences))
      names(df2) <- c("Count")
      df2$Read <- "R2"
      
      df <- rbind(df1, df2)
    }
    return(df)
  }, options = list(scrollX = TRUE, scrollY = TRUE))
  
  output$trim_overrepresented <- renderDataTable({
    req(input$trimming, input$trimmed_fastq_file)
    df <- clickValues$trimmed_fastq
    json <- fromJSON(file = df[df$SampleName == input$trimmed_fastq_file, ]$Summary)
    
    if (read_type() == "Single") {
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
  
  output$trim_base_quality <- renderPlot({
    req(input$trimming, input$trimmed_fastq_file)
    df <- clickValues$trimmed_fastq
    json <- fromJSON(file = df[df$SampleName == input$trimmed_fastq_file, ]$Summary)
    curvePlot(json)
  })
  
  output$trim_gc_content <- renderPlot({
    req(input$trimming, input$trimmed_fastq_file)
    df <- clickValues$trimmed_fastq
    json <- fromJSON(file = df[df$SampleName == input$trimmed_fastq_file, ]$Summary)
    curvePlot(json, curves = "content_curves")
  })  

  ##################################### FASTA File #####################################
  output$ensembl_table <- renderDataTable(ensembl_table, filter = "top", selection = "none",options = list(pageLength = 10, scrollX = TRUE))
  output$select_ensembl <- renderUI({pickerInput("ensembl", "Select Ensembl Database", choices = ensembl_table$name, options = list(`live-search` = TRUE))})
  
  ensembl <- reactive({input$ensembl})
  assembly_type <- reactive({input$assembly_type})
  
  observeEvent(input$ensembl_fasta, {
    showModal(modalDialog("Loading FASTA From Ensembl ... ", footer=NULL))
    all_fasta <- getGenome(db = "ensembl",
                           organism = ensembl(),
                           release = 106,
                           path = session$token,
                           assembly_type = assembly_type())
    clickValues$full_fasta <- Biostrings::readDNAStringSet(all_fasta)
    output$select_fasta_chrom <- renderUI({
      pickerInput("fasta_chrom", "Select Chromosomes for FASTA File", 
                  choices = names(clickValues$full_fasta), multiple = TRUE, 
                  options = pickerOptions(actionsBox = TRUE, liveSearch = TRUE))
    })
    removeModal()
    updateActionButton(session, "ensembl_fasta", label = "Done Loading!", icon = icon("check"))
  })
  
  observeEvent(input$load_fasta, {
    full_fasta <- clickValues$full_fasta
    fasta_chrom <- input$fasta_chrom
    Biostrings::writeXStringSet(full_fasta[fasta_chrom],append=FALSE,filepath = paste0(session$token, ".fasta"))
    show_alert(title = "Generate Fasta", text = "Fasta File from Ensembl Successfully Created!", type = "success", btn_colors = "#73C6B6")
    updateActionButton(session, "load_fasta", label = "Fasta Created!", icon = icon("check"))
  })
  
  fasta_file <- reactive({
    if (input$fasta_how == "Upload") {
      return(input$fasta_file$datapath)
    } else {
      return(paste0(session$token, ".fasta"))
    }
  })
  
  ##################################### GTF File #####################################
  observeEvent(input$ensembl_txdb, {
    showModal(modalDialog("Creating TxDb From Selected Chromosomes ... ", footer=NULL))
    clickValues$full_txdb <- makeTxDbFromEnsembl(organism = ensembl())
    output$select_gtf_chrom <- renderUI({
      pickerInput("gtf_chrom", "Select Chromosomes for GTF File",
                  choices = seqlevels(clickValues$full_txdb), multiple = TRUE, 
                  options = pickerOptions(actionsBox = TRUE, liveSearch = TRUE))
    })
    removeModal()
    updateActionButton(session, "ensembl_txdb", label = "Done Loading!", icon = icon("check"))
  })
  
  observeEvent(input$load_gtf, {
    seqlevels(clickValues$full_txdb) <- input$gtf_chrom
    show_alert(title = "Generate GTF", text = "TxDb from Ensembl Successfully Created!", type = "success", btn_colors = "#73C6B6")
    updateActionButton(session, "load_gtf", label = "Txdb Created!", icon = icon("check"))
  })
  
  txdb <- reactive({
    if (input$gtf_how == "Upload") {
      x <- makeTxDbFromGFF(input$gtf_file$datapath)
      return(x)
    } else {
      return(clickValues$full_txdb)
    }
  })
  
  ##################################### Alignment #####################################
  sampleFile <- reactive({
    if (input$fastq_align) {
      return(paste0(session$token, sep = "/", "trimmed.txt"))
    } else {
      return(paste0(session$token, sep = "/", "untrimmed.txt"))
    }
  })
  paired <- reactive({input$paired})
  splicedAlignment <- reactive({input$splicedAlignment})
  aligner <- reactive({input$aligner})
  maxHits <- reactive({input$maxHits})
  selectReadPosition <- reactive({input$selectReadPosition})
  orientation <- reactive({input$orientation})
  useRead <- reactive({input$useRead})
  includeSpliced <- reactive({input$includeSpliced})
  includeSecondary <- reactive({input$includeSecondary})
  mapqMin <- reactive({input$mapqlty[1]})
  mapqMax <- reactive({input$mapqlty[2]})
  
  observeEvent(input$align, {
    showModal(modalDialog("Aligning ... ", footer=NULL))
    extract_splice_sites(txdb(), outfile = paste0(session$token, "/", "splice_site.txt"))
    
    if (splicedAlignment() & aligner() == "Rhisat2") {
      proj <- qAlign(sampleFile = sampleFile(),
                     genome = fasta_file(),
                     paired = paired(),
                     splicedAlignment = splicedAlignment(),
                     aligner = aligner(),
                     maxHits = maxHits(),
                     alignmentParameter = paste0("--known-splicesite-infile ", session$token, "/", "splice_site.txt"),
                     alignmentsDir = paste0(session$token, "/alignments"))
    } else {
      proj <- qAlign(sampleFile = sampleFile(),
                     genome = fasta_file(),
                     paired = paired(),
                     splicedAlignment = splicedAlignment(),
                     maxHits = maxHits(),
                     aligner = aligner(),
                     alignmentsDir = paste0(session$token, "/alignments"))
    }
    removeModal()
    
    # Feature counts: Gene level
    output$gene_counts <- renderDataTable(server = FALSE,{
      datatable(
        qCount(proj, 
               txdb(), 
               reportLevel = "gene",
               selectReadPosition = selectReadPosition(),
               orientation = orientation(),
               useRead = useRead(),
               includeSpliced = includeSpliced(),
               includeSecondary = includeSecondary(),
               mapqMin = mapqMin(),
               mapqMax = mapqMax()), 
        rownames = TRUE,
        extensions = 'Buttons',
        options = list(paging = TRUE, searching = TRUE, fixedColumns = TRUE, 
                       autoWidth = TRUE,  ordering = TRUE, pageLength=15,
                       dom = 'Bfrtip', 
                       buttons = 
                         list('copy', 'print', list(
                           extend = 'collection',
                           buttons = c('csv', 'excel', 'pdf'),
                           text = 'Download'
                         )))
        )
    })
    
    # Feature counts: Exon level
    output$exon_counts <- renderDataTable(server = FALSE,{
      datatable(
        qCount(proj, 
               txdb(), 
               reportLevel = "exon",
               selectReadPosition = selectReadPosition(),
               orientation = orientation(),
               useRead = useRead(),
               includeSpliced = includeSpliced(),
               includeSecondary = includeSecondary(),
               mapqMin = mapqMin(),
               mapqMax = mapqMax()), 
        rownames = TRUE,
        extensions = 'Buttons',
        options = list(paging = TRUE, searching = TRUE, fixedColumns = TRUE, 
                       autoWidth = TRUE,  ordering = TRUE, pageLength=15,
                       dom = 'Bfrtip', 
                       buttons = 
                         list('copy', 'print', list(
                           extend = 'collection',
                           buttons = c('csv', 'excel', 'pdf'),
                           text = 'Download'
                         )))
      )
    })
    
    # Feature counts: Junction level
    output$junction_counts <- renderDataTable(server = FALSE,{
      datatable(
        data.frame(
          qCount(proj, 
                 query = NULL,
                 reportLevel = "junction",
                 selectReadPosition = selectReadPosition(),
                 orientation = orientation(),
                 useRead = useRead(),
                 includeSpliced = includeSpliced(),
                 includeSecondary = includeSecondary(),
                 mapqMin = mapqMin(),
                 mapqMax = mapqMax())
        ),
        rownames = TRUE,
        extensions = 'Buttons',
        options = list(paging = TRUE, searching = TRUE, fixedColumns = TRUE, 
                       autoWidth = TRUE,  ordering = TRUE, pageLength=15,
                       dom = 'Bfrtip', 
                       buttons = 
                         list('copy', 'print', list(
                           extend = 'collection',
                           buttons = c('csv', 'excel', 'pdf'),
                           text = 'Download'
                         )))
      )
    })
    
    # Alignment Statistics
    clickValues$full_qc <- qQCReport(proj, pdfFilename = paste0(session$token, "/", "qc_alignment.pdf"), useSampleNames = TRUE)
    
    # Plot Alignment Coverage
    qExportWig(proj)
    wig_file <- list.files(pattern = ".wig.gz$")
    hist_col <- colorspace::qualitative_hcl(n = length(wig_file), palette = "Dark 3")
    
    dxTrack <- c()
    for (i in 1:length(wig_file)) {
      x <- DataTrack(range = wig_file[i], 
                     name = gsub(pattern = ".wig.gz", replacement = "", x = basename(wig_file[i])),
                     type = "h",
                     col = hist_col[i])
      dxTrack <- c(dxTrack, x)
    }
    clickValues$dxTrack <- dxTrack
  })
  
  output$align_stat_plot <- renderPlot({align_plot(clickValues$full_qc, input$align_stat_type)})
  
  ##################################### View Coverage #####################################
  output$select_plot_chrom <- renderUI({
    available_chrom <- seqlevels(txdb())
    pickerInput("plot_cov_chrom", "Select Chromosome To View Coverage", choices = available_chrom)
  })
  plot_cov_chrom <- reactive({input$plot_cov_chrom})
  
  output$select_plot_range <- renderUI({
    if (input$gtf_how == "Upload") {
      genes <- genes(txdb(),filter=list(tx_chrom = plot_cov_chrom()))
      df <- data.frame(genes)
      min_range <- min(df$start)
      max_range <- max(df$end)
    } else {
      min_range <- 1
      max_range <- seqlengths(txdb())[[plot_cov_chrom()]]
    }
    q1 <- floor((max_range - min_range)/4 + min_range)
    q2 <- floor((max_range - min_range)/2 + min_range)
    sliderInput("plot_cov_range", "Range to plot", min = min_range, max = max_range, value = c(q1, q2))
  })
  plot_cov_from <- reactive({input$plot_cov_range[1]})
  plot_cov_to <- reactive({input$plot_cov_range[2]})
  
  output$selected_ranges <- renderDT({
    gr <- GRanges(seqnames = plot_cov_chrom(),
                  ranges = IRanges(start = plot_cov_from(), plot_cov_to()))
    
    x <- switch (input$sel_ranges_type,
                 "gene" = genes(txdb()),
                 "transcript" = transcripts(txdb()),
                 "exon" = exons(txdb())
    )
    
    return(data.frame(subsetByOverlaps(x, gr)))
  }, filter = "top", selection = "none", options = list(pageLength = 15))
  
  extend_right <- reactive({input$extend_right})
  extend_left <- reactive({input$extend_left})
  
  txTrack_tr <- reactive({GeneRegionTrack(txdb(), name = "Transcripts", showId = TRUE, 
                                          background.title = "brown", background.panel = "#FFFEDB")})
  txTrack_ge <- reactive({
    txTrack <- GeneRegionTrack(txdb(), name = "Genes", showId = TRUE, collapseTranscripts = TRUE,
                               background.title = "darkblue", background.panel = "#FFFEDB")
    symbol(txTrack) <- gene(txTrack)
    return(txTrack)
  })
  
  observeEvent(input$view_cov, {
    axTrack <- GenomeAxisTrack()
    
    output$cov_plot <- renderPlot({
      if (input$cov_plot_type == "Transcripts") {
        plotTracks(as.list(c(axTrack, clickValues$dxTrack, txTrack_tr())),
                   chromosome = plot_cov_chrom(),
                   from = plot_cov_from(), 
                   to = plot_cov_to(),
                   extend.right = extend_right(),
                   extend.left = extend_left())
      } else {
        plotTracks(as.list(c(axTrack, clickValues$dxTrack, txTrack_ge())),
                   chromosome = plot_cov_chrom(),
                   from = plot_cov_from(), 
                   to = plot_cov_to(),
                   extend.right = extend_right(),
                   extend.left = extend_left())
      }
    })
  })
  
  ##################################### Alternative Splicing #####################################
  output$download_targets <- downloadHandler(
    filename = function() {
      paste("targets-", Sys.Date(), ".txt", sep="")
    },
    
    content = function(file) {
      sample_names <- clickValues$sample_name
      bam_files <- list.files(paste0(session$token, "/alignments"), pattern = ".bam$", full.names = TRUE)
      targets <- data.frame(row.names = sample_names, bam = bam_files)
      write.table(x = targets, append = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    }
  )
  
  targets <- reactive({input$targets_file})
  minReadLength <- reactive({input$minReadLength})
  maxISize <- reactive({input$maxISize})
  minAnchor <- reactive({input$minAnchor})
  threshold <- reactive({input$threshold})
  
}

shinyApp(ui, server)



