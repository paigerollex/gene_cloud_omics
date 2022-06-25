ui <- tagList(
  shinyjs::useShinyjs(),
  useShinydashboard(),
  navbarPage(id = "navbar", theme = shinytheme("flatly"), title = "",
    navbarMenu("Quality Control",
      ######################################### FASTQC ########################################################
       tabPanel("FastQC",
         sidebarLayout(
           sidebarPanel(
             tabsetPanel(
               tabPanel("Upload",
                  div(style = "margin-top: +20px"),
                  radioButtons("read_type", "Type of FASTQ Read", choices = c("Single", "Paired"), inline = TRUE),
                  
                  conditionalPanel(condition = "input.read_type == 'Single'",
                                   fileInput("raw_fastq_se", "Choose Fastq Zip File", accept = ".zip")),
                  
                  conditionalPanel(condition = "input.read_type == 'Paired'",
                                   fileInput("raw_fastq_pe_r1", "Choose Fastq Zip File (Read 1)", accept = ".zip"),
                                   fileInput("raw_fastq_pe_r2", "Choose Fastq Zip File (Read 2)", accept = ".zip")),
                  
                  textInput("suffix", "Fastq Suffix"),
                  div(style = "margin-top: -10px"),
                  helpText(HTML("<font size=-1>Example: .fastq.gz/ _1.fastq.gz / _R1.fastq.gz</font>")),
                  
                  div(style = "margin-top: +25px"),
                  actionButton("run_fastqc", "Run Fastqc")
               ),
               
               tabPanel("View Plot",
                  div(style = "margin-top: +20px"),
                  radioButtons("fqc", "Type of FASTQC", choices = c("Multi", "Single"), inline = TRUE),
                  
                  conditionalPanel(condition = "input.fqc == 'Multi'",
                                   radioButtons("multiqc_plot_type", "Type of Multi-QC Plot", choices = multiqc.plot.type)),
                  
                  conditionalPanel(condition = "input.fqc == 'Single'",
                                   uiOutput("select_raw_fastq"),
                                   radioButtons("singleqc_plot_type", "Type of Single-QC Plot", choices = singleqc.plot.type))
               )
             )
           ),

           mainPanel( h3("FASTQC Plots"),
             tabsetPanel(
               tabPanel("Multi-QC",
                 conditionalPanel(condition = "input.multiqc_plot_type == 'Per Base Sequence Quality'",
                                  h3(strong("Per Base Sequence Quality")),
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
                                  )),

                 conditionalPanel(condition = "input.multiqc_plot_type == 'Mean Sequence Quality Per Read'",
                                  h3(strong("Per Sequence Quality Scores")),
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
                                  )),

                 conditionalPanel(condition = "input.multiqc_plot_type == 'Per Base Sequence Content'",
                                  h3(strong("Per Base Sequence Content")),
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
                                  )),

                 conditionalPanel(condition = "input.multiqc_plot_type == 'GC Content'",
                                  h3(strong("Per Sequence GC Content")),
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
                                  )),

                 conditionalPanel(condition = "input.multiqc_plot_type == 'N Content'",
                                  h3(strong("Per base N content")),
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
                                  )),

                 conditionalPanel(condition = "input.multiqc_plot_type == 'Sequence Duplication Levels'",
                                  h3(strong("Sequence Duplication Levels")),
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
                                  )),

                 conditionalPanel(condition = "input.multiqc_plot_type == 'Adapter Content'",
                                  h3(strong("Adapter Content")),
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
                                  )),

                 conditionalPanel(condition = "input.multiqc_plot_type == 'Overrepresented Sequences'",
                                  h3(strong("Overrepresented Sequences")),
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
                                  )),

                 div(style = "margin-top: +50px"),

                 plotOutput("multiqc_plot")
               ),

               tabPanel("Single-QC",
                conditionalPanel(condition = "input.singleqc_plot_type == 'Per Base Sequence Quality'",
                                 h3(strong("Per Base Sequence Quality")),
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
                                 )),

                conditionalPanel(condition = "input.singleqc_plot_type == 'Mean Sequence Quality Per Read'",
                                 h3(strong("Per Sequence Quality Scores")),
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
                                 )),

                conditionalPanel(condition = "input.singleqc_plot_type == 'Per Base Sequence Content'",
                                 h3(strong("Per Base Sequence Content")),
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
                                 )),

                conditionalPanel(condition = "input.singleqc_plot_type == 'GC Content'",
                                 h3(strong("Per Sequence GC Content")),
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
                                 )),

                conditionalPanel(condition = "input.singleqc_plot_type == 'N Content'",
                                 h3(strong("Per base N content")),
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
                                 )),

                conditionalPanel(condition = "input.singleqc_plot_type == 'Sequence Duplication Levels'",
                                 h3(strong("Sequence Duplication Levels")),
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
                                 )),

                conditionalPanel(condition = "input.singleqc_plot_type == 'Adapter Content'",
                                 h3(strong("Adapter Content")),
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
                                 )),

                conditionalPanel(condition = "input.singleqc_plot_type == 'Overrepresented Sequences'",
                                 h3(strong("Overrepresented Sequences")),
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
                                 )),

                 div(style = "margin-top: +50px"),

                 plotOutput("singleqc_plot")
               )
             )
           )
         )
       ),
      
      ########################################### Trimming #############################################
      tabPanel("Trimming",
        sidebarLayout(
          sidebarPanel(
            h4(strong("Fastp Trim Settings")),
            div(style = "margin-top: +30px"),

            div(dataTableOutput("table_ts"), style = "font-size:80%; height:400px; overflow-y: scroll"),
            div(style = "margin-top: +30px"),

            actionButton("trimming", "Run Trimming")
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

                                       conditionalPanel(condition = "input.read_type == 'Paired'",
                                                        radioButtons("adapter2", "Adapter detection (R2)", choices = c("auto", "manual"), inline = TRUE),
                                                        conditionalPanel(condition = "input.adapter2 == 'manual'",
                                                                         textInput("adapter2_manual", "Adapter Sequence (R2)", value = ""),
                                                                         div(style = "margin-top: -15px"),
                                                                         helpText(HTML("<font size=-1>Leave blank if adapter sequence for R1 & R2 are the same.</font>")))))
                    ),

                    box(width = 4,
                      title = "Quality", status = "primary", solidHeader = TRUE,

                      radioButtons("qualityFiltering", "Filter by quality", choices = c("yes", "no"), inline = TRUE),

                      conditionalPanel(condition = "input.qualityFiltering == 'yes'",
                                       numericInput("qualityFilterPhred", "Min base quality score", value = 15),
                                       div(style = "margin-top: -15px"),
                                       helpText(HTML("<font size=-1>Min quality phred score that a base is qualified.</font>")),

                                       div(style = "margin-top: +20px"),

                                       numericInput("qualityFilterPercent", "Max % of bases < Min quality", value = 40),
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

                      numericInput("maxNfilter", "Max N", value = 5),
                      div(style = "margin-top: -15px"),
                      helpText(HTML("<font size=-1>Max no. of N allowed in the sequence, else discarded.</font>"))

                    ),

                    tabBox(width = 4,
                      title = "Per Read Cutting", side = "right", selected = "Enable",

                      tabPanel("Settings", value = "Settings",
                               conditionalPanel(condition = "input.cutLowQualTail == 'yes'",
                                                numericInput("cutTailWindowSize", "Window size: Cut tail", value = 4),
                                                numericInput("cutTailMeanQual", "Mean quality: Cut tail", value = 20)),

                               conditionalPanel(condition = "input.cutSlideWindowRight == 'yes'",
                                                numericInput("cutSlideWindowSize", "Window size: Cut right", value = 4),
                                                numericInput("cutSlideWindowQual", "Mean quality: Cut right", value = 20)),

                               conditionalPanel(condition = "input.cutLowQualFront == 'yes'",
                                                numericInput("cutFrontWindowSize", "Window size: Cut front", value = 4),
                                                numericInput("cutFrontMeanQual", "Mean quality: Cut front", value = 20))
                      ),

                      tabPanel("Enable", value = "Enable",
                        radioButtons("cutLowQualTail", "Cut tail", choices = c("no", "yes"), inline = TRUE),
                        div(style = "margin-top: -15px"),
                        helpText(HTML("<font size=-1>
                                        Move a sliding window from tail (3') to the front, drop the bases in the window if its
                                        mean quality < threshold, stop otherwise. This is similar to Trimmomatic TRAILING method.
                                      </font>")),

                        div(style = "margin-top: +30px"),

                        radioButtons("cutSlideWindowRight", "Cut right", choices = c("no", "yes"), inline = TRUE),
                        div(style = "margin-top: -15px"),
                        helpText(HTML("<font size=-1>
                                        Move a sliding window from front to tail, if meet one window with mean quality < threshold,
                                        drop the bases in the window and the right part, and then stop. This is simliar to Trimmomatic
                                        SLIDINGWINDOW method.
                                      </font>")),

                        div(style = "margin-top: +30px"),

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

                      conditionalPanel(condition = "input.read_type == 'Paired'",
                                       numericInput("trimFrontRead2", "Trim front R2", value = 0),
                                       numericInput("trimTailRead2", "Trim tail R2", value = 0)),

                      div(style = "margin-top: -15px"),
                      helpText(HTML("<font size=-1>No. of bases to hard trim from the front/tail.</font>")),
                      div(style = "margin-top: +30px"),

                      numericInput("maxLengthRead1", "Max length (R1 for PE)", value = 0),

                      conditionalPanel(condition = "input.read_type == 'Paired'",
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

                      div(style = "margin-top: +30px"),

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

                      conditionalPanel(condition = "input.read_type == 'Single'",
                                       helpText("Only for PE reads")),

                      conditionalPanel(condition = "input.read_type == 'Paired'",
                                       radioButtons("correctionOverlap", "Run correction overlap", choices = c("no", "yes"), inline = TRUE),
                                       div(style = "margin-top: -15px"),
                                       helpText(HTML('<font size=-1>
                                                        If a proper overlap between a pair of reads is found, it can correct mismatched base pairs or if 1 base is with
                                                        high quality while the other is with ultra low quality. If a base is corrected, the quality of its paired base
                                                        will be assigned to it so that they will share the same quality.
                                                     </font>')),

                                       conditionalPanel(condition = "input.correctionOverlap == 'yes'",
                                                        numericInput("minOverlapLength", "Min overlap length", value = 30),
                                                        div(style = "margin-top: -15px"),
                                                        helpText(HTML("<font size=-1>Min length to detect overlapped region of PE reads.</font>")),

                                                        numericInput("maxOverlapMismatch", "Max mismatched bases", value = 5),
                                                        div(style = "margin-top: -15px"),
                                                        helpText(HTML("<font size=-1>Max number of mismatched bases to detect overlapped region of PE reads.</font>")),

                                                        numericInput("maxOverlapMismatchPercentage", "Max mismatch % in overlap", value = 20),
                                                        div(style = "margin-top: -15px"),
                                                        helpText(HTML("<font size=-1>Max percentage of mismatched bases to detect overlapped region of PE reads.</font>"))))

                    ),

                    box(width = 4,
                      title = "UMI", status = "warning", solidHeader = TRUE,
                      radioButtons("umi", "Preprocess UMI", choices = c("no", "yes"), inline = TRUE),

                      conditionalPanel(condition = "input.umi == 'yes'",
                                       selectInput("umiLoc", "Location of UMI", choices = c("read1", "read2", "index1", "index2")),
                                       div(style = "margin-top: -15px"),
                                       helpText(HTML("<font size=-1>
                                                        <li>read1: Head of read1 is used as UMI. If PE, used for both reads.</li>
                                                        <li>read2: Head of read2 is used as UMI. If PE, used for both reads.</li>
                                                     </font>")),

                                       div(style = "margin-top: +30px"),

                                       conditionalPanel(condition = "input.umiLoc == 'read1'| input.umiLoc == 'read2'",
                                                        numericInput("umiLength", "Length of UMI in read1/read2", value = 0),
                                                        numericInput("umiSkipBaseLength", "No. of bases to skip", value = 0),
                                                        div(style = "margin-top: -15px"),
                                                        helpText(HTML("<font size=-1>Skip n bases after UMI to trim the UMI separator and A/T tailing in read1/read2.</font>"))),

                                       conditionalPanel(condition = "input.umiLoc == 'index1'| input.umiLoc == 'index2'",
                                                        fileInput("index1Filter", "Filter Index 1"),
                                                        div(style = "margin-top: -15px"),
                                                        fileInput("index2Filter", "Filter Index 2"),
                                                        div(style = "margin-top: -25px"),
                                                        helpText(HTML("<font size=-1>A file which contains a list of barcodes of Index 1/2 to be filtered out, one barcode per line.</font>"))),

                                       textInput("umiPrefix", "UMI Prefix", value = ""),
                                       div(style = "margin-top: -15px"),
                                       helpText(HTML("<font size=-1>String prefiz to label UMI sequence.</font>")))
                    )
                  )
                )
              ),

              tabPanel("Results",
                h3(strong("Trim Result")),
                uiOutput("select_trimmed_fastq"),

                dashboardSidebar(disable = TRUE),
                dashboardBody(
                  fluidRow(
                    box(width = 6, height = "400px",
                      title = "Summary", status = "primary", solidHeader = TRUE,
                      tableOutput("table_trim_summary")
                    ),
                  
                    tabBox(width = 6, side = "right", height = "400px",
                      title = "Overrepresented Sequences",
                       tabPanel("After",
                                conditionalPanel(condition = "input.overrepresentationAnalysis == 'yes'",
                                  div(dataTableOutput("trim_overrepresented"), style = "font-size:80%; height:330px; overflow-y: scroll; overflow-x: scroll")
                                ),
                                
                                conditionalPanel(condition = "input.overrepresentationAnalysis == 'yes'",
                                  helpText(HTML("Only when conducting over-representation analysis."))
                                )
                       ),
                       
                       tabPanel("Before",
                                conditionalPanel(condition = "input.overrepresentationAnalysis == 'yes'",
                                  div(dataTableOutput("raw_overrepresented"), style = "font-size:80%; height:330px; overflow-y: scroll; overflow-x: scroll")
                                ),
                                
                                conditionalPanel(condition = "input.overrepresentationAnalysis == 'yes'",
                                                 helpText(HTML("Only when conducting over-representation analysis."))
                                )
                       )                 
                    )

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
    
    # 2nd nav menu
    navbarMenu("Aligning",
               
      ########################################### Aligning #############################################
      tabPanel("Generate Counts",
        sidebarLayout(
          sidebarPanel(
            tabsetPanel(
              tabPanel("Files",
                div(style = "margin-top: +30px"),
                
                prettySwitch("fastq_align", strong("Use Trimmed Fastq Files?"), bigger = TRUE, fill = TRUE, status = "success", value = TRUE),
                div(style = "margin-top: -15px"),
                helpText(HTML("<font size=-1>
                                Choose whether to use raw fastq files (uploaded in FASTQC) or trimmed fastq files.
                              </font>")),
                
                uiOutput("select_ensembl"),

                div(style = "margin-top: +30px"),

                radioButtons("fasta_how", "FASTA File", choices = c("Upload", "Generate"), inline = TRUE),

                conditionalPanel(condition = "input.fasta_how == 'Upload'",
                                 fileInput("fasta_file", "Choose FASTA File")),

                conditionalPanel(condition = "input.fasta_how == 'Generate'",
                                 div(style = "margin-top: -10px"),
                                 helpText(HTML("<font size=-1>
                                                  1. Click the button to get FASTA from Ensembl.<br>
                                                  Once complete, you will see a dropdown of available chromosomes.<br>
                                                  Note: toplevel - For all except Homo Sapiens & Mouse
                                                </font>")),
                                 
                                 radioButtons("assembly_type", "Assembly Type", choices = c("toplevel", "primary_assembly"), inline = TRUE),
                                 
                                 div(style = "margin-top: +15px"),
                                 
                                 actionBttn("ensembl_fasta", "Get Ensembl Fasta", style = "jelly", icon = icon("dna"), color = "primary", size = "sm"),

                                 div(style = "margin-top: +20px"),

                                 helpText(HTML("<font size=-1>
                                                  2. From the dropdown, select the chromosomes you need. <br>
                                                  Once complete, click the button to get FASTA files for alignment.
                                                </font>")),
                                 uiOutput("select_fasta_chrom"),
                                 actionBttn("load_fasta", "Load Fasta File", style = "jelly", icon = icon("dna"), color = "primary", size = "sm")),

                div(style = "margin-top: +30px"),

                radioButtons("gtf_how", "GTF File", choices = c("Upload", "Generate"), inline = TRUE),

                conditionalPanel(condition = "input.gtf_how == 'Upload'",
                                 fileInput("gtf_file", "Choose GTF File")),

                conditionalPanel(condition = "input.gtf_how == 'Generate'",
                                 div(style = "margin-top: -10px"),
                                 helpText(HTML("<font size=-1>
                                                  1. Click the button to get database from Ensembl (~2-3 min) <br>
                                                  Once complete, you will see a dropdown of available chromosomes.
                                                </font>")),
                                 actionBttn("ensembl_txdb", "Get Ensembl TxDb", style = "jelly", icon = icon("dna"), color = "primary", size = "sm"),

                                 div(style = "margin-top: +20px"),

                                 helpText(HTML("<font size=-1>
                                                  2. From the dropdown, select the chromosomes you need. <br>
                                                  Once complete, click the button to get Txdb for alignment.
                                                </font>")),
                                 uiOutput("select_gtf_chrom"),
                                 actionBttn("load_gtf", "Load TxDb", style = "jelly", icon = icon("dna"), color = "primary", size = "sm"))
              ),

              tabPanel("Settings",
                div(style = "margin-top: +30px"),

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
                actionButton("align", "Align")
              )
            )
          ),

          mainPanel(
            tabsetPanel(
              tabPanel("Ensembl Database",
                dataTableOutput("ensembl_table")
              ),

              tabPanel("Count: Genes",
                       div(style = "margin-top: +20px"),
                       withLoader(dataTableOutput("gene_counts"), type = "html", loader = "loader5")),

              tabPanel("Count: Exons",
                       div(style = "margin-top: +20px"),
                       withLoader(dataTableOutput("exon_counts"), type = "html", loader = "loader5")),

              tabPanel("Count: Junctions",
                       div(style = "margin-top: +20px"),
                       withLoader(dataTableOutput("junction_counts"), type = "html", loader = "loader5"))
            )
          )
        )
      ),
      
      tabPanel("Genome Browser",
        sidebarLayout(
          sidebarPanel(
            tabsetPanel(
              tabPanel("Alignment Stats",
                div(style = "margin-top: +20px"),
                radioButtons("align_stat_type", "Pick Alignment Statistic to Plot",
                             choices = c("Quality Score", "Nucleotide Content by Cycle",
                                         "Duplicate Level", "Mapping Statistics",
                                         "Library Complexity", "Mismatch Frequency",
                                         "Mismatch Types", "Fragment Size"))
              ),

              tabPanel("View Coverage",
                div(style = "margin-top: +20px"),
                uiOutput("select_plot_chrom"),
                uiOutput("select_plot_range"),

                div(style = "margin-top: +30px"),

                numericInput("extend_right", "Extend plot from the right: ", value = 0),
                numericInput("extend_left", "Extend plot from the left:", value = 0),
                actionBttn("view_cov", "View Coverage", style = "unite", size = "sm", color = "success")
              )
            )
          ),

          mainPanel(
            tabsetPanel(
              tabPanel("Alignment Summary",
                       conditionalPanel(condition = "input.align_stat_type == 'Quality Score'",
                                        h3(strong("Quality Score")),
                                        h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                                        HTML("<p><b>
                                              Plot shows the distribution of base quality values as a box plot for each position in the input sequence.
                                              </b></p>")),

                       conditionalPanel(condition = "input.align_stat_type == 'Nucleotide Content by Cycle'",
                                        h3(strong("Nucleotide Content by Cycle")),
                                        h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                                        HTML("<p><b>
                                              Plot shows the frequency of A, C, G, T and N bases by position in the read.
                                              </b></p>")),

                       conditionalPanel(condition = "input.align_stat_type == 'Duplicate Level'",
                                        h3(strong("Duplicate Level")),
                                        h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                                        HTML("<p><b>
                                              Plot shows for each sample the fraction of reads observed at different duplication levels.
                                              </b></p>")),

                       conditionalPanel(condition = "input.align_stat_type == 'Mapping Statistics'",
                                        h3(strong("Mapping Statistics")),
                                        h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                                        HTML("<p><b>
                                              Shows fractions of reads that were (un)mappable to the reference genome.
                                              </b></p>")),

                       conditionalPanel(condition = "input.align_stat_type == 'Library Complexity'",
                                        h3(strong("Library Complexity")),
                                        h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                                        HTML("<p><b>
                                              Shows fractions of unique read(-pair) alignment positions, as a measure of the complexity in the sequencing library.</br>
                                              *Note: This measure is not independent from the total number of reads in a library, and is best compared between </br>
                                              libraries of similar sizes.
                                              </b></p>")),

                       conditionalPanel(condition = "input.align_stat_type == 'Mismatch Frequency'",
                                        h3(strong("Mismatch Frequency")),
                                        h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                                        HTML("<p><b>
                                              Shows the frequency and position (relative to the read sequence) of mismatches in the alignments against the reference genome.
                                              </b></p>")),

                       conditionalPanel(condition = "input.align_stat_type == 'Mismatch Types'",
                                        h3(strong("Mismatch Types")),
                                        h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                                        HTML("<p><b>
                                             Shows the frequency of read bases that caused mismatches in the alignments to the reference genome, separately for each genome base.
                                             </b></p>")),

                       conditionalPanel(condition = "input.align_stat_type == 'Fragment Size'",
                                        h3(strong("Fragment Size")),
                                        h4(span("Description:", style = "color:#73C6B6;font-weight: bold;")),
                                        HTML("<p><b>
                                             Shows the distribution of fragment sizes inferred from aligned read pairs.
                                             </b></p>")),

                       div(style = "margin-top: +20px"),

                       withLoader(plotOutput("align_stat_plot", height = "2000px"), type = "html", loader = "dnaspin")
              ),

              tabPanel("Select Ranges in Current View",
                       div(style = "margin-top: +30px"),
                       awesomeRadio("sel_ranges_type", "Genomic Feature Type", choices = c("gene", "transcript", "exon"), inline = TRUE, checkbox = TRUE),
                       withLoader(DTOutput("selected_ranges"), type = "html", loader = "loader5")),

              tabPanel("View Plot Coverage",
                       div(style = "margin-top: +30px"),
                       awesomeRadio("cov_plot_type", "Select Coverage Plot Type", choices = c("Transcripts", "Genes"), inline = TRUE, checkbox = TRUE),
                       withLoader(plotOutput("cov_plot", height = "1400px"), type = "html", loader = "dnaspin"))
            )
          )
        )
      )
    ),
    
    # Specific Workflows
    navbarMenu("Specific Workflows",
      tabPanel("Alternative Splicing",
        sidebarLayout(
          sidebarPanel(
            downloadButton("download_targets", "Download Targets Template"),
            
            div(style = "margin-top: +20px"),
            
            fileInput("targets_file", "Upload Targets Template with Conditions"),
            
            div(style = "margin-top: +20px"),
            
            numericInput("minReadLength", "Minimum Read Length", value = 100),
            div(style = "margin-top: -10px"),
            helpText(HTML("<font size=-1>
                            Min read length of sequenced library. It is used for computing  E1I and IE2 read summarization.
                          </font>")),
            
            div(style = "margin-top: +20px"),
            
            numericInput("maxISize", "Maximum Intron Size", value = 50000),
            div(style = "margin-top: -10px"),
            helpText(HTML("<font size=-1>Junctions longer than this size will be dicarded</font>")),
            
            div(style = "margin-top: +20px"),
            
            sliderInput("minAnchor", "Min % of Read to be Anchored", min = 0, max = 100, value = 10, post = "%"),
            div(style = "margin-top: -10px"),
            helpText(HTML("<font size=-1>
                            Min % of read that should be aligned to an exon-intron boundary.</br>
                            An intronic junction must overlap completely and at least an minAnchor% into the exon region and the intron region.</br>
                            The regions can be exon1-intron or intron-exon2.
                          </font>")),
            
            div(style = "margin-top: +20px"),
            
            numericInput("threshold", "Min No. of Reads Supporting Junctions", value = 5)
          ),
          
          mainPanel()
        )
      ),
      
      tabPanel("CHIP-Seq")
    )
  )
)