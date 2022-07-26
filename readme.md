
## Key Features
1. **FastQC**: Provides quality control for fastq files resulting from high throughput sequencing. Resulting plots can help to identify potential errors during sequencing and guide trimming configuration.
   
2. **ngsReports**: R package to parse FastQC reports and display aggregated plots for multiple fastq files.
   
3. **Rfastp**: R wrapper of fastp which supports quality and length filtering, adapter detection, window and global trimming, overrepresentation analysis as well as UMI sequence/id handling of fastq files.

4. **Rhisat2**: R wrapper of HISAT2 for splice-aware alignment
   
5. **Rbowtie**: R wrapper of Bowtie for fasta and memory-efficient short read alignment (unmapped alignment, not splice-aware) as well as SpliceMap for spliced alignments.
   
6. **QuasR**: R package which provides a direct pipeline from aligning with Rhisat2 and Rbowtie to generating feature counts, alignment quality reports and wig files for coverage plots.
   
7. **ramwas**: R package mainly used for methylome-wide association studies, which helps to plot alignment statistics from bam files generated after alignment.
   
8. **Gviz**: R graphics package to display a genome browser showing coverage as well as gene and transcript annotation within selected range and chromosome.


## Steps
1. <a href="https://github.com/paigerollex/gene_cloud_omics/blob/main/guide/Uploading%20Fastq.pdf">Upload fastq files</a>
   
2. <a href="https://github.com/paigerollex/gene_cloud_omics/blob/main/guide/FastQC.pdf">Run FastQC</a>
   
3. <a href="https://github.com/paigerollex/gene_cloud_omics/blob/main/guide/Trim%20with%20Rfastp.pdf">Trim fastq with Rfastp</a>
   
4. <a href="https://github.com/paigerollex/gene_cloud_omics/blob/main/guide/Aligning.pdf">Upload/Download annotation files and align</a>
   
5. <a href="https://github.com/paigerollex/gene_cloud_omics/blob/main/guide/Alignment%20Summary.pdf">Alignment statistics and Genome Browser</a>


## Interface
### Step 1: Upload fastq
<p float="left">
  <img src="https://github.com/paigerollex/gene_cloud_omics/blob/main/screenshots/1.png" width="100%" />
</p>

### Step 2: FastQC
<p>
  <img src="https://github.com/paigerollex/gene_cloud_omics/blob/main/screenshots/2a.png" width="49%" />
  <img src="https://github.com/paigerollex/gene_cloud_omics/blob/main/screenshots/2b.png" width="49%" /> 
</p>

### Step 3: Trimming
<p>
  <img src="https://github.com/paigerollex/gene_cloud_omics/blob/main/screenshots/3a.png" width="49%" />
  <img src="https://github.com/paigerollex/gene_cloud_omics/blob/main/screenshots/3b.png" width="49%" /> 
</p>
<p>
  <img src="https://github.com/paigerollex/gene_cloud_omics/blob/main/screenshots/3c.png" width="49%" />
  <img src="https://github.com/paigerollex/gene_cloud_omics/blob/main/screenshots/3d.png" width="49%" /> 
</p>

### Step 4: Align
<p>
  <img src="https://github.com/paigerollex/gene_cloud_omics/blob/main/screenshots/4a.png" width="49%" />
  <img src="https://github.com/paigerollex/gene_cloud_omics/blob/main/screenshots/4b.png" width="49%" /> 
</p>

### Step 5: Plot Alignment
<p align="centre">
  <img src="https://github.com/paigerollex/gene_cloud_omics/blob/main/screenshots/5a.png" width="49%" />
</p>

<p>
  <img src="https://github.com/paigerollex/gene_cloud_omics/blob/main/screenshots/5b.png" width="49%" />
  <img src="https://github.com/paigerollex/gene_cloud_omics/blob/main/screenshots/5c.png" width="49%" /> 
</p>


## References
* Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at:
  <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>
* Ward C, To T, Pederson S (2019). “ngsReports: An R Package for managing FastQC reports
  and other NGS related log files.” _Bioinformatics_. doi:10.1093/bioinformatics/btz937
  <https://doi.org/10.1093/bioinformatics/btz937>
* Wang W, Carroll T (2022). _Rfastp: An Ultra-Fast and All-in-One Fastq Preprocessor
  (Quality Control, Adapter, low quality and polyX trimming) and UMI Sequence Parsing)._. R
  package version 1.6.0.
* Soneson C (2022). _Rhisat2: R Wrapper for HISAT2 Aligner_. R package version 1.12.0,
  <https://github.com/fmicompbio/Rhisat2>.
* Langmead B, Trapnell C, Pop M, Salzberg SL. Ultrafast and memory-efficient alignment of
  short DNA sequences to the human genome. Genome Biology 10(3):R25 (2009).
* Au KF, Jiang H, Lin L, Xing Y, Wong WH. Detection of splice junctions from paired-end
  RNA-seq data by SpliceMap. Nucleic Acids Research, 38(14):4570-8 (2010).
* Hahne F, Lerch A, Stadler MB. Rbowtie: An R wrapper for bowtie and SpliceMap short read
  aligners. (unpublished)
* Gaidatzis D, Lerch A, Hahne F, Stadler MB. QuasR: Quantification and annotation of short
  reads in R. Bioinformatics 31(7):1130-1132 (2015).
* Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory
  requirements. Nat Methods, 12(4):357-60 (2015).
* Andrey A. Shabalin, Mohammad W. Hattab, Shaunna L. Clark, Robin F. Chan, Gaurav Kumar,
  Karolina A. Aberg, and Edwin J.C.G. van den Oord. RaMWAS: fast methylome-wide association
  study pipeline for enrichment platforms, Bioinformatics 2018,
  <http://dx.doi.org/10.1093/bioinformatics/bty069>
* Hahne F, Ivanek R. Visualizing Genomic Data Using Gviz and Bioconductor. Methods Mol
  Biol. 1418:335-51 (2016).
  
  
