
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
<div class="row">
   <img src="https://github.com/paigerollex/gene_cloud_omics/blob/main/screenshots/1.png" alt="Fig 1. Upload fastq" style="width:100%">
   <figcaption>Fig 1. Upload fastq.</figcaption>
</div>

