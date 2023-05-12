# CellRanger pipeline for processing FASTQ files
## 0: Converting BAM files back into FASTQ files 
In some cases, only the BAM files can be downloaded from the SRA/ENA (e.g. https://www.ebi.ac.uk/ena/browser/view/PRJEB23051?show=related-records). In these situations, `CellRanger` has a tool that is specifically made to convert BAM files back into FASTQ files for the rest of the `CellRanger` pipeline. Refer to [here](https://support.10xgenomics.com/docs/bamtofastq) for more information. To run the script: 

```
cellranger bam2fastq --nthreads 16 /path/to/bam/file.bam /path/to/output/folder
```

Note that the output folder cannot exist or else this will crash. 

## 1: CellRanger count
Information for this section can be found [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct) on the 10X Genomics site. 

Prior to starting, ensure you have the following files: 
1. The FASTQ files for your dataset. 
2. The reference transcriptome that will be used to align the reads in the FASTQ files to. 
The reference files as well as test data can be found [here](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) on the 10X Genomics site. 

To run the pipeline, enter the following command into the terminal: 
```
cellranger count --id=id_name --fastqs=/path/to/fastq/ --sample=fastq_prefix --transcriptome=/path/to/transcriptome
```
`--id` is the name that the output folder will be given. `--sample` is the prefix (delimited by underscores) the FASTQ files that your sample has. If not specified, all of the files will be included in the count. `--transcriptome` is the path to the reference transcriptome folder that you downloaded and decompressed earlier. 

The output of the command is to the directory you are currently in. 