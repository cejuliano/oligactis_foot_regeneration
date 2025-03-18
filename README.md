# oligactis_foot_regeneration
Transcriptomic timecourse analysis of foot regeneration in Hydra oligactis

# RNA-seq data analysis for: Campos *et al.* 2025, Wnt signaling restores evolutionary loss of regenerative potential in Hydra

This repository contains the scripts and files required to analyze a time course of regeneration in *Hydra oligactis* from RNA-seq data as presented in our manuscript. The repository is divided into two sections based on the type of computer environment in which the scripts in each section can be executed. The `cluster` folder contains scripts meant to be executed on a computing cluster running Ubuntu that uses the slurm workload manager. The `local` contains scripts meant to be executed on macOS Sonoma (14.4.1) running R (4.3.3) on RStudio (2024.04.1). `Cluster` scripts perform raw read trimming, mapping and expression estimation. `Local` scripts calculate differential gene expression analysis, coexpression networks and differential gene expression patterns.

### Accessing required files not available on the Github repository

Due to file size limitations on Github, this repository does not contain RNA-seq raw reads which can be found at the NCBI (https://www.ncbi.nlm.nih.gov) accession: PRJNA1231128.

Descriptions for the files contained in this repository can be found in the `fileDescriptions.txt` file.

#### Publicly available resources

[RSEM](https://github.com/deweylab/RSEM) - required in the `cluster/resources` folder.

[trimmomatic v 0.36](http://www.usadellab.org/cms/?page=trimmomatic) - required in 
the `cluster/resources` folder.

[Java] (https://www.java.com/en/download/help/mac_install.html) - required in the `cluster/resources` folder 

```
sudo apt install default-jre
```

## Cluster Analysis

#### Setup

The first step is downloading the `cluster` folder from this repository. Raw RNA-seq reads from this study need to be downloaded to the `cluster` directory as gzipped fastq files. Raw data for this study can be found at NCBI with the accession [ ] .

The analysis requires the following software pacakges [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.html), [rsem](https://deweylab.github.io/RSEM/), [samtools](http://www.htslib.org/), [trimmomatic v 0.36](http://www.usadellab.org/cms/?page=trimmomatic), [jdk](https://www.oracle.com/java/technologies/downloads/) and [Python3](https://www.python.org/download/releases/3.0/) installed and added to the user's path.

#### Installing python packages

The `cluster` analysis requires a number of python packages (listed in the `cluster/resources/python_requirements.txt` file). To create a virtual environment in which to download the required packages and run the analysis execute the following commands:

```         
python3 -m venv resources/venv
source resources/venv/bin/activate
pip -r resources/python_requirements.txt
deactivate
```

#### Indexing reference sequences

Read mapping requires indexing the provided fasta files using either bowtie2 or rsem. If your computing cluster uses slurm, this can be done by executing the following command:

```         
sbatch -p [partition name] resources/slurm_prepare_references.sh
```

On other linux systems, you can execute:

```         
cd resources
./prep_references.sh
```

#### RNA-seq data processing

The RNA read mapping pipeline first removes index sequences and low quality basecalls using trimmomatic, and then maps reads and quantifies transcript read counts using rsem. It will also generate fastqc reports from both the raw and trimommatic-filtered data.

If your computing cluster uses slurm, you can execute:

```         
sbatch -p [partition name] resources/slurm_RNA_processing.sh
```

For a system without slurm, you can run the pipeline on each individual sample using the `RNA_Mapping_Pipeline.sh` script like so:

```         
./resources/RNA_processing.sh [insert prefix here]
```

After the pipeline has finished running for all samples, read count matrices need to be calculated for the data generated in this study and for the data from Wenger et al. you can execute the code within your terminal session by running:

```         
cd resources
rsem-generate-data-matrix *RNA.genes.results > OLI.mat.txt
```

## Local Analysis

#### Setup

The following section describes the `local` analysis, which was designed to be run on a personal computer running MacOS.

For the RNA-seq analyses, we have located `OLI_mat.txt` file produced by the code in the `cluster` folder to the `local/resources`folder for its analysis.

#### Differential Gene Expression

We use edgeR to identify differentially expressed genes between compared time points using a negative binomial generalized linear model over the fold change values (glmTreat). These scripts output the results both as a RData object and as individual results tables in a csv format for *H. oligactis* and *H. vulgaris*.

```         
Rscript DGE_analysis_OLI.R
Rscript DGE_analysis_AEP.R
```

#### OrthoClust Analysis

Normalized counts for *H. oligactis* and *H. vulgaris* can then be used to analyze co-expression of orthologous genes using OrthoClust. To do this, we first determine orthology by reciprocal blasting the protein sequences from H. oligactis and H. vulgaris in fasta format which is done using a MacOS terminal.

Translated amino acid sequences for *H. oligactis* reference transcriptome were obtained using from [https://hydratlas.unige.ch/](https://hydratlas.unige.ch/blast/blast_link.cgi){.uri}, and the counterpart for *H. vulgaris* transcriptomic reference used in Cazet, *et al.*, 2021, was downloaded from the *Hydra* 2.0 Genome Project Portal (<https://research.nhgri.nih.gov/hydra/>). The last step of reciprocal blasting requires the fasta sequences of both species to be concatenated into the same file, which was done using MacOS terminal

```         
cat ./resources/OLI.AA.fasta ./resources/AEP.AA.fasta > ./resources/OLI-AEP.AA.fasta
```

NCBI tools were downloaded from [https://www.ncbi.nlm.nih.gov/](https://www.ncbi.nlm.nih.gov/books/NBK569861/){.uri} to be used on MacOS terminal. A Blastp is performed prior to the identification of orthologous genes in the `resources` folder.

```         
### Build BLAST database  ===============================================================
cd processed_data
makeblastdb -in OLI-AEP.AA.fasta \
            -out OLI-AEP.AA.blastdb \
            -dbtype prot \
            -logfile makeblastdb.log

### Perform BLAST analysis  ===============================================================
# this step is performed in "processed_data" folder
blastp -evalue 0.00001 \
  -outfmt 6 -db OLI-AEP.AA.blastdb \
  -query OLI-AEP.AA.fasta > OLI-AEP.blastout 
```

The result from the Blast above can be now used with the `ReciprocalBlastHit.py` script in the `local/resources` folder on a MacOS terminal with python installation. The syntax to use this script requires to name the initials for the ID names of each species:

```         
python ReciprocalBlastHit.py OLI-AEP.blastout R Sc OLI-AEP.RBH.txt
```

We have included the `OLI-AEP.blastout` and `OLI-AEP.RBH.txt` files in the `local/resources` folder to be used directly on the next part of the OrthoClust analysis pipeline.

Before running the OrthoClust analysis the OrthoClust function requires to be installed in R. We have provided the "tarball" file `OrthoClust_1.0.tar.gz` to do this in the `local/resources` folder:

```         
install.packages("resources/OrthoClust_1.0.tar.gz", repos = NULL, type = "source")
```

The OrthClust analysis requires expression data in the form of fpkm for both *H. oligactis* and *H. vulgaris.* The fpkm file for *H. oligactis* was obtained by modifying RSEM function `rsem-generate-data-matrix` by chaning "my \$offsite = 4;" to "my \$offsite = 5;" within the function. For *H. vulgaris* the fpkm was calculated with the provided Rscript `FPKM_calcuculation.R` using the normalized counts from DGE analysis and the `dovetail_length.csv` file containing the length associated to each *H. vulgaris* transcript ID. This script outputs the `fpkm.AEP.csv` file provided in the `local/resources` folder.

```         
Rscript FPKM_calculation.R
```

OrthoClust analysis can be run directly with our modified Rscript with the files provided in the `local/resources` folder. This script outputs a summary list containing cluster IDs and number of genes and gene ratio for each species as a .csv. In addition, it generates a table containing association of each transcript ID from both species to a specific cluster ID.

```         
Rscript OrthoClust_analysis.R
```

#### maSigPro Analysis

To determine differential gene expression patterns between foot regenerating tissue wit alsterpaullone treatment and foot regeneration controls with DMSO, we used maSigPro pipeline in R. Running this script requires `mat.ALP.OLI.csv` a file containing counts per million of alsterpaullone and DMSO treated samples and `design.ALP.OLI.csv` a file with details about expression time, replicates, and treatments. Both of this files have been provided in the `local/resources` folder. This analysis can be run in either Rstudio or by running `DGE_pattern_anlysis.R` in terminal.

```         
Rscript DGE_pattern_analysis.R
```

