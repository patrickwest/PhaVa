# PhaVa
PhaVa is an approach for finding potentially **Pha**se **Va**riable invertible regions, also referred to as invertons, in long-read seqeuncing data

## Dependencies
Versions listed are the versions PhaVa has been tested on.
+ python (v3.9+)
+ [EMBOSS](http://emboss.open-bio.org/html/use/ch02s07.html) (v. 6.5.7) [einverted](https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/einverted.html)
+ [minimap2](https://github.com/lh3/minimap2) (v. 2.17)
+ [pysam](https://github.com/pysam-developers/pysam) (v. 0.17.0)
+ [Biopython](https://biopython.org/) (v. 1.81)
+ [mmseqs2](https://github.com/soedinglab/MMseqs2) 

PhaVa is developed and tested on Linux operating systems (CentOS Linux 7), however it should compatible with Mac OSX and Windows

## Usage
The PhaVa workflow is divided into three steps: locate, create, and ratio. 
```
phava locate -i genome.fasta -d out_dir
phava create -d out_dir
phava ratio -r long_reads.fastq -d out_dir
```
Alternatively, all three steps can be run in a single command via variation_wf
```
phava variation_wf -i genome.fasta -r long_reads.fastq -d out_dir
```
Output from each step is centered around a output directory (-d) and should be the same directory for the entire workflow 
The locate and create steps only need to be performed once for a given genome or metagenome, and ratio can then be run on long-read samples using the same output directory (-d)

Any invertons with at least 1 read aligning in the reverse orientation will be found in the output. However, it is strongly recommended to fruther filter based on a minimum reverse read count and minimum % reverse of all reads cutoff (3 and 3% are recommended, respectively)

Expected output:
![Expected output](https://github.com/patrickwest/PhaVa/blob/main/PhavaExpectedOutput-01.png?raw=true)

## Installation
Installing is likely easiest with conda:
```
conda install phava -c bioconda
```
## Testing
PhaVa install can be tested on a small simulated dataset, typically in <1 minute, with:
```
phava test
```
## Tutorial
As an example run through of the PhaVa workflow, we have three isolate long-read sequencing datasets from B. theta, grown in different conditions, we would like to check for invertons. The first step is to identify IRs in the B. theta genome (and find which genes they overlap, if a gene annotation file is provided). In this case, we do want gene overlap information, so we first run prodigal gene prediction on the reference genome
```
prodigal -i GCA_000011065.1_ASM1106v1_genomic.fna -f gff -o GCA_000011065.1_ASM1106v1_genomic.fna.gff
```
Next, we run the phava locate and create steps to identify IRs in the genome, specifying the prodigal output as our gene information
```
phava locate -i GCA_000011065.1_ASM1106v1_genomic.fna -d btheta_inv
phava create -i GCA_000011065.1_ASM1106v1_genomic.fna --genesFormat --genes GCA_000011065.1_ASM1106v1_genomic.fna.gff -d btheta_inv
```
Note that the same output directory (-d) is used for all of our commands in this analysis. A different output directory would only be used if you performed analysis on another reference genome.
Finally, we run the ratio step for each long-read datsaset to map our long-reads and pull out invertons that have reads mapping to the inverted version
```
phava ratio -r condition1.fastq -d btheta_inv
phava ratio -r condition2.fastq -d btheta_inv
phava ratio -r condition3.fastq -d btheta_inv
```
In the output folder, there will be three tsv files containing the putative invertons and the evidence leading to their detection
```
condition1.fastq_vs_GCA_000011065.1_ASM1106v1_genomic.fna
condition2.fastq_vs_GCA_000011065.1_ASM1106v1_genomic.fna
condition3.fastq_vs_GCA_000011065.1_ASM1106v1_genomic.fna
```

Note: instead of running the locate, create, and ratio commands separately, the variation_wf command can be used to perform all three in one step and get the same result
```
phava variation_wf -i GCA_000011065.1_ASM1106v1_genomic.fna --genesFormat --genes GCA_000011065.1_ASM1106v1_genomic.fna.gff -r condition1.fastq -d btheta_inv
phava variation_wf -i GCA_000011065.1_ASM1106v1_genomic.fna --genesFormat --genes GCA_000011065.1_ASM1106v1_genomic.fna.gff -r condition2.fastq -d btheta_inv
phava variation_wf -i GCA_000011065.1_ASM1106v1_genomic.fna --genesFormat --genes GCA_000011065.1_ASM1106v1_genomic.fna.gff -r condition3.fastq -d btheta_inv
```

## Citation
