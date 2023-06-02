# PhaVa
PhaVa is an approach for finding potentially **Pha**se **Va**riable invertible regions, also referred to as invertons, in long-read seqeuncing data

## Dependencies
+ [EMBOSS](http://emboss.open-bio.org/html/use/ch02s07.html) [einverted](https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/einverted.html)
+ [minimap2](https://github.com/lh3/minimap2)
+ [pysam](https://github.com/pysam-developers/pysam)
+ [Biopython](https://biopython.org/)

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

## Installation
```
git clone https://github.com/patrickwest/PhaVa
```
## Citation
