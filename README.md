# PhaVa
PhaVa is an approach for finding potentially **Pha**se **Va**riable invertible regions, also referred to as invertons, in long-read seqeuncing data

## Dependencies
+ [einverted](https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/einverted.html)
+ [minimap2](https://github.com/lh3/minimap2)
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
Output from each step is centered around a output directory (-d) and should be the same directory for each command.

## Installation
