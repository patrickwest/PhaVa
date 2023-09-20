# PhaVa
PhaVa is an approach for finding potentially **Pha**se **Va**riable invertible regions, also referred to as invertons, in long-read seqeuncing data

## Dependencies
Versions listed are the versions PhaVa has been tested on.
+ [EMBOSS](http://emboss.open-bio.org/html/use/ch02s07.html) (v. 6.5.7) [einverted](https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/einverted.html)
+ [minimap2](https://github.com/lh3/minimap2) (v. 2.17)
+ [pysam](https://github.com/pysam-developers/pysam) (v. 0.17.0)
+ [Biopython](https://biopython.org/) (v. 1.81)

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

## Installation
Beyond installing dependencies, PhaVa install is:
```
git clone https://github.com/patrickwest/PhaVa
```
## Testing
PhaVa install can be tested on a small simulated dataset with pytest and a pytest module located in the 'tests' subdirectory:
```
pytest phava_test.py
```
## Citation
