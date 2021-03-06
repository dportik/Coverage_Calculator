# Coverage_Calculator

python Coverage_Calc.py [directory with bam files] [directory with associated bed files]

The 'Sample_indexXX__sorted.bam' files are located in the 'intarget_assemblies' folder resulting 
from script 6 (6-TransExonCapPhylo 'contig' option) of the 
https://github.com/MVZSEQ/denovoTargetCapturePhylogenomics workflow.

The 'Sample_indexXX_coding.bed' style files are found in the 'assemblies/In_target' folder resulting
from script 5 (5-FindingTargetsPhylo).

Naming of files should follow the general format:
'JMPD004_index20_coding.bed'
or
'JMPD004_index20_flanking_ONLY.bed'
(or you can use flanking and coding bed)
and
'JMPD004_index20_sorted.bam'

The "samplename_indexNo_" naming scheme for files is important for all components of this script.

COPY THE BED FILES OF INTEREST to a new directory before you start!  Else this will run over
every bed file in the directory per sample, which is defaulted as FIVE!

-Performs coverage calculation across entire bed file for each sample

-Writes output file per sample

-Iterates over all output files, calculates an average coverage per CONTIG (not EXON per CONTIG)

-Appends averages for individual contigs across all sample to 'All_Sample_Coverages.txt', format:

JMPD002_index10	Contig809	31.8

JMPD002_index10	Contig81	254.1

JMPD002_index10	Contig810	23.7

JMPD002_index10	Contig811	67.7

Can be used to create meaningful boxplots in R, where factor is sample or contig name.

written for Python 2.7.3

DEPENDENCIES:

numpy - Numerical Python

samtools - needs to be in path to call from command line

------------------------

**Citation Information:**

These scripts were published as part of the following work:

+ *Portik, D.M., Smith, L.L., and K. Bi. 2016. An evaluation of transcriptome-based exon capture for frog phylogenomics across multiple scales of divergence (Class: Amphibia, Order: Anura). Molecular Ecology Resources 16: 1069–1083. https://doi.org/10.1111/1755-0998.12541*

If you use or modify these scripts for your own purposes, please cite our publication.


**Contact:**

Daniel Portik, PhD

Postdoctoral Researcher

University of Arizona

daniel.portik@gmail.com
