import sys
import os
import subprocess as sp
import shutil
import numpy as np
import time
'''
python Coverage_Calc_v1.py [directory with bam files] [directory with associated bed files]

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
Use with associated R script: Coverage_Calculator_v1.R

##############
DEPENDENCIES:
numpy - Numerical Python
##############
------------------------
written for Python 2.7.3
Dan Portik
daniel.portik@berkeley.edu
October 2015
------------------------
'''

#go to directory with bam files
bam_dir = sys.argv[1]
os.chdir(bam_dir)

#initiate empty set for sample names
prefix_set = set()

#iterate over files in bam directory
for fl in os.listdir('.'):
    if fl.endswith('_sorted.bam'):
        print "Found {} and adding to list...".format(fl)
        #split the file name
        flnames = fl.split('_')
        #reconstruct the sample prefix
        prefix = flnames[0]+'_'+flnames[1]
        #add prefix to set
        prefix_set.add(prefix)

#convert set into list
prefix_list = list(prefix_set)

#move to bed directory
bed_dir = sys.argv[2]
os.chdir(bed_dir)

#begin sample counter
x = int(1)

#iterate across all sample names
for pref in prefix_list:
    print '\n'
    #for each sample name, scan files in bed directory
    for fl in os.listdir('.'):
        #split file name
        flsplit = fl.split('_')
        match = flsplit[0]+'_'+flsplit[1]
        #and see if matches sample name, if it does:
        if match == pref:
            #create appropriately named SAMPLE_INDEX_coverage.txt output file
            out_name = match+'.coverage.txt'
            #reconstuct path to bam file matching this bed file
            bam_name = bam_dir+'/'+pref+'_sorted.bam'
            #samtools depth -b  bed  sorted.bam > per-site-coverage.txt
            print "Sample {0} - {1}: Now calculating coverage for {2}".format(x,pref,fl)

            #system call to samtools depth function, write to output designated above
            cov_string = "samtools depth -b {0} {1} > {2}".format(fl, bam_name, out_name)
            proc_cov = sp.call(cov_string, shell=True)
    x+=1
    print time.asctime( time.localtime(time.time()) )

#quick function for getting an average of a list of numbers as a string    
def quickstats(x):
    x_array = np.asarray(x,dtype=np.float64)
    x_avg = np.average(x_array)
    x_avg = str(np.around(x_avg, decimals = 1))
    return x_avg

#create master output file
out_file = 'All_Sample_Coverages.txt'
fh_out = open(out_file, 'a')

#for all files in the bed directory, where outputs are writted
for fl in os.listdir('.'):
    #find coverage output files
    if fl.endswith('.coverage.txt'):
        print "Starting calculations at", time.asctime( time.localtime(time.time()) )
        print "Working to calculate averages for all contigs in {}:".format(fl)
        #create empty temp set to fill with contig names
        temp_set = set()
        #open output file
        fh_temp = open(fl,'r')
        #make list of file lines
        lines = fh_temp.readlines()
        #for every line in this file
        for line in lines:
            #clean lines and split by tabs
            line = line.strip()
            line = line.split('\t')
            
            #add the contig name to set
            temp_set.add(line[0])

            #split this line further to reconstruct sample name again
            subs = line[0].split('_')
            sample = subs[0]+'_'+subs[1]

        #turn set into list, go contig by contig across file to calc averages        
        temp_list = list(temp_set)
        temp_list = sorted(temp_list)
        number = len(temp_list)
        print '\t', "There are {} contigs to sift through...".format(number)
        for contig in temp_list:
            print '\t', '\t', "-avg for {}".format(contig)
            #temp list to append coverages across every base of this contig
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
                    newsubs = line[0].split('_')
                    cont_name = newsubs[2]
            #perform averaging for this contig
            contig_avg = quickstats(concovlist)
            #write sample name and average coverage for each contig
            fh_out.write(sample+'\t'+cont_name+'\t'+contig_avg+'\n')
        fh_temp.close()
        print '-----------------------------------------------------------------', '\n'

fh_out.close()

print '\n', "Process complete! That took a long while...", '\n'
