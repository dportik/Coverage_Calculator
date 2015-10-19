rm(list=ls())

#change path to your 'All_Sample_Coverages.txt' file.

covdata <- read.delim("/Volumes/Portik_Storage/ENDS_COVERAGE/All_coding_beds/All_Sample_Coverages.txt", header = FALSE, sep = "\t")

ls(covdata)

#sample name info
summary(covdata$V1)
#contig name info
summary(covdata$V2)
#coverage
summary(covdata$V3)

######################################
#boxplots by Sample name

boxplot(covdata$V3~covdata$V1, data=covdata, main="Average coverage per contig for all samples", xlab="Sample", 
        ylab="Coverage", col='grey', outline=FALSE, las = 2, cex.axis = 0.5)


######################################
#boxplots by contig name
#****be warned this could be too big to plot meaningfully!!!!!!!!!!!!!

boxplot(covdata$V3~covdata$V2, data=covdata, main="Average coverage per sample for all contigs", xlab="Contig", 
        ylab="Coverage", col='grey', outline=FALSE, las = 2, cex.axis = 0.5)
