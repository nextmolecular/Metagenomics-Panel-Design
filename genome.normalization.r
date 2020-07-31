
#########################################################################################
## Purpose: To Merge a multi-fasta file of whole genome fasta and contig fasta files from
##          NCBI assembly to create a multi-fasta file with whole genome files only 
##
## Input:   Multi-fasta species file with mixed whole genomes and contigs
##   
## Output:  Multi-fasta species with whole genomes only
##
#########################################################################################


library(Biostrings)
library(tidyverse)

setwd("../Metagenomics-genomes")


######## Accessory Functions####################
################################################

write.fasta = function(an,filename){

  for (i in 1:nrow(an)){
       
       write(paste(">",an$names[i],sep=""),filename,append=TRUE)
       write(paste(an$seq[i]),filename,append=TRUE)

  }

 }



combine.seq = function(seq){

  paste(seq, collapse = '')
}


###############################################
###############################################



##Process an individual set of DNA files for an organism



fastafile = "Bartonella.bacillofromis.all.fasta"
an2 = readDNAStringSet(fastafile)
an3 = str_extract(names(an2),"(NZ|NC)_([:Alpha:]+|[:Digit:]+)")



an4 =  tibble(as.data.frame(an2)) %>%
                  select(seq = x) %>%
                  add_column(names = an3) %>%  
		  group_by(names) %>%
		  summarize(seq = combine.seq(seq)) %>%
		  ungroup() %>%
		  mutate(width = nchar(seq))
    

write.fasta(an4,"anaplasmap.merged.fasta")





##for B.Burgdoferi will use only the primary chromosome, not the plasmids

an5 = an2[which(width(an2)>800000)]
writeXStringSet(an5,"BBurg.all.merged.fasta")


##for B.Microti will use only the primary chromosome, not the plasmids

an5 = an2[which(width(an2)>1000000)]
writeXStringSet(an5,"Babesia.microti.merged.fasta")










