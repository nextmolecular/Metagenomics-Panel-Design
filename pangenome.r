



get.stop.start = function(sequence){
  ##Purpose: Identify start and stop sites in a DNAString
  
  #positive strand
  stop.start = PDict(DNAStringSet(c("atg", "taa", "tag", "tga")))
  res =  matchPDict(stop.start,sequence) 
  start = as.data.frame(res[[1]]) %>% add_column(codon = rep("start",nrow(.)))
  stop1 = as.data.frame(res[[2]]) %>% add_column(codon = rep("stop",nrow(.)))
  stop2 = as.data.frame(res[[3]]) %>% add_column(codon = rep("stop",nrow(.)))
  stop3 = as.data.frame(res[[4]]) %>% add_column(codon = rep("stop",nrow(.)))
  rr1 = tibble(rbind(start,stop1,stop2,stop3)) %>% add_column(strand = rep("1",nrow(.)))
  
  
  #negative strand
  stop.start = PDict(reverseComplement(DNAStringSet(c("atg", "taa", "tag", "tga"))))
  res =  matchPDict(stop.start,sequence) 
  start = as.data.frame(res[[1]]) %>% add_column(codon = rep("start",nrow(.)))
  stop1 = as.data.frame(res[[2]]) %>% add_column(codon = rep("stop",nrow(.)))
  stop2 = as.data.frame(res[[3]]) %>% add_column(codon = rep("stop",nrow(.)))
  stop3 = as.data.frame(res[[4]]) %>% add_column(codon = rep("stop",nrow(.)))
  rr2 = rbind(start,stop1,stop2,stop3) %>% add_column(strand = rep("-1",nrow(.)))
  
  ss = rbind(rr1,rr2)
  ss = ss %>% mutate(frame = start %% 3)
}




get.orfs = function(stop.start){
##Purpose: Generate predicted open-reading-frames from stop.start data

res2.0 = stop.start %>% filter(frame == 0)
res2.1 = stop.start %>% filter(frame == 1)
res2.2 = stop.start %>% filter(frame == 2)


#positive strand

#frame 1
starttemp = -1;
res3 = res2.0 %>% filter(strand ==1) %>% arrange(start)
county = function(start,codon){ if(codon == "start") {if(starttemp == -1){starttemp <<- start} ;return(0)}; if(starttemp == -1)return(0); rval = start-starttemp; starttemp <<- -1; return(rval)}
res0 = res3 %>% rowwise() %>% mutate(x=county(start,codon))
res0 = res0 %>% filter(x > 0) %>% mutate(ostart = start -x) %>% mutate(oend = end) %>% mutate(owidth = oend - ostart + 1) %>%
                                  select(ostart,oend,owidth,strand,frame)
						  

#frame 2
starttemp = -1;
res3 = res2.1 %>% filter(strand ==1) %>% arrange(start)
county = function(start,codon){ if(codon == "start") {if(starttemp == -1){starttemp <<- start} ;return(0)}; if(starttemp == -1)return(0); rval = start-starttemp; starttemp <<- -1; return(rval)}
res1 = res3 %>% rowwise() %>% mutate(x=county(start,codon))
res1 = res1 %>% filter(x > 0) %>% mutate(ostart = start -x) %>% mutate(oend = end) %>% mutate(owidth = oend - ostart + 1) %>%
                                  select(ostart,oend,owidth,strand,frame)

#frame 3
starttemp = -1;
res3 = res2.2 %>% filter(strand ==1) %>% arrange(start)
county = function(start,codon){ if(codon == "start") {if(starttemp == -1){starttemp <<- start} ;return(0)}; if(starttemp == -1)return(0); rval = start-starttemp; starttemp <<- -1; return(rval)}
res2 = res3 %>% rowwise() %>% mutate(x=county(start,codon))
res2 = res2 %>% filter(x > 0) %>% mutate(ostart = start -x) %>% mutate(oend = end) %>% mutate(owidth = oend - ostart + 1) %>%
                                  select(ostart,oend,owidth,strand,frame)

orf.pos = res0 %>% add_row(res1) %>% add_row(res2) %>% na.omit() %>% arrange(desc(owidth))




#negative strand
 
#frame1
starttemp = -1;
res3 = res2.0 %>% filter(strand == -1) %>% arrange(desc(end))
county = function(end,codon){ if(codon == "start") {if(starttemp == -1){starttemp <<- end} ;return(0)}; if(starttemp == -1)return(0); rval = starttemp-end+2; starttemp <<- -1; return(rval)}
res3 = res3 %>% rowwise() %>% mutate(x=county(end,codon))
res4 = res3 %>% filter(x > 0) %>% mutate(ostart = start) %>% mutate(oend = start+x) %>% mutate(owidth = oend -ostart + 1) %>%
                                  select(ostart,oend,owidth,strand,frame)

#frame2
starttemp = -1;
res3 = res2.1 %>% filter(strand == -1) %>% arrange(desc(end))
county = function(end,codon){ if(codon == "start") {if(starttemp == -1){starttemp <<- end} ;return(0)}; if(starttemp == -1)return(0); rval = starttemp-end+2; starttemp <<- -1; return(rval)}
res3 = res3 %>% rowwise() %>% mutate(x=county(end,codon))
res5 = res3 %>% filter(x > 0) %>% mutate(ostart = start) %>% mutate(oend = start+x) %>% mutate(owidth = oend - ostart + 1) %>%
                                  select(ostart,oend,owidth,strand,frame)

#frame3
starttemp = -1;
res3 = res2.2 %>% filter(strand == -1) %>% arrange(desc(end))
county = function(end,codon){ if(codon == "start") {if(starttemp == -1){starttemp <<- end} ;return(0)}; if(starttemp == -1)return(0); rval = starttemp-end+2; starttemp <<- -1; return(rval)}
res3 = res3 %>% rowwise() %>% mutate(x=county(end,codon))
res6 = res3 %>% filter(x > 0) %>% mutate(ostart = start) %>% mutate(oend = start+x) %>% mutate(owidth = oend - ostart + 1) %>%
                                  select(ostart,oend,owidth,strand,frame)




orf.neg = res4 %>% add_row(res5) %>% add_row(res6) %>% na.omit() %>% arrange(desc(owidth))


orf.all = orf.pos %>% add_row(orf.neg)


return(tibble(orf.all))

}






filter.orfs = function(orf){


 orf %>% filter(owidth > 75)


}



get.orf.to.DNAStringSet = function(orf,sequence,sequencename){

bgv = DNAStringSet(Views(sequence,start=orf$ostart,end=orf$oend))
names(bgv) = rep(sequencename,length(bgv))
names(bgv) = paste(names(bgv),1:length(names(bgv)),sep=":")



return(bgv)

}



 

get.core.genome = function(orf.fasta,species.fasta,orf.core.fasta){

#Purpose:Compute the core set of ORFs for a given set of strains in a species
#Input:
#  orf.fasta: name of the fasta file holding the species ORFs 
#  species.fasta: name of the fasta file holding the species genomes
#  orf.core.fasta: name for the fasta file to be output for the core species ORFs
#Output:
#  orf.core.fasta: fasta file of the core species ORFs
#  orf.genome.blast (returned): tibble holding the result of blasting ORFs against species genomes

#command line for making species blastdb
print("Creating Species Database ...")
system(paste("makeblastdb -in ",species.fasta," -dbtype nucl"))


#command line for making orf blastdb
print("Creating ORF Database ...")
system(paste("makeblastdb -in ",orf.fasta," -dbtype nucl"))


#command line for blasting
print("Blasting ORF vs Species Database ...")
system(paste("megablast -d ",species.fasta," -i ",orf.fasta," -D 3 > e1.txt")) 

#parse the megablast result file
#and find core orfs shared by all strains
orf.genome.blast = read_tsv("e1.txt",comment="#",col_names=F)


	   
e2 = orf.genome.blast %>%
     group_by(X1) %>%
	 summarize(target.genomes = n_distinct(X2)) %>%
	 filter(target.genomes ==9) %>%
	 mutate(orf = as.numeric(str_split_fixed(X1,":",n=3)[,3])) %>%
	 mutate(source.organism.a  = str_split_fixed(X1,":",n=3)[,1]) %>%
	 mutate(source.organism.b  = str_split_fixed(X1,":",n=3)[,2]) %>%
	 mutate(source.organism = str_c(source.organism.a,source.organism.b,sep= " ")) %>%
	 arrange(source.organism,orf) %>%
	 select(source.organism,orf,X1) 
	 
	 
fres =readDNAStringSet(orf.fasta)	 
f2 = fres[names(fres) %in% e2$X1]
writeXStringSet(f2,orf.core.fasta)


return(orf.genome.blast)



}





get.orf.clusters = function(orf.genome.blast,centers=2){
#Purpose:
#Input:
#Output:

#find cluster of strains
e3 = orf.genome.blast %>%
     select(X1,X2) %>%
	 distinct(.) %>%
	 mutate(strain = as.vector(str_split_fixed(X1,":[:digit:]+$",n=2)[,1])) %>%
	 select(source.strain = X2, target.strain = strain) %>%
	 group_by(source.strain,target.strain) %>%
	 summarize(n()) %>%
	 select(source.strain=source.strain,target.strain = target.strain,count='n()') %>%
     spread(target.strain,count) %>%
     column_to_rownames("source.strain")	 
	  
e3.kmeans  = kmeans(e3,centers=centers)

return(e3.kmeans$cluster)	  

}	  
	  





get.oligos.from.core.genome = function(fastafile,oligofile){



}








filter.oligos = function(oligofile){



}




find.unique.oligos = function(oligofile,databaselocation){




}



get.primers.from.oligos = function(oligofile,databaselocation){





}




find.unique.primers = function(primerfile,databaselocation){




}







library(tidyverse)
library(Biostrings)

setwd("../Metagenomics-genomes")



####on Babesia Microti genomes


bg = readDNAStringSet("Ehr.fasta")


res2 = get.stop.start(bg[[1]])
res3 = get.orfs(res2) 
res4 = filter.orfs(res3)
res5 = get.orf.to.DNAStringSet(res4,bg[[1]],names(bg[1]))
fres=res5





for (i in 2:length(bg)){

 res2 = get.stop.start(bg[[i]])
 res3 = get.orfs(res2) 
 res4 = filter.orfs(res3)
 res5 = get.orf.to.DNAStringSet(res4,bg[[i]],names(bg[i]))
 

 fres = c(fres,res5)
}




writeXStringSet(fres,"ehrlicia.1.orf.fasta")




res = get.core.genome("ehrlicia.1.orf.fasta","Ehr.fasta","ehrlicia.1.core.orf.fasta")
clusters = get.orf.clusters(res,centers=2)








##############################################################
#  Find common kmers  ########################################
##############################################################



#load species fasta
bg = readDNAStringSet("Ehr.fasta")
f1 = readDNAStringSet("ehrlicia.core.orf.fasta")

f2 = f1[1:1185]



sequence = f2[[224]]


sequencelength = length(sequence)
kmerlength = 150
start1 = 1 
start2 = sequencelength - kmerlength + 1
end1 = kmerlength
end2 = sequencelength



res = DNAStringSet(Views(sequence,start=start1:start2,end=end1:end2))
r2 = PDict(res)


r3 = vcountPDict(r2,bg)
	
###

	 

##first strain orfs	




get.common.oligo = function(sequence,strain.sequences,name,kmerlength=150){
   
   sequencelength = nchar(sequence)
   start1 = 1 
   start2 = sequencelength - kmerlength + 1
   end1 = kmerlength
   end2 = sequencelength

   res = DNAStringSet(Views(sequence,start=start1:start2,end=end1:end2))
   
   r2 = PDict(res)
   r3 = vcountPDict(r2,bg)
   f3 = tibble(as.data.frame(r3))	
   f3[f3 > 0] =1

   
   r4 = tibble(as.data.frame(res)) 
   r5 = r4 %>% mutate(names=rep(name,nrow(r4))) %>% filter(rowSums(f3)>5)
   
   
   return(r5)
   
   
	
}






kmer.collection = f1[1:1185]
f2 = tibble(as.data.frame(kmer.collection)) %>% 
     mutate(name=names(kmer.collection)) %>%
	 rowwise() %>%
	 summarize(res= get.common.oligo(x,bg,name) ) 
	 
	 
	 
	 
	 
	 
	 
f3 = tibble(as.data.frame(f2$res))	
f3[f3 > 0] =1

km2 = kmer.collection[rowSums(f3) > 4]



f4 = f3 %>%
 	 mutate(mm = rowSums(.)) %>%
	 filter(mm > 6)
      	 
	 
	 
	 
 
 
 













sequence = f2[[102]]
	
	











	
	 
	 
	 
	 
 

























 
 
 
	 
	 





































