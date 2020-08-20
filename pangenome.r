



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
county = function(start,codon){ if(codon == "start") {if(starttemp == -1){starttemp <<- end} ;return(0)}; if(starttemp == -1)return(0); rval = starttemp-end+2; starttemp <<- -1; return(rval)}
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



get.orf.to.DNAStringSet = function(orf,sequence){

bgv = DNAStringSet(Views(sequence,start=orf$ostart,end=orf$oend))
return(bgv)

}






make.species.database = function(fastafile){

system(paste("makeblastdb -in ",fastfile," -dbtype nucl" ))




}


get.core.genome = function(fastafile,species.database){





}





get.unique.core.genome = function(fastafile,databaselocation){


}





get.oligos.from.unique.core.genome = function(fastafile,oligofile){



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
res5 = get.orf.to.DNAStringSet(res4,bg[[1]])
names(res5) = rep(names(bg[1]),length(res5))
names(res5) = paste(names(res5),1:length(names(res5)),sep=":")

fres = res5





for (i in 2:length(bg)){

 res2 = get.stop.start(bg[[i]])
 res3 = get.orfs(res2) 
 res4 = filter.orfs(res3)
 res5 = get.orf.to.DNAStringSet(res4,bg[[i]])
 names(res5) = rep(names(bg[i]),length(res5))
 names(res5) = paste(names(res5),1:length(names(res5)),sep=":")


 fres = c(fres,res5)
}




writeXStringSet(fres,"ehrlicia.1.orf.fasta")


 


#command line for making blastdb
system("makeblastdb -in ehrlicia.1.orf.fasta -dbtype nucl")


#command line for blasting
system("megablast -d  Ehr.fasta -i ehrlicia.1.orf.fasta -D 3 > e1.txt")


#parse the megablast result file




e1 = read_tsv("e1.txt",comment="#",col_names=F)

      
	   
	   
	   
e2 = e1 %>%
     group_by(X1) %>%
	 summarize(target.genomes = n_distinct(X2)) %>%
	 filter(target.genomes ==9) %>%
	 mutate(orf = as.numeric(str_split_fixed(X1,":",n=3)[,3])) %>%
	 mutate(source.organism.a  = str_split_fixed(X1,":",n=3)[,1]) %>%
	 mutate(source.organism.b  = str_split_fixed(X1,":",n=3)[,2]) %>%
	 mutate(source.organism = str_c(source.organism.a,source.organism.b,sep= " ")) %>%
	 select(source.organism,orf,X1) %>%
	 arrange(source.organism,orf) 
	 
	 
	 
	 e2 %>% group_by(source.organism) %>% summarize(orfnum = n())
	 
	 
	 
	 
	 
     
	 
	 





































