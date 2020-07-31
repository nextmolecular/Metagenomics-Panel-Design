



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
  stop.start = PDict(DNAStringSet(c("atg", "taa", "tag", "tga")))
  res =  matchPDict(stop.start,reverseComplement(sequence)) 
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
res3 = res2.0 %>% filter(strand == -1) %>% arrange(start)
county = function(start,codon){ if(codon == "start") {if(starttemp == -1){starttemp <<- start} ;return(0)}; if(starttemp == -1)return(0); rval = start-starttemp; starttemp <<- -1; return(rval)}
res3 = res3 %>% rowwise() %>% mutate(x=county(start,codon))
res4 = res3 %>% filter(x > 0) %>% mutate(ostart = start -x) %>% mutate(oend = end) %>% mutate(owidth = oend - ostart + 1) %>%
                                  select(ostart,oend,owidth,strand,frame)

#frame2
starttemp = -1;
res3 = res2.1 %>% filter(strand == -1) %>% arrange(start)
county = function(start,codon){ if(codon == "start") {if(starttemp == -1){starttemp <<- start} ;return(0)}; if(starttemp == -1)return(0); rval = start-starttemp; starttemp <<- -1; return(rval)}
res3 = res3 %>% rowwise() %>% mutate(x=county(start,codon))
res5 = res3 %>% filter(x > 0) %>% mutate(ostart = start -x) %>% mutate(oend = end) %>% mutate(owidth = oend - ostart + 1) %>%
                                  select(ostart,oend,owidth,strand,frame)

#frame3
starttemp = -1;
res3 = res2.2 %>% filter(strand == -1) %>% arrange(start)
county = function(start,codon){ if(codon == "start") {if(starttemp == -1){starttemp <<- start} ;return(0)}; if(starttemp == -1)return(0); rval = start-starttemp; starttemp <<- -1; return(rval)}
res3 = res3 %>% rowwise() %>% mutate(x=county(start,codon))
res6 = res3 %>% filter(x > 0) %>% mutate(ostart = start -x) %>% mutate(oend = end) %>% mutate(owidth = oend - ostart + 1) %>%
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

system(paste("makeblastdb -in",fastfile)




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


bg = readDNAStringSet("Babesia.microti.merged.fasta")
res2 = get.stop.start(bg[[1]])
res3 = get.orfs(res2) 
res4 = filter.orfs(res3)


res5 = get.orf.to.DNAStringSet(res4,bg[[1]])


























