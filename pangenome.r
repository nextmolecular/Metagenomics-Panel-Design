
library(tidyverse)
library(Biostrings)




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
}






get.orfs = function(start.stop){
##Purpose: In-silico identification of orfs from start and stop sites

       orfs = function(sstart,end,width,codon,stops){
          orf =   stops %>%
                  filter(end > sstart) %>%
	  	          mutate(start = sstart,width=end-start+1,strand=1) %>%
		          select(start,end,width,strand)

          return (orf)
       }
	   
       rvorfs = function(sstart,end,width,codon,stops){
          orf =   stops %>%
                  filter(end < sstart) %>%
	  	          mutate(start = sstart,width=start-end+1,strand=-1) %>%
		          select(start,end,width,strand)

          return (orf)
       }

      starts.f = start.stop %>% filter(codon == "start") %>%  filter(strand == 1)
      stops.f = start.stop %>% filter(codon == "stop") %>% filter(strand == 1)
      starts.r = start.stop %>% filter(codon == "start") %>%  filter(strand == -1)
      stops.r = start.stop %>% filter(codon == "stop") %>% filter(strand == -1)

      res =  starts.f %>%
            rowwise() %>%
	    summarize(orfs(start,end,width,codon,stops.f)) 

      res2 =  starts.r %>%
            rowwise() %>%
	    summarize(rvorfs(start,end,width,codon,stops.r)) 

      #combine postive and negative strand orfs and filter by integer triplets
      ret = rbind(res,res2)  %>% filter(width %% 3 == 0)
   
}




seq = DNAString("atgaaaatgcagtaacccatgccc")
ss = get.stop.start(seq)
orfs = get.orfs(ss)




####on reverse strand



res2 = get.stop.start(reverseComplement(seq))
orfs2 = get.orfs(res2,reverse=TRUE)


























