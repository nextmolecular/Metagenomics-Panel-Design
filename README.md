# Metagenomics-Panel-Design
Improved workflow design for designing the Metagenomics panel

* Genome annotation of each species using protka
  
  Example:
  prokka --kingdom Bacteria --outdir prokka_GCA_000008285 --genus Listeria --locustag GCA_000008285 GCA_000008285.1_ASM828v1_genomic.fna
 
* Derivation of species-specific regions using Roary on genomes of known strains

  Example:
  roary -f ./demo -e -n -v ./gff/*.gff

* Selection of regions for conserved protein-coding sequences and not rRNA

* Decomposition of species-specific regions into 200-mers

* Screening of 200-mers for 40-60% GC content and for minimizing runs of homopolymers

* Screening of 200-mers by lack of homology to hg38 and the Bacterial Pan-genome

* Selection of primers to generate 100 Amplicons/species-target that are uniformly
  arrayed accross the genome of the organism.

* Screening of primers from 200-mers by lack of homology to h38 and the Bacterial Pan-Genome

