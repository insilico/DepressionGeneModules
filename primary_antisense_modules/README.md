# Identification and replication of RNA-Seq gene network modules associated with depression severity.
Trang T. Le, Jonathan Savitz, Hideo Suzuki, Masaya Misaki, T. Kent Teague,  Bill C. White, Julie H. Marino, Graham Wiley, Patrick M. Gaffney, Wayne C. Drevets, Jerzy Bodurka, 
and  Brett A. McKinney

`primary_antisense_modules/` 

To reproduce the gene module analysis, run the R scripts with prefixes 0-4 in order. 

`0.8genes.filtered.corrected.Rdata` contains the filtered and batch corrected antisense log counts per million data matrix, `rnaSeq`. 

`nofilter.cpm.corrected.Rdata` contains the unfiltered and non-logged batch corrected antisense counts per million data matrix, also called `rnaSeq`. 

`enrichmentProfile.Rdata` contains the data matrix of filtered antisense count data collapsed onto 23 module predictors: `collapse.genes`. 

`R01_madrs_data.csv` contains demographic and phenotype data. 
