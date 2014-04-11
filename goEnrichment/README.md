Trying to pull in biological meaning by checking Embryonic Stem Cell Specifc genes against the data
========================================================

**This was ultimately ignored for the final project**

I tried to look for Embryonice Stem Cell (ESC) specific genes, and most of this work was done around 2003. The main paper that identified 918 ESC genes no longer had a working link to their data. So, I decided to make a list of ESC genes from what I could find in the literature. I came up with a list of 97 genes that are upregulated in ESC, that can be seen here [`hESC_expressedGenes.txt`](goEnrichment/hESC_expressGenes.txt)

The analysis that I did to find intersecting genes can be found here [`Analise_GOAnalysis.R`](goEnrichment/Analise_GOAnalysis.R)

Furthermore, some heatmaps looking at the expression profiles of the hESC that are found in the differentially expressed/methylated hit lits can be seen here ['hESC_specific_gene_analysis'](plots/hESC_specifc_gene_analysis)

As you can see from these plots there are some of the genes that have a very clear difference between the Stem cells and Somatic cells, and some that are not so convincing.

I also tried to use [DAVID] (http://david.abcc.ncifcrf.gov/)  to look for GO Term enrichment, and with all the different lists of gene intersects and hits from the FDR cut off, I did not find anything very convincing and this part of the analysis was discarded. 



