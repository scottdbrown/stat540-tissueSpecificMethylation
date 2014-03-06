stat540-tissueSpecificMethylation
=================================

Final Proposal: Tissue Specific DNA Methylation and Gene Expression

Previous papers have found a negative correlation between gene expression and DNA methylation. We will start with differential methylation analysis of all genes between 4 neural cell lines and 4 non-neural cell lines, using methylation data from the Illumina 450K Methylation array available through the ENCODE project. The null hypothesis is that there is no tissue specific methylation between neural and other cell types. We will try a few methods of assigning a methylation state to a gene, including looking at probes within the body of the gene, in the gene promoter, and in the upstream CpG Island.

Once we have a list of genes with tissue specific methylation between the ‘other’ and ‘neural’ groups, we will look at the expression levels of those genes using RNA-seq data (two experimental replicates for each cell line, 16 total gene expression files). We will test whether they also show tissue specific expression, and if they do, whether there is a positive or negative correlation between methylation and expression. We expect our hit-list of genes to have neural specific methylation and expression levels, thus should show some overlap with our list of neural-specific marker genes (synapsin 1&2, MAP1B (Microtubule-associated protein 1B), netrin, neuron-specific enolase (NSE), vimentin, etc.).

Ultimately, our hypothesis is that differential methylation between neural cells and other tissues will be correlated with differential gene expression, and the hit list will be especially enriched with neural specific genes. 

With extra time, we propose to validate our hit list with other non-neuronal and neuronal cell lines if available, and add in other epigenetic markers.


Confirmed Cell Lines (No treatment, bold is non-neuronal cell line): 
--------------------------------------
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethyl450/ - Methylation Data
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibRnaSeq/ - HAIB RNA-seq
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/ - CSHL RNA-seq

BE2-C (HAIB-RNASeq + HAIB 450K)

PFSK-1 (HAIB + HAIB)

SK-N-SH (HAIB + HAIB)

U87 (HAIB + HAIB)

**Jurkat (HAIB + HAIB) *T Cell**

**PANC-1 (HAIB + HAIB) *Pancreatic Epithelial**

**GM12878 (CSHL + HAIB) *B Cell**

**HepG2 (CSHL + HAIB) *Liver Epithelial**

Methods to use:
---------------

To perform the methylation analysis, we will take advantage of the annotations for probes from the IlluminaHumanMethylation450K.db package. This will allow us to determine which probes to assign to which gene, and which part of the gene it belongs to. wateRmelon may be used to perform quality control on the methylation data.

For differential methylation analysis, we will use limma. Differential expression analysis using the RNA-seq data will use limma-voom.

Break-down of member responsibilities:
-------------------------
**Jessica:** Differential gene expression analysis (neuronal cells vs non-neuronal cells)  
**Analise:** Data clean-up & exploratory analysis (batch effect btw replicates, remove outliers)  
**Scott:** Differential methylation analysis  
**Howard:** Data normalization  
**Nathaniel:** Differential methylation analysis  