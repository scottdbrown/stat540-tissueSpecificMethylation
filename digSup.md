An Exploration of DNA Methylation in Stem vs. Somatic Cells, and How It Relates to Gene Expression
====================================
STAT 540 - Digital Supplement
====================================

Poster
---------
A pdf of the [printed poster](poster/gsat540_v2.pdf) is available.

Data Acquisition
---------------------
We obtained data from the [Nazor et al. Cell Stem Cell 2012](http://www.ncbi.nlm.nih.gov/pubmed/22560082)^^1 study.  

All data was availble on [GEO](http://www.ncbi.nlm.nih.gov/geo/). 450K Methylation Array data was available under [`GSE31848`](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31848), and was acquired via the [`getMethylation.R`](dataAcquisition/getMethylation.R) script. HT12v3 Expression Array data was available under [`GSE30652`](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30652), and was acquired via the [`getExpression.R`](dataAcquisition/getExpression.R) script.

Data Analysis and Output
----------------------------

### Methylation
We chose two methods of assigning methylation status to a gene. After obtaining gene annotations for each probe from the `IlluminaHumanMethylation450k.db` R package, and removing probes not annotated to a gene, we averaged all probes in a sample that had an annotation to each gene. We also performed the same averaging, but for probes which were annotated to the promoter region of a gene only. The code for this analysis is in [`assignMethylationToGene.R`](methylation/assignMethylationToGene.R).

To normalize the methylation values...**TO BE FILLED IN**

To perform differential methylation...**TO BE FILLED IN**

### Expression
As the expression array had more than 1 probe for a given gene in some cases, we averaged the value from all probes for a given gene, much the same way we did for methylation values. This analysis was taken care of in [`RNAPreProc.R`](expression/RNAPreProc.R).

To normalize the expression values...**TO BE FILLED IN**

To perform differential expression...**TO BE FILLED IN**

### Correlation
The [`exploreCorrelation.R`](correlation/exploreCorrelation.R) script takes the TopTable results from the differential methylation and differential expression analysis as inputs, as well as the normalized data for both. Pearson correlations are calculated for genes present on each platform (the intersect), using the log2 transformed expression values and normalized M-values for methylation.

A random sample was chosen to illustrate the correlation between methylation M values and log2 expression values, seen [here](plots/example_correlation.pdf). 

The main figure showing the differences in correlations between tissue type, gene group, and methylation metric is [here](plots/correlations_by_cell_group_CpG2.pdf). Pair-wise t-tests were performed for all meaningful pairs (keeping two variables constant while changing the third), and the results are summarized in [ttest_table1.xlsx](correlation/ttest_table1.xlsx). Grey boxes are the groups being compared, while coloured boxes indicate the status of the other variables.

Venn diagrams were made to show [A)](plots/gene_overlap_venn.pdf) the overlap between genes covered by the different platforms and metrics, and [B)](plots/differential_gene_overlap_venn.pdf) the overlap between differentially expressed and methylated genes on the different platforms and metrics. Note that differential analysis was done on the total set of genes available on that platform, and not only on the intersecting genes between platforms.


### GO Analysis

We used GOrilla^^2 ...takes a ranked gene list...our lists were ranked by differnetial p-value...intersects used...rough number of genes...**TO BE FILLED IN**

We used [GOrilla] (http://www.biomedcentral.com/1471-2105/10/48)^^2, a web-based GO-enrichment tool on our analysis results to evaluate for biological relevance. GOrilla requires only ranked gene lists and has no known gene number limits, which was both suitable for our use and generated results at a very fast rate. (Comparatively, [DAVID] (http://david.abcc.ncifcrf.gov/) only accepts at most 3000 HUGO/Official gene symbols, which is unsuitable for our work, as some of our lists exceed 5000 genes).

**TBC**


Future Work
-------------
TO BE FLESHED OUT MORE
- Regarding methylation, looking at CpG Islands, Gene Body separate from promoter, etc.
- binarizing expression and methylation prior to looking at correlation (up/downregulated expression and hyper/hypomethylation)


References
-------------

1. Nazor, Kristopher L., et al. "Recurrent variations in DNA methylation in human pluripotent stem cells and their differentiated derivatives." Cell stem cell 10.5 (2012): 620-634. [Link](http://www.ncbi.nlm.nih.gov/pubmed/22560082)
2. Eden, Eran, et al. "GOrilla: a tool for discovery and visualization of enriched GO terms in ranked gene lists." BMC bioinformatics 10.1 (2009): 48. [Link](http://www.biomedcentral.com/1471-2105/10/48)