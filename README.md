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

Meeting Minutes (Meeting with Paul)
====================================
-xyplot, with rank correlation? linear correlation? look at every pair. Pearson correlation? point for every gene… (gene + methylation pair) 
   null: no relationship between gene + methylation pair, cor.test() will give pair p-values, can claim that the overall correlation across all genes between methylation and gene expression is negative. How to make sure it is not an artifact? -- randomly arrange gene with methylation of diff gene, if hypothesis is correct the made up pairs should be zero. (want a different result) 
- want to be able to explain one thing by the other, and the better the explanation is the parameter is telling me something
-artificial grouping (need to defend why..) 
-work in the GO terms of why some genes follow the hypothesis, are they on certain chromosomes, do they have certain functional groups 
-if we use another tool we need to know how the statistical method works
-if there are 100 genes (which ones have a negative vs. positive correlation might be interesting, look at which tissues the pattern is more pronounced in) 
-cluster the genes, look at contrasts in differential analysis, and look at different groupings
-residual, what is left when you take out what is accounted for by the model, and the residual is what is left after the math is done, error is what is actually going on in the data. 
-lmFit object you can call residuals on it 

Look for a larger data set, with match expression and methylation, maybe tissue analysis.  

-Look for more data to combine

-Brain-span data set, take some of it

-NIH epigenetics road map, may have this data as well 

-10 tissues 
- 450K data want condensed to one gene, taking the average would be ok for this project

New Data
=========
We will be using data from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3348513/ for data on somatic tissues. 

For discussion on March 20th 2014: 
-----------------------------------
Supplementary Table 1 of this paper has fairly detailed information about the samples.

Of the 336 samples (at the patient/cell line level), 228 had both expression and methylation analysis performed. 87 of those had the methylation analysis on 450K (the rest were the older 27K array).

Of these 87 samples, there are 26 somatic tissue (taken from a patient, see table below of counts from each tissue), 20 somatic primary (cell lines, see counts below), 21 iPS cells (different passages), and 30 hESCs.

Somatic Tissues & Counts:
Adipose  2
Adrenal	4
Bladder	2
Kidney	2
Liver	1
Lung	4
Lymph.node	2
Skeletal.muscle	2
Stomach	4
Tongue	1
Ureter	2

Somatic Primary & Counts:
Fetal Lung Fibroblast	1
Fibroblast, lung	1
Human Foreskin Fibroblast	1
Human_Bladder_Smooth_Muscle_Cells	1
Human_Brain_Vascular_Smooth_Muscle_Cells	1
Human_Cardiac_Fibroblasts	1
Human_Cardiac_Fibroblasts_Adult_Atrial	1
Human_Dermal_Fibroblast	2
Human_Dermal_Fibroblasts_Fetal	1
Human_Esophageal_Smooth_Muscle_Cells	1
Human_Intrahepatic_Biliary_Epithelial_Cells	1
Human_Periodontal_Ligament_Fibroblasts	1
Human_Skeletal_Muscle_Cells	1
Keratinocyte	2
Preadipocytes_visceral	1
Renal Proximal Tubular Epithelial_Cells	1
Renal_Cortical_Epithelial	1
Renal_Epithelial_Cells	1



Potential analysis (not a todo list, just possible ideas):
1. Attempt different methods of assigning a single methylation value to a gene (mentioned in the above message).
2. Differential methylation between our old neuronal samples and these new somatic primary cells (cell line to cell line comparison), or somatic tissue (cell line to normal tissue comparison)
3. Differential methylation between somatic tissue and somatic primary (normal tissue to cell line comparison)
4. Differential methylation between somatic cells and stem cells (mix to cell line comparison)
5. Differential methylation between hESCs and iPSCs (cell line to cell line comparison) <- this may actually be interesting...what differences are there, if any, between a stem cell and an iPS cell?
6. With whatever genes we found that have differential methylation, do they have differential expression? Is there any correlation with expression?
7. If multiple ways of looking at methylation for a gene, which correlates with expression the best?