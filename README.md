# ROS

Code used for the next-generation RNA sequencing analysis accompanying the manuscript "**Endogenous Redox levels delineate functional heterogeneity and responsiveness of adult hippocampal stem cells in mice**".
Two experiments are described; 
* A comparison of the two neurogenic niches of the adult mouse hippocampus.
* Characterisation of the Nestin-positive neural precursor cells sorted on the basis of ROS (reactive oxygen species) content.
All analysis has been performed using R/BioConductor.

Raw sequencing data can be found at [GEO](https://www.ncbi.nlm.nih.gov/geo/) under the SuperSeries accession [GSE124095](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124095).

## Included files
There are two scripts containing helper code;
* *Extras.R* several useful custom functions that simplify the main code
* *Signatures.R* a workflow for calculating 'enriched' and 'signature' genes and a standard functional enrichment pipeline we have assembled

### Data
There are two .RData files each containing a list with the following data;
* *enriched* a list of genes classed as 'enriched' in each of the sample groups (see Methods below for details)
* *signatures* a list of genes classed as 'signatures' of each of the sample groups (see Methods below for details)
* *top.tags* the full output of the *topTags()* command from *edgeR*
* *cpm* counts per million after filtering
* *Group* a factor describing the sample groups
* *group.colours* some colours we used to identify the groups
* *headers* a data.frame containing the full set of factors and metadata describing tha samples
* *counts* the raw counts
* *report" some additional information about the number of transcripts before and after quality control as well as detailed information about the session in which this data freeze was produced


### Figures
Scripts to perform the analyses in figures 2 and 4 (and some supplementary data) are presented in their own folders.

### Resources
*
*

## Methods
### RNA extraction 
Two separate RNA sequencing experiments were performed in this study. In the first experiment, tissue was microdissected from the subventricular zone of the lateral ventricle and the dentate gyrus of the same animal following the protocol of Walker et al. (2014). RNA was extracted using the RNeasy micro kit (Qiagen). In the second experiment, Nes-GFP positive cells from the DG were gated into ROS classes, based on their intracellular ROS content. 400 GFP positive cells from each ROS class were FACsorted into a PCR tube containing 8.5 μl of a hypotonic reaction buffer and RNA was prepared using the SMARTer Ultra Low RNA HV Kit (Takara Bio) according to the manufacturer’s protocol. For both experiments, cDNA of polyadenylated mRNA was synthesized from RNA of the lysed cells using SmartScribe reverse transcriptase, a universally tailed poly-dT primer and a template switching oligonucleotide (Takara Bio). This was followed by 12 cycles of amplification of the purified cDNA with the Advantage 2 DNA Polymerase. After ultrasonic shearing of the amplified cDNA (Covaris S2), samples were subjected to standard Illumina fragment library preparation using the NEBnext Ultra DNA library preparation chemistry (New England Biolabs). In brief, cDNA fragments were end-repaired, A-tailed and ligated to indexed Illumina Truseq adapters. Resulting libraries were PCR-amplified for 15 cycles using universal primers, purified using XP beads (Beckman Coulter) and then quantified with the Fragment Analyzer. Final libraries were equimolarly pooled and subjected to 75-bp-single-end sequencing on the Illumina HiSeq2500 platform, providing ~35 (24–60) million reads per sample. Reads were mapped to the latest mouse genome build (mm10) using the STAR algorithm (Dobin et al., 2013) and counts per ENSEMBL gene model prepared using the RSubread package {Liao:2013fo} for R/BioConductor.

### Quality control and differential expression
RNASeq counts were filtered to have at least 1 count per million reads (CPM) in a minimum of 75 % of the samples from at least one cell population. CPM were calculated using the function cpm from the edgeR package in R/BioConductor (Robinson et al., 2010). Samples were clustered by plotting the first two principal components and by unsupervised hierarchical clustering. One sample from each of the two experiments showed reduced sequencing depth and complexity and did not cluster with replicates. In both cases, these samples came from preparations with very low input RNA so we decided to remove them from further analyses. A filter for differential expression was then performed using the edgeR functions lmFit and topTags and only significantly (adjusted p < .05) differentially expressed transcripts were used for enrichment analyses. Transcripts were classified as 'enriched' in an experimental group if the mean expression in that group was significantly above the average expression over all groups. Enriched genes were further filtered by hierarchical clustering into two clusters based on inter-group t-statistics and those where only one group clustered separately were termed 'signature' genes (for the DG/SVZ comparison, these are obviously equivalent to the 'enriched' gene set).

### Functional enrichment and expression profiles
Enrichment for Gene Ontology terms was performed using the R package topGO (Alexa et al., 2018) with the DG/SVZ ‘enriched’ gene lists as query sets and all genes passing the CPM filter (see above) as background. Expression profiles of curated gene lists were calculated by identifying all genes in the current data that corresponded to the genes in the gene list of interest and calculating the first principal component to collapse their expression into an ‘eigengene’. The eigengene values were then used for plotting. A similar approach was used to reanalyse the Shin et al., dataset where the genes corresponding to the ‘signature’ genes for each ROS class were identified in the single-cell dataset and the first principal component of these used to create an eigengene. Smooth spline interpolation of the eigengene was then performed to produce the plotted values.


## Additonal packages required

* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) - Differential expression of sequencing data
* [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html) - Functional enrichment of Gene Ontology terms
* [ReactomePA](http://bioconductor.org/packages/release/bioc/html/ReactomePA.html) - Functional enrichment of Reactome terms
* [XLConnect](https://cran.r-project.org/web/packages/XLConnect/index.html) - Writing Excel spreadsheets for convenient output of
functional enrichment results
* [readxl](https://cran.r-project.org/web/packages/readxl/index.html) - Efficient reading of Excel spreadsheets (but not writing)
* [org.Mm.eg.db](http://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html) - Mouse gene annotation data
* [AnnotationDbi](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html) - Annotation tools
* [VennDiagram](https://cran.r-project.org/web/packages/VennDiagram/index.html) - Venn diagrammes
