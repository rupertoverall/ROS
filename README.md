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

## Methods


## Additonal packages required

* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) - Differential expression of sequencing data
* [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html) - Functional enrichment of Gene Ontology terms
* [ReactomePA](http://bioconductor.org/packages/release/bioc/html/ReactomePA.html) - Functional enrichment of Reactome terms
* [XLConnect](https://cran.r-project.org/web/packages/XLConnect/index.html) - Writing Excel spreadsheets for convenient output of functional enrichment results
* [VennDiagram](https://cran.r-project.org/web/packages/VennDiagram/index.html) - Venn diagrammes

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

