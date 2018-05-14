# eRum 2018

## Interactivity meets Reproducibility: the `ideal` way of doing RNA-seq analysis

 My presentation for the eRum 2018 conference in Budapest.
 
 Submitted abstract:
 
    Next generation sequencing technologies, such as RNA-Seq, generate tens of 
    millions of reads to quantify the expression levels of the features of interest.
    A wide number and variety of software packages have been developed for accommodating
    the needs of the researchers, mostly in the R/Bioconductor framework. 
    
    Many of these focus on the identification of differentially expressed (DE) genes
    (DESeq2, edgeR, ...) to discover quantitative changes between experimental groups,
    while other address alternative splicing, discovery of novel transcripts, or RNA editing.
    
    Moreover, Exploratory Data Analysis is a common step to all these workflows, and 
    despite its importance for generating highly reliable results, it is often neglected,
    as many of the steps involved might require a considerable proficiency of the user 
    in the programming languages. Principal Components Analysis (PCA) is used often 
    to obtain a dimension-reduced overview of the data.
    
    Our proposal addresses the two steps of Exploratory Data Analysis and Differential
    Expression analysis with two different Bioconductor packages, pcaExplorer and ideal.
    We developed web applications in the Shiny framework also including support for 
    reproducible analyses, thanks to an embedded text editor and a template report, 
    to seamlessly generate HTML documents as a result of the user's exploration.
    
This repository contains the source (Rmd and md) and html version of the presentation, together with the css file to customize the appearance (making it somewhat Bioconductor-ish).

You can display this presentation here [`https://federicomarini.github.io/erum2018/`](https://federicomarini.github.io/erum2018/)
