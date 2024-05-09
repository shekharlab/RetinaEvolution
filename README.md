# Summary
This repository contains code written to analyze single cell RNA sequencing data for 17 vertebrate species associated with the paper //include citation here//. The analysis heavily relies on the R package [Seurat](https://github.com/satijalab/seurat) and may be useful to users who wish to reproduce the results of the paper or try alternative analysis strategies. The raw and processed sequencing data produced in this work are available via the Gene Expression Omnibus (GEO) under accession GSE237215. The species-specific datasets are available via the subseries accession numbers GSE237202–GSE237214. Previously published data utilized in this paper were downloaded from GEO repositories with accession numbers GSE81905, GSE137400, GSE152842, GSE148077, GSE15910 and GSE236005.

Please direct any questions associated with this repository to Joshua Hahn ([joshhahn@berkeley.edu](mailto:joshhahn@berkeley.edu)) or Karthik Shekhar ([kshekhar@berkeley.edu](mailto:kshekhar@berkeley.edu)). 

## Class and Type Annotation: An Exemplary Notebook Detailing the Analysis Pipeline
This notebook is designed to guide users through the process of annotating retinal cell classes, as well as annotating retinal ganglion cell (RGC) and bipolar cell (BC) types, using snRNA-seq data collected from tree shrew retinas. Steps to annotate cell classes include loading the count matrices, setting up the Seurat object, initial clustering, and annotation of clusters using canonical markers. Once annotated into cell classes, the data is split into separate RGC and BC objects and further processed using data integration, removal of contaminant cell classes, and cluster visualization. Many of these tasks are performed using functionalities within the R package Seurat.

## Species_Scripts: Individual Analysis for each Species
This folder contains notebooks used to analyze each of the 17 vertebrate species, with one notebook per species. The analysis pipeline for each species closely resembles that present in the Class and Type Annotation notebook, but we present them here for those who would like to reproduce the results.  

## Figure Scripts: Notebooks for Figure Reproducibility
This folder contains notebooks used to generate the figures in //Include citate here//. This includes generation of the cross-species integrated objects.

## utils: Custom Functions
This folder contains three scripts containing custom functions used for analysis.

### plottingFxns.R
Contains a variety of functions used for plotting and figure generation.
### utilFxns.R
Contains functions to condense portions of the analysis.
### xgboost_train.R
Contains functions for implementing the xgboost algorithm.


