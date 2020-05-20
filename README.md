# Project Whalehead

## Introduction

Project Whalehead (arbitrarily named for the shoebill genus) is a GitHub repository for a project comparing phylogenetic outputs of Bayesian and parsimony analysis of morphological character matrices using stratigraphic congruence metrics. The authors are [April Wright](https://github.com/wrightaprilm) and [Graeme Lloyd](https://github.com/graemetlloyd) and this repository also serves as supplementary information for the final published version.

## What is here?

Below is a summary of the contents:

* Data (Folder containing the data used)
    - Ages (Tip age text files generated for each data set primarily using calls to the Paleobiology Database API)
    - Boxplots (Boxplots showing the distribution of stratigraphic congruence metrics for Bayesian and parsimony trees calculated in the strap R package for each data set)
    - Combined trees (The combined Bayesian and parsimony trees in NEXUS format)
    - kMeans_MIG (k-means clusters based on MIG scores)
    - kMeansTreespace (kmeans clusters based on Robinson-Foulds score)
    - MIGDifferences (Bivariate plots showing distance in MDS space verus difference in MIG for just the first 2 MDS axes and for all MDS axes)
    - MPTs (The most parsimonious trees generated from TNT searches for each data set in Newick format)
    - NEXUS (The original morphological character matrix in NEXUS format as taken from graemetlloyd.com)
    - RevBayesOutput (The raw RevBayes output for each data set; NB: Missing from GitHub as too large to upload, but called in scripts)
    - StratCongruenceResults (The stratigraphic congruence results output from the strap R package for each data set)
    - TreeSpacePlots (Treespace visualisations using modified RWTY R package code for each data set)
    - ViolinPlots (Violin plots showing the distribution of stratigraphic congruence metrics for Bayesian and parsimony trees calculated in the strap R package for each data set)
    - XML (XML files for each data set containing metadata used to generate pseudo-independent data sets and reconcile OTU names with Paleobiology Database taxa)
* Figures (Final and draft versions of manuscript figures)
* R (R package folder for replicating study)
* man (R package folder for replicating study)
* rb_scripts (RevBayes scripts used in Bayesian analysis)
* vignettes (Vignettes for replicating analysis)
* DESCRIPTION (R package file for replicating study)
* NAMESPACE (R package file for replicating study)
* README.md (This file)
* TableS1.md (Supplementary Table 1 summarising full results for each data set)
* refs.bib (BibTex file containing references cited)
