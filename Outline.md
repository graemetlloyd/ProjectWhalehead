---
title: "Bayesian Analyses In Paleontology: Interpreting the Posterior Sample"
output: html_document
author:
- Graeme T. LLoyd
- April M. Wright
bibliography: ../refs.bib
---

# Introduction:

Inferring the tree of life, and using phylogenetic trees to understand organismal form and function [broader? Like anything above the species-level needs a tree right?] is a prime challenge in biology.
Phylogenetic trees provide a historical perspective on relationships between taxa, and performing comparative biological studies without a phylogenetic context is misleading (cite Felsenstein).
While many researchers in comparative biology have moved to using molecular data to estimate phylogenetic trees, fossils provide the only direct evidence for historical trends in many groups (Cobbett et al. 2007; Koch and Parry in review) and excluding fossils from trees has been demonstrated to bias the inferences we draw from phylogenetic comparative analyses (cite Slater).
Therefore, how to best include fossils in phylogenetic analyses remains a vital and thriving research topic.

Phylogenetic analyses of discrete morphological data have traditionally been conducted under the maximum parsimony criterion.
However, in recent times, Bayesian analyses have become more common for estimating phylogenetic trees, as well as more complex analyses, such as divergence time estimation.
Bayesian methods assume an explicit model of evolution, and allow the researcher to set priors on each parameter in their model.
In effect, this allows researchers to incorporate prior information about the value of parameters into their analysis.
The most basic and restrictive model of morphological evolution is the Mk model, proposed in 2001 by Paul Lewis.
This model assumes equal rates of change between each state at a character.
In practice, these assumptions are not much different than the basic assumptions of an equal-weights parsimony step matrix.
However, these methods calculate a rate-based branch length, allowing for the accommodations of multiple changes along a branch (e.g., 0 -> 1 -> 0), something parsimony cannot do.
More elaborate models that allow for relaxation of the core assumptions of the Mk model have been published and evaluated (Wright et al. 2016; OTHERS?).

Comparing the results of a Bayesian analysis and a parsimony analysis can be difficult.
Bayesian methods estimate what is called the _posterior sample_.
The posterior sample is a set of trees and associated model parameters that is plausible given the data, the model, and the prior.
Phylogenetic trees are often estimated using Markov-Chain Monte Carlo sampling.
Under this algorithm, a tree and model parameters are proposed and evaluated given the model and priors.
In general, if the tree and model parameters have a better score than the previous tree and model parameters, the new solution is kept.
MCMC is considered a memoryless process; that is, the next proposed tree and model parameters is not chosen based on previously-sampled trees.
A good tree may therefore be visited many times in a phylogenetic estimation.
In fact, how often a particular tree or set of relationships appears in the posterior sample is considered to be a measure of support for that tree.
Usually, the researcher will compute a summary tree for publication, though tools now exist to visualize sets of phylogenetic trees.

This is in stark contrast to parsimony.
The aim of parsimony is to estimate the most parsimonious tree.
The most parsimonious tree is the one that minimizes the number of character changes in the dataset implied by the tree.
Under the criterion of maximum parsimony, a tree is proposed, and the changes implied in the morphological character matrix by that tree are counted.
A tree is considered better than another if it has a lower parsimony score.
Many datasets will have one tree that minimizes the parsimony score.
In this case, a point estimate of the phylogeny is typically published by the author.
However, in some datasets due to either lack of signal or conflicting signal, several most parsimonious trees may be returned.
In this case, a summary tree is usually created.
The most common of these is a majority-rule consensus tree (is it? I feel ike strict is common too), which displays all phylogenetic relationships that are present in more than 50\% of the set of most parsimonious trees.

How to compare a Bayesian posterior sample and a set of most parsimonious trees is a fraught topic.
There are many aspects of parsimony analyses that are not strictly comparable to Bayesian analyses.
For example, synapomorphies are interpreted differently in a Bayesian analysis, as multiple changes can occur along a branch.
Bayesian analyses consider branch lengths (in substitutions per site for an undated tree) to be parameters of the model.
Therefore, tree summarizations take these values into account, not solely the topology.
Support values for bipartitions in the tree are calculated as part of a Bayesian estimation, being the number of times that a set of relationships is contained in the posterior sample.
The assumption is that the sample, not simply the `best` tree, contains vital information.

Formulating a comparison between parsimony and Bayesian trees is difficult.



3) How can we evaluate a Bayesian method?
  a) Intrinsic measures: those measures that require a model
    i) Simulation. We can simulate under the model and see if our data look like real data
    ii) Bayes factor model fitting. Of the models we have, which has the greatest suppoort
  b) Extrinsic: Measure that compare the estimated tree to some external information
    i) Stratigraphic congruence
    ii) What did they do in the Sansom paper?

# Methods

## Dataset Filtering

Datasets were downloaded from graemetlloyd.com/matr.html. In the repository, there are some datasets that are derived from other datasets. For example, a matrix in the repository may have derivative matrices in which taxon or characters sampling was expanded from. For this study, if there were mutliple dataset derivatives, we chose the largest or most expansive dataset using the same criteria applied by Wright et al. (2016). We then excluded molecular data sets and data sets with polymorphisms or uncertainties.

## Phylogenetic Analysis

We analysed each dataset under the Mk model in the software RevBayes.
We partitioned each dataset according to the number of states per character in order to specify an appropriate _Q_-matrix.
Among-character rate variation was parameterized according to a Gamma distribution with four discrete categories.
Branch lenghts were drawn from an exponential distribution.
Datasets were run for one million generations, and then assessed visually for convergence.

We did not re-estimate maximum parsimony trees for the data matrices. (Oh, I totally did do this. But in theory they are the same.)
For these datasets, we re-used the maximum parsimony trees published by the original study authors.

## Stratigraphic Congruence

We calculated stratigraphic congruence values for each tree in the Bayesian posterior sample, and each most parsimonious tree in the dataset.
For each dataset, we used the `PaleobiologyDBOccurrenceQuerier` function in the `metatree` R package (github.com/graemetlloyd/metatree) to query the Paleobiology Database for minimum and maximum ages of each tip in the analysis.
Using the R package `strap` (Bell and Lloyd 2015), we calculated several stratigraphic congruence metrics: Minimum Implied Gap (MIG), Stratigraphic Consistency Index (SCI), Relative Completeness Index (RCI), Gap Excess Ratio (GER), Manhattan Stratigraphic Measure (MSM*), and Wills' modifications of GER (GERt and GER*). An explanation of these metrics can be seen in Table 1.

| Metric | Meaning | Range |
|--------|---------|-------|
| Stratigraphic consistency index (SCI) | Proporion of nodes for which the oldest descendent of that node is younger than the oldest descendent of that node's ancestor| 0 to 1, with one being perfectly consistent |
| Minimum Implied Gap (MIG) | The sum of the branch legnths without terminal tip ranges | Positive numbers in millions of years |
|  Relative Completeness Index (RCI) | MIG score proportional to the total length of tip intervals | All real numbers |
| Manhattan Stratigraphic Measure (MSM*)| MIG for the most stratigraphically consistent possible tree divided by the actual MIG | 0 to 1, with one being the most consistent |
| Gap Excess Ratio (GER) | MIG minus the best possible stratigraphic fit, scaled by the contrast between the best and worst fit values | 0 to 1, with one being the most consistent |

For each set of topologies and stratigraphic congruence results, we plotted a treespace plot in the R package `RWTY`.
`RWTY` calculates the topological differences between the trees in the posterior sample.
Then, it uses a multi-dimensional scaling algorithm to plot these differences to 2-D space.
We modified `RWTY` to color points according to their stratigraphic consistency score.
We visualized the set of most parsimonious trees and the Bayesian posterior sample on the same treespace plots.

The full set of code and data used are available at: github.com/graemetlloyd/ProjectWhalehead.

# Results


# Discussion

1) How do we summarize the posterior?
  a) Distilling to a point estimate would be misleading
  b) Random samples (e.g. Sansom et al.) give a sense of variation, but lose the data related to proportional representation
  of good solutions

2) What does this mean for stratigraphic congruence - is it higher in Bayesian trees?

3) Plea for scientists coding reproducible figures and taking a clear look at variation in their data.

# References

Bell, M. A. and Lloyd, G. T., 2015. strap: an R package for plotting phylogenies against stratigraphy and assessing their stratigraphic congruence. Palaeontology, 58, 379-389.

Cobbett, A., Wilkinson, M. and Wills, M. A., 2007. Fossils impact as hard as living taxa in parsimony analyses of morphology. Systematic Biology, 56, 753-766.

Koch, N. M. and Parry, L. A., in review. Death is on our side: paleontological data drastically modify phylogenetic hypotheses. https://www.biorxiv.org/content/10.1101/723882v1

Wright, A. M., Lloyd, G. T. and Hillis, D. M., 2016. Modeling character change heterogeneity in phylogenetic analyses of morphology through the use of priors. Systematic Biology, 65, 602-611.

