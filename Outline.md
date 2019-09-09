---
title: "Bayesian Analyses In Paleontology: Interpreting the Posterior Sample"
output: html_document
author:
- Graeme T. LLoyd
- April M. Wright
bibliography: ../refs.bib
---

# Introduction:

1) Phylogenetic analyses of discrete morphological data have traditionally been conducted under the maximum parsimony criterion.
However, in recent times, Bayesian analyses have become more common for estimating phylogenetic trees, as well as more complex analyses, such as divergence time estimation.



Many researchers have attempted to formulate comparisons beween parsimony and likelihood tress.
  a) Maximum parsimony has typically been used to estimate phylogeny from discrete morphology
  b) A lot of baggage comes with this that does not strictly translate to Bayesian setting. For example, do not need an outgroup if tip-dating, synapomorphies are not strictly relevant as characters can change state and back again along branches, trees have a sampling frequency, consensus methds are therefore very different, there are no menanigful branch lengths in parsimony etc.

2) How does a Bayesian method differ from parsimony?
  a) Explicit assumption of a model
  b) Sample of trees is the goal, not something that sometimes occurs
    i) How does a posterior sample differ from a set of MPTs? What does this mean for
    interpretation?

3) How can we evaluate a Bayesian method?
  a) Intrinsic measures: those measures that require a model
    i) Simulation. We can simulate under the model and see if our data look like real data
    ii) Bayes factor model fitting. Of the models we have, which has the greatest suppoort
  b) Extrinsic: Measure that compare the estimated tree to some external information
    i) Stratigraphic congruence
    ii) What did they do in the Sansom paper?

# Methods

## Dataset Filtering

Datasets were downloaded from GraemeTLloyd.net. In the repository, there are some datasets that are derived from other datasets. For example, a matrix in the repository may have derivative matrices in which taxon or characters sampling was expanded from. For this study, if there were mutliple dataset derivatives, we chose the largest or most expansive dataset.


## Phylogenetic Analysis

We analysed each dataset under the Mk model in the software RevBayes.
We partitioned each dataset according to the number of states per character in order to specify an appropriate _Q_-matrix.
Among-character rate variation was parameterized according to a Gamma distribution with four discrete categories.
Branch lenghts were drawn from an exponential distribution.
Datasets were run for one million generations, and then assessed visually for convergence.

We did not re-estimate maximum parsimony trees for the data matrices.
For these datasets, we re-used the maximum parsimony trees published by the original study authors.

## Stratigraphic Congruence

We calculated stratigraphic congruence values for each tree in the Bayesian posterior sample, and each most parsimonious tree in the dataset.
For each dataset, we used the `PaleoDBQuerier` function in the `Claddis` R package to query the Paleobiology Database for minimum and maximum ages of each tip in the analysis.
Using the R package `Strap`, we calculated several stratigraphic congruence metrics: Minimum Implied Gap (MIG), Stratigraphic Consistency Index (SCI), Relative Completeness Index (RCI), Gap Excess Ratio (GER), Manhattan Stratigraphic Measure (MSM), and Wills' modifications of GER. An explanation of these metrics can be seen in Table 1.

| Metric | Meaning | Range |
|--------|---------|-------|
| Stratigraphic consistency index (SCI) | Proporion of nodes for which the oldest descendent of that node is younger than the oldest descendent of that node's ancestor| 0 to 1, with one being perfectly conssitent |
| Minimum Implied Gap (MIG) | The sum of the branch legnths without terminal tip ranges | Positive numbers |
|  Relative Completeness Index (RCI) | MIG score proportional to the total length of tip intervals | All real numbers |
| Manhattan Stratigraphic Measure (MSM)| MIG for the most stratigraphically consistent possible tree divided by the actual MIG | 0 to 1, with one being the most consistent |
| Gap Excess Ratio (GER) | MIG minus the best possible stratigraphic fit, scaled by the contrast between the best and worst fit values | 0 to 1, with one being the most consistent |

For each set of topologies and stratigraphic congruence results, we plotted a treespace plot in the R package `RWTY`.
`RWTY` calculates the topological differences between the trees in the posterior sample.
Then, it uses a multi-dimensional scaling algorithm to plot these differences to 2-D space.
We modified `RWTY` to color points according to their stratigraphic consistency score.
We visualized the set of most parsimonious trees and the Bayesian posterior sample on the same treespace plots.

# Results


# Discussion

1) How do we summarize the posterior?
  a) Distilling to a point estimate would be misleading
  b) Random samples (e.g. Sansom et al.) give a sense of variation, but lose the data related to proportional representation
  of good solutions

2) What does this mean for stratigraphic congruence - is it higher in Bayesian trees?

3) Plea for scientists coding reproducible figures and taking a clear look at variation in their data.
