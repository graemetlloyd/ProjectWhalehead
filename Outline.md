---
title: "Bayesian Analyses In Paleontology: Interpreting the Posterior Sample"
output: html_document
author:
- Graeme T. LLoyd
- April M. Wright
bibliography: ../refs.bib
---

# Introduction:

Inferring the tree of life, and using phylogenetic trees to understand organismal form and function is a prime challenge in biology.
Phylogenetic trees provide a historical perspective on relationships between taxa, and performing comparative biological studies without a phylogenetic context is misleading (Felsenstein 1985).
While many researchers in comparative biology have moved to using molecular data to estimate phylogenetic trees, fossils
remain the sole source of character data for many extinct groups.
Fossils provide the only direct evidence for historical trends in many groups (Cobbett et al. 2007; Koch and Parry in review).
 Excluding fossils from trees has been demonstrated to bias the inferences we draw from phylogenetic comparative analyses (Slater et al. 2012).
Therefore, how to best include fossils in phylogenetic analyses remains a vital and thriving research topic.

Phylogenetic analyses of discrete morphological data have traditionally been conducted under the maximum parsimony criterion.
However, in recent times, Bayesian analyses have become more common for estimating phylogenetic trees, as well as more complex analyses, such as divergence time estimation.
Bayesian methods assume an explicit model of evolution, and allow the researcher to set priors on each parameter in their model.
In effect, this allows researchers to incorporate prior information about the value of parameters into their analysis.
The most basic and restrictive model of morphological evolution is the Mk model, proposed in 2001 by Paul Lewis (Lewis 2001).
This model assumes equal rates of change between each state at a character.
In practice, these assumptions are not much different than the basic assumptions of an equal-weights parsimony step matrix.
However, these methods calculate a rate-based branch length in expected substitutions per site, allowing for the accommodations of multiple changes along a branch (e.g., a `0` state to a `1` state and back), something parsimony cannot do.
More elaborate models that allow for relaxation of the core assumptions of the Mk model have been published and evaluated (Nylander et al. 2004 [I assume this is the ref you mean?], Wright et al. 2016; OTHERS?).

Comparing the results of a Bayesian analysis and a parsimony analysis can be difficult.
Bayesian methods estimate what is called the _posterior sample_.
The posterior sample is a set of trees and associated model parameters that is plausible given the data, the model, and the prior.
Phylogenetic trees are often estimated using Markov-Chain Monte Carlo sampling.
Under this algorithm, a tree and model parameters are proposed and evaluated given the model and priors.
In general, if the tree and model parameters have a better score than the previous tree and model parameters, the new solution is kept.
MCMC is considered a memoryless process; that is, the next proposed tree and model parameters is not chosen based on previously-sampled trees.
A good tree may therefore be visited many times in a phylogenetic estimation.
In fact, how often a particular tree or set of relationships appears in the posterior sample is considered to be a measure of support for that tree.
Usually the researcher will compute a summary tree for publication, though tools now exist to visualize sets of phylogenetic trees.

This is in stark contrast to parsimony.
The aim of parsimony is to estimate the most parsimonious tree.
The most parsimonious tree is the one that minimizes the number of character changes in the dataset implied by the tree.
Under the criterion of maximum parsimony, a tree is proposed, and the changes implied in the morphological character matrix by that tree are counted.
A tree is considered better than another if it has a lower parsimony score.
Many datasets will have one tree that minimizes the parsimony score.
In this case, a point estimate of the phylogeny is typically published by the author.
However, in some datasets due to either lack of signal or conflicting signal, several most parsimonious trees may be returned.
In this case, a summary tree is usually created.
The most common of these is a consensus tree , which displays all phylogenetic relationships that are present in some proportion of the most parsimonious trees.
Common variations include a strict consensus tree, in which all clades on the summary tree are represented in all parsimonious trees, and majority rule consensus trees, in which clades on the summary tree are represented in more than 50\% of the set of most parsimonious trees.

How to compare a Bayesian posterior sample and a set of most parsimonious trees is a fraught topic.
There are many aspects of parsimony analyses that are not strictly comparable to Bayesian analyses.
For example, synapomorphies are interpreted differently in a Bayesian analysis, as multiple changes can occur along a branch.
Bayesian analyses consider branch lengths (in substitutions per site for an undated tree) to be parameters of the model.
Therefore, tree summarizations take these values into account (cite BEAST2, Puttick), not solely the topology.
Support values for bipartitions in the tree are calculated as part of a Bayesian estimation, being the number of times that a set of relationships is contained in the posterior sample.
The assumption is that the sample, not simply the `best` tree, contains vital information.

Because of these differences, formulating a comparison between parsimony and Bayesian trees is difficult.
Most studies to date have focused on intrinsic comparisons, those comparisons about the tree itself.
For example, most simulation studies to date have simulated data along a given tree or set of trees.
Then, from the simulation data, a tree has been estimated under both parsimony and Bayesian methods.
Finally, a summary tree for each method is computed and compared to the tree under which it was simulated.
Often, this focuses on the behavior of the researcher, comparing a Bayesian consensus tree to a parsimony consensus tree (self-flagellate).
Most phylogenetic estimates in published articles are presented as point estimates.
Because parsimony trees are most commonly published as majority-rule or strict consensus trees, computing this type of summary for both treatments (parsimony and Bayesian) and comparing them is fairly straightforward.
Comparisons have typically involved tree-based metrics, such as the Robinson-Foulds, which supplies the number of bifurcations that differ between two trees.
While this approach makes a degree of sense [shade: cast], it also does not include or account for most of the results associated with Bayesian estimation (the posterior sample).

Sansom et al (2017) used stratigraphic congruence to compare parsimony and Bayesian summary trees to assess which analytical method produces trees that best fit the fossil record.
The most attractive feature of this approach is that it employs extrinsic criteria, evaluating how well the tree fits data external to that used to infer the tree.
Stratigraphic congruence methods (Table One) use various measures to determine how well an estimated tree fits first and last appearance data for the fossil tips on this tree.
Sansom estimated trees for empirical datasets under both Maximum Parsimony and Bayesian estimation, computed consensus trees, and calculated stratigraphic fit for a sample of 500 trees from the posterior, then averaged those stratigraphic congruences.
This is a novel way to independently compare the fit of trees produced under different optimality criteria, although previouslty it has been used to compare competing hypotheses of relationships (Brochu and Norell 2000).
In this manuscript, we extend this approach to evaluate the whole posterior sample for stratigraphic fit and, critically, plot stratigraphic fit of individual trees in the posterior or the set of most parsimonious trees in "treespace" (sensu Hillis et al. 2006).
We conclude that averaging a posterior sample to a summary statistic, and encourage the use of more sophisticated visual tools, such as RWTY (Warren et al. 2017) to incorporate the variation in any one metric.

# Methods

## Dataset Filtering

Empirical datasets were downloaded from graemetlloyd.com/matr.html, the same starting repository used by Sansom et al (2017).
We initially excluded any data sets that were molecular, phenetic, ontogenetic or meta-analytical as these do not represent data sets intended for morphological phylogenetic inference (our focus here).
In the repository, and more broadly amongst morphological phylogenetics, many datasets are derived from older datasets with little or no modification.
For example, a matrix in the repository may have derivative matrices in which taxon or characters sampling was expanded from (see Figure S1 of Hartman et al. 2019 for an excellent illustration of the how complex this issue can become).
For this study, if there were multiple dataset derivatives, we chose the largest or most expansive dataset using the same criteria applied by Wright et al. (2016), employing the appropriate metadata available in the XML files associated with each data set.
The same XML files contain taxonomic reconciliations - links between OTU names used in the trees with taxa entered into the Paleobiology Database, the source used for our (and Sansom's) age data.
Using the API (Peters and McClennen 2015) we further filtered the data down to just those where all taxa are both reconciled and have age data available in the database (i.e., at least one fossil occurrence is attributable to each OTU).
The final matrices represent 128 different data sets.

## Phylogenetic Analysis

Each data set was analysed under both Maximum Parsimony and Bayesian criteria.

For the Bayesian approach we analysed each dataset under the Mk model in the software RevBayes.
We partitioned each dataset according to the number of states per character in order to specify an appropriate _Q_-matrix.
Among-character rate variation was parameterized according to a Gamma distribution with four discrete categories.
Branch lengths were drawn from an exponential distribution.
Datasets were run for one million generations, and then assessed visually for convergence using Tracer.

Maximum Parsimony trees were re-estimated for each data matrix using TNT (Goloboff and Catalano 2016).
Implicit enumeration was used where there were 24 or fewer tips.
For larger tip counts 20 replicates of new technology searches followed by a round of tree bisection-reconnection were applied with maxtreees capped at 100,000.

## Stratigraphic Congruence

We calculated stratigraphic congruence values for each tree in the Bayesian posterior sample, and each most parsimonious tree in the dataset.
For each dataset, we used the `PaleobiologyDBOccurrenceQuerier` function in the `metatree` R package (github.com/graemetlloyd/metatree) to query the Paleobiology Database for minimum and maximum ages of each tip in the analysis.
Using the R package `strap` (Bell and Lloyd 2015), we calculated several stratigraphic congruence metrics: Minimum Implied Gap (MIG), Stratigraphic Consistency Index (SCI), Relative Completeness Index (RCI), Gap Excess Ratio (GER), Manhattan Stratigraphic Measure (MSM*), and Wills' modifications of GER (GERt and GER*). An explanation of these metrics can be seen in Table 1.

| Metric | Meaning | Range |
|--------|---------|-------|
| Stratigraphic consistency index (SCI) | Proportion of nodes for which the oldest descendent of that node is younger than the oldest descendent of that node's ancestor| 0 to 1, with one being perfectly consistent |
| Minimum Implied Gap (MIG) | The sum of the branch lengths without terminal tip ranges | Positive numbers in millions of years |
| Relative Completeness Index (RCI) | MIG score proportional to the total length of tip intervals | All real numbers |
| Manhattan Stratigraphic Measure (MSM*)| MIG for the most stratigraphically consistent possible tree divided by the actual MIG | 0 to 1, with one being the most consistent |
| Gap Excess Ratio (GER) | MIG minus the best possible stratigraphic fit, scaled by the contrast between the best and worst fit values | 0 to 1, with one being the most consistent |

For each set of topologies and stratigraphic congruence results, we plotted a treespace plot in the R package `RWTY`.
`RWTY` calculates the topological differences between the trees in the posterior sample using the Robinson-Foulds metric.
This metric calculates the number of splits that are on one tree, but not the other.
Then, `RWTY` uses a multi-dimensional scaling algorithm to plot these differences to 2-D space.
We modified `RWTY` to color points according to their stratigraphic consistency score.
We visualized the set of most parsimonious trees and the Bayesian posterior sample on the same treespace plots.

The full set of code and data used are available at: github.com/graemetlloyd/ProjectWhalehead.

# Results

For direct comparison with Sansom et al (2017), we made boxplots of the distributions of several stratigraphic congruence measures
We made these figures comparing the posterior sample and the set of most parsimonious trees for each dataset, and each stratigraphic congruence metric shown on Table One.
Note that a major difference between a Bayesian sample and a parsimony sample is that the former can include duplicate trees (those visited multiple times by the MCMC, whereas a parsimony analysis only returns unique topolog(ies).
One such boxplot can be seen in Fig. 1A.
The spread of stratigraphic congruence metrics is much higher in the Bayesian analysis and this is typical of the majority of data sets.
For most metrics, the median of the parsimony sample is lower than the median of the Bayesian sample.
However, in most instances the stratigraphic congruence of the most parsimonious tree or set of trees is contained within the quartiles of the boxplot.

![](Figures/Fig1Labelled.pdf)

> Fig. 1: Plots showing Minimum Implied Gap score (MIG) for the posterior sample and set of most parsimnonious trees for the Yates (2003) dataset. Panel A shows boxplots of the distributios of MIG scores, as in Sansom (2018). Panel B shows the same data, but as a violin plot with points overlain. Panel C shows a treespace visualization, with points colored by MIG score. Large points indicate parsimony trees.

We also visualized these data as treespace plots.
An example treespace plot can be seen in Fig. 1C, and the full set of treespace plots is available in the supplemental material.
These multi-dimensional scaling graphs demonstrate that the Bayesian posterior samples contain many more trees, including many more trees that are substantially different from one another, than the parsimony estimations do.
In most estimations, the posterior sample contains the parsimony trees, in addition to other solutions plausible under the model.
As shown on Fig. 1C, it is very possible (and even common) in both Bayesian and maximum parsimony estimation to have topologies with good stratigraphic fit plotted near trees with poor stratigraphic fit (and for dissimilar trees to have similar or even identical stratigraphic fits).
This indicates that in some treespaces, there may be little topological difference between a tree that is quite good with respect to stratigraphic fit and a tree with poor stratigraphic fit.
In other words, the landscapes of stratigraphic congruence is generally highly volatile - something boxplots fail to capture.
As in the boxplots, the treespace plots largely indicate that the parsimony trees fall within the range of solutions explored in a Bayesian search.

![](Figures/Fig2Labelled.pdf)

> Fig. 2: Plots showing Minimum Implied Gap score (MIG) for the posterior sample and set of most parsimnonious trees for the Demar (2013) dataset. Plot types and labels are the same as in Fig. 1. In contrast to Fig. 1, there are defined regions of treespace in which trees with better MIG scores are found.

# Discussion

## Summarizing a posterior sample

The aim of a Bayesian analysis is not a single point estimate of a solution.
Rather, it is to examine solutions and outcomes that are plausible given the model and the data.
This is particularly important in phylogenetics, where we are estimating lineage divergences that occurred tens or hundreds of millions of years ago.
We do this from scarce, and likely biased data, using models that may or may not adequately capture reality.
To responsibly present a solution under these conditions must mean incorporating uncertainty in that solution.
It is prudent to avoid using a single point estimate summary of a posterior sample, whether a summary tree or an average value computed across several trees.

The authors of this work are not immune to using a point estimate of topology from a posterior sample to compare to a point estimate of a parsimony topology (see Wright and Hillis 2014).  
But the properties of a Bayesian analysis are different than a parsimony analysis.
A Bayesian analysis represents solutions that are plausible under a model.
A parsimony analysis, on the other hand, pre-filters trees by whether or not they are the best under the maximum parsimony criterion.
This means that the Bayesian sample will almost certainly contain more sub-optimal trees than the parsimony criterion.
This has the effect, when stratigraphic congruence measures are averaged for the posterior sample, of pulling down the mean of the distribution of stratigraphic consistency values.

There is an additional complication to interpreting a Bayesian posterior distribution: that the distribution itself is meaningful.
Because a Markov Chain Monte Carlo analysis is memoryless, and solutions are sampled in proportion to their likelihood and the prior, a good solution will be visited many times.
In fact, popular metrics of phylogenetic support, such as the posterior probability are calculated directly from the posterior sample for this reason.
That the distribution of solutions itself contains information means that approaches, such as randomly sampling trees from the posterior, may give some sense of the variation.
But they can also fail to approximate the true distribution of the posterior, particularly given that the variance in a posterior sample is unlikely to be normally distributed.

However, it is also not reasonable for researchers to examine every single tree in a posterior sample.
How can a researcher draw reasonable conclusions from their posterior sample?
Over the past five years, modern and open-source treespace visualization software has become available.
As demonstrated in Fig. 1C (a treespace plot made from the posterior sample estimated from the Yates (2003) dataset), treespace visualizations show the distance between trees in the sample, and enable researchers to color points by other factors (formally generating so-called tree "landscapes"; St John 2017).
As can be seen in Fig. 1c, both the Bayesian and parsimony analyses predominantly sampled two islands of tree space.
The points in this graphic show that topologically similar trees can have wildly different Minimum Implied Gap scores.
For this dataset, it does not appear that certain regions of treespace produce better stratigraphic congruence.
This is in contrast to Fig. 2C (a treespace plot made from the posterior sample estimated from the Demar (2013) dataset), in which some areas of treespace clearly contain trees that have better fit to the stratigraphic record.

## Stratigraphic Congruence and optimality criterion

Which analytical method produces better stratigraphic congruences is an all together murkier question when the fullness of the posterior sample is considered.
A table of averages for Yates (2003) and Demar (2013) is provided on Table Two.
As can be seen on this table for both the Demar (2013) and Yates (2003) datasets, some average metrics are better for the Bayesian posterior sample, some for the set of most parsimonious trees.
Overall, 54 of 128 datasets had lower average SCI using Bayesian methods (Table S1). [So this section is I guess a somewhat random accounting of the strat congruence metrics. The "big four" are probably SCI, MSM, RCI and GER. I dislike SCI as it treats the data as discrete rather than continuous so all gaps are "equal" when really they aren't. The other three are all just driven by MIG and so this would be my preferred metric, even though it isn't strictly a metric.][#AMW: can you write this in to the discussion? You know SC metrics much better than I]
36 had better SCI values under parsimony, and in 38 datasets the values were identical.
On the other hand, for the RCI, 62 datasets had better average scores under parsimony.

An average may not be providing a good accounting of the variation in the results for each dataset.
In Fig. 2A and 2B, we show a dataset from Demar et al (2013).
The treespace for this dataset can be seen in Fig. 2C.
This is a highly structured treespace: because Bayesian MCMC samples in proportion to the posterior probability, we can infer from the shape of this treespace visualization that there are two peaks to this distribution that contained fairly good trees.
There are three most-parisimonious trees; two are sampled from one peak with poor stratigraphic fit, one is sampled from the peak with fairly good stratigraphic fit.
In this case, averaging is unlikely to produce a value that represents either peak adequately, and the average in this case is a value that doesn't belong to a tree that was sampled in the analysis at all.

In the case of the Yates (2003) dataset, treespace is sampled much more evenly.
In this case, taking an average is likely to represent the sample a bit better.
Even so, there are a number of problems with comparing means.
As shown in the boxplots of Fig. 1, the mean and variance for the parsimony estimates are almost wholly subsumed in the Bayesian posterior sample.
This is largely expected: the Bayesian posterior sample encompasses all trees sampled in the analysis (typically thinned by some proportion as they are exported from the tree estimation software).
The set of most parsimonious trees contains a sample of the best trees according to an optimality criterion (fewest implied evolutionary steps).
Assuming both estimation criteria are using the same data to sample the same treespace, we would expect that averaging among the best trees should produce better stratigraphic congruence.
However, looking across all the datasets, only a slight improvement is seen from this biased averaging (Table S1, Sansom et al. Fig. 1).

Which analytical method produces better stratigraphic congruence is the wrong question.
A better question to ask may be "How can researchers explore and quantify variation in their sets of solutions?".
Here we argue a key but underappreciated visualization tool (treespace) may be especially useful.
Implementations of such tools have been in place for well over a decade (Hillis et al. 2005), but we are unaware of them being used in any published analysis of samples of trees produced from morphological data.
We also note that we see many different "kinds" of treespace across the data sets examined here.
For example, single or multiple-tree islands, smooth or volatile gradients, parsimony samples enveloped by Bayesian samples, or parsimony and Bayesian samples occupying different parts of the space.

In writing this manuscript, we used and modified a set of tools available in the R programming language. [This seems like a swerve so made it a separate paragraph. Might want a header as well though - Conclusion?]
These tools, such as the package `RWTY`, enabled us to read in many large posterior samples, and to calculate treespace plots across 128 empirical datasets.
The code and data available for this are freely available on GitHub. [Can prolly delete as we say this earlier?]
We would like to close this manuscript by noting that the tools to perform complex formatting and subsetting datasets, including large ones such as a Bayesian posterior sample, are more accessible than ever.
We look forward to many years of interesting analyses about how different methods explore and sample phylogenetic tree space.




Table Two [I would maybe cut this down. I feel like MIG is all that really matters and maybe we can just summarise the range/varaince of this?]

|Tree | Mean MIG   |  Max MIG | Min MIG | 
|-----|------------|----------|---------|
| Demar - Bayesian |  219     |    248  |  177 | Demar - Parsimony |  197 |    255 |    161 
| Yates - Bayesian  | 149  |    302 |    57 
| Yates - Parsimony | 178  | 305 | 128 |



# References

Bell, M. A. and Lloyd, G. T., 2015. strap: an R package for plotting phylogenies against stratigraphy and assessing their stratigraphic congruence. Palaeontology, 58, 379-389.

Brochu, C. A. and Norell, M. A., 2000. Temporal congruence and the origin of birds. Journal of Vertebrate Paleontology, 20, 197-200.

Cobbett, A., Wilkinson, M. and Wills, M. A., 2007. Fossils impact as hard as living taxa in parsimony analyses of morphology. Systematic Biology, 56, 753-766.

Felsenstein, J., 1985. Phylogenies and the comparative method. American Naturalist, 125, 1-15.

Goloboff, P. A. and Catalano, S. A., 2016. TNT version 1.5, including a full implementation of phylogenetic morphometrics. Cladistics, 32, 221-238.

Hartman, S., Mortimer, M., Wahl, W. R., Lomax, D. R., Lippincott, J. and Lovelace, D. M., 2019. A new paravian dinosaur from the Late Jurassic of North America supports a late acquisition of avian flight. PeerJ, 7, e7247.

Hillis, D. M., Heath, T. A. and St John, K., 2005. Analysis and visualization of tree space. Systematic Biology, 54, 471–482.

Koch, N. M. and Parry, L. A., in review. Death is on our side: paleontological data drastically modify phylogenetic hypotheses. https://www.biorxiv.org/content/10.1101/723882v1

Lewis, P. O., 2001. A likelihood approach to estimating phylogeny from discrete morphological character data. Systematic Biology, 50, 913-925.

Nylander, J. A., Ronquist, F., Huelsenbeck, J. P., Nieves-Aldrey, J. L., 2004. Bayesian phylogenetic analysis of combined data. Systematic Biology, 53, 47-67.

Peters, S. E. and McClennen, M., 2015. The Paleobiology Database application programming interface. Paleobiology, 42, 1-7.

Sansom, R. S., Choate, P. G., Keating, J. N. and Randle, E., 2018. Parsimony, not Bayesian analysis, recovers more stratigraphically congruent phylogenetic trees. Biology Letters, 14, [????].

Slater, G. J., Harmon, L. J. and Alfaro, M. E., 2012. Integrating fossil with molecular phylogenies improves inference of trait evolution. Evolution, 66, 3931-3944.

St. John, K., 2017. The shape of phylogenetic treespace. Systematic Biology, 66, e83-e94.

Warren, D. L., Geneva, A. J. and Lanfear, R., 2017. RWTY (R We There Yet): an R package for examining convergence of Bayesian phylogenetic analyses. Molecular Biology and Evolution, 34, 1016-1020.

Wright, A. M., Lloyd, G. T. and Hillis, D. M., 2016. Modeling character change heterogeneity in phylogenetic analyses of morphology through the use of priors. Systematic Biology, 65, 602-611.
