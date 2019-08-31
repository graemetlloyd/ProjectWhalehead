---
title: "Bayesian Analyses In Paleontology: Interpreting the Posterior Sample"
output: html_document
author:
- Graeme T. LLoyd
- April M. Wright
bibliography: ../refs.bib
---

# Introduction:

1) Historical perspective on discrete morphology analyses
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

1) Reanalysis under the Mk model in RevBayes
2) Stratigraphic congruence across the whole sample
  a) Number of unique solutions is interesting
  b) But so if if the posterior is biased towards or away from good SC values. A good
  solution will be hit many times
3) SC on summary trees

# Discussion

1) How do we summarize the posterior?
  a) Distilling to a point estimate would be misleading
  b) Random samples (e.g. Sansom et al.) give a sense of variation, but lose the data related to proportional representation
  of good solutions

2) What does this mean for stratigraphic congruence - is it higher in Bayesian trees?

3) Plea for scientists coding reproducible figures and taking a clear look at variation in their data.
