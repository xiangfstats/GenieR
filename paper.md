---
title: 'genieR: An R package for inference of demographic history of phylogenies'
authors:
- affiliation: 1
  name: Fei Xiang
  orcid: 0000-0003-3644-7114
- affiliation: 2
  name: Bethany Dearlove
  orcid: 0000-0003-3653-4592
- affiliation: 1
  name: Simon Frost
  orcid: 0000-0002-5207-9879
date: "07 January 2018"
output:
  pdf_document: default
  html_document: default
bibliography: paper.bib
tags:
- R
- demographic history
- phylogeny
- coalescent models
affiliations:
- index: 1
  name: Department of Veterinary Medicine, University of Cambridge, Cambridge, United
    Kingdom
- index: 2
  name: US Military HIV Research Program, Rockville, MD, United States
---

# Summary

genieR is an R package (@R) for the inference of demographic history from reconstructed molecular phylogenies (see @Bethanycluster2017 for one application), modeled after the C++ package GENIE [@PybusRambautGENIE]. The package processes phylogenetic trees in Newick format, in which taxa may be sampled either at the same time (isochronous) or at different times (heterochronous); the latter are particularly common in the analysis of pathogen sequence data. In addition to extracting basic information from a particular phylogeny such as sampling and coalescent times, genieR also fits a coalescent model (@kingman1982) under several demographic scenarios, including constant population size, exponential growth, expansion growth, logistic growth, with both continuous and piecewise variants. The fit of different models can be compared using Akaike's Information Criterion (AIC) [@AIC1973] values. genieR uses C++ code as an option for better computational performance, and models can be fit either using maximum likelihood or Markov Chain Monte Carlo [@Hastings1970]. genieR can also simulate phylogenies under different demographic models.

# References
