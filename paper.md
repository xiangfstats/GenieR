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
output: pdf_document
bibliography: paper.bib
tags:
- R
- demographic history
- phylogeny
- coalescent models
affiliations:
- index: 1
  name: Department of Veterinary Medicine, University of Cambridge, Cambridge, United Kingdom
- index: 2
  name: US Military HIV Research Program, Rockville, MD, United States
---

# Summary

genieR is a R package (@R) for the inference of demographic history from reconstructed molecular phylogenies, seeing @Bethanycluster2017 for one application. The package processes phologentic trees in .tre format and allows users to extract basic information from a particular phylogeny such as sampling and coalescent times. Furthermore, genieR covers a few functions designed in GENIE, @PybusRambautGENIE, e. g.  the parameter estimation of a coalescent model (@kingman1982) when fitting on a particular phylogeny. One can fit 7 coalescent models including const (constant population size), expo (exponetial growth),expan (expansion growth), log (logistic growth), step (piecewise constant), pexpan (piecewise expansion growth) and plog (piecewise logistic growth) for his phylogenetic tree and compare models through AIC (@AIC1973) values. The package provides c++ code as an option for users to infer parameters in an efficient manner and also enables MCMC (Markov Chain Monte Carlo, @Hastings1970) inference on coalescent models. genieR can also simulate coalescent times for isochronous/heterochronous phylogenies.

# References
