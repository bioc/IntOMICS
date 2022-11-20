# IntOMICS: an R package for integrative analysis of multi-omics data to infer regulatory networks.

IntOMICS is an efficient integrative framework based on Bayesian networks.
IntOMICS systematically analyses gene expression (GE), DNA methylation (METH), copy number variation (CNV) and biological prior knowledge (B) to infer regulatory networks. 
IntOMICS complements the missing biological prior knowledge by so-called empirical biological knowledge (empB), estimated from the available experimental data. 
An automatically tuned MCMC algorithm ([Yang and Rosenthal, 2017](http://probability.ca/jeff/ftpdir/jinyoung1.pdf)) estimates model parameters and the empirical biological knowledge.
Conventional MCMC algorithm with additional Markov blanket resampling (MBR) step ([Su and Borsuk, 2016](https://jmlr.org/papers/volume17/su16a/su16a.pdf)) infers resulting regulatory network structure consisting of three types of nodes: GE nodes refer to gene expression levels, CNV nodes refer to associated copy number variations, and METH nodes refer to associated DNA methylation probe(s).

IntOMICS takes as input: 
* gene expression matrix (required), 
* associated copy number variation matrix sampled from the same individuals (optional), 
* associated DNA methylation matrix of beta-values sampled from the same individuals (optional), and 
* the biological prior knowledge with information on known interactions among molecular features (optional, highly recommended). 

The resulting regulatory network structure contains the edge weights $`w_i`$ representing the empirical frequency of given edge over samples of network structures from two independent MCMC simulations.

<p align="center">
  <img src="vignettes/figures/IntOMICS_framework_METH_empB_modules.png" width="400" height="680" alt="IntOMICS framework">
</p>
     
For further details about the IntOMICS algorithm, its performance and benchmark analysis, see manuscript [Pacinkova \& Popovici, 2022](https://assets.researchsquare.com/files/rs-1291540/v1_covered.pdf?c=1643735189). 

# Installation
```ruby
# bioconductor install
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("IntOMICS")

# install the newest (development) version from GitHub
# install.packages("remotes")
remotes::install_github("anna-pacinkova/IntOMICS")

```
