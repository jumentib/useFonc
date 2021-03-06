% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mediation.R
\name{r_mediation_real_cell_type}
\alias{r_mediation_real_cell_type}
\title{r_mediation_real_cell_type : function to simulate DNA methylation data for mediation analyzis (and cell type)}
\usage{
r_mediation_real_cell_type(CT = NULL, CT.l = NULL, n = nrow(CT),
  p = ncol(CT.l), K = 2, K.ct = 5, freq = NULL,
  prop.causal.x = 0.01, prop.causal.y = 0.01, prop.causal.ylx = 0.5,
  prop.variance.y = 0.1, prop.variance.x = 0.1, rho = 0.1,
  sigma = 1, sd.A = 0.1, mean.A = 1, sd.B = 0.1, mean.B = 1,
  sd.U = 1, sd.V = 1, strength = 1, sd.ct = 1, alpha = NULL,
  prob.bin = NULL)
}
\arguments{
\item{CT}{: Cell type proportion. null by default. You can use real data (n * K.ct matrix)}

\item{CT.l}{: Cell type loading. null by default. You can use real data (K.ct * p matrix)}

\item{n}{: number of individuals}

\item{p}{: number of cpg variables}

\item{K}{: number of latent factors}

\item{K.ct}{: number of cell type}

\item{freq}{: (vector) mean methylation values (if NULL, set randomly)}

\item{prop.causal.x}{: proportion of causal cpg M -> x}

\item{prop.causal.y}{: proportion of causal cpg M -> y}

\item{prop.causal.ylx}{: proportion of causal y in causal x}

\item{prop.variance.y}{: proportion of phenotypic variance explained by latent structure (intensity of confounding)}

\item{prop.variance.x}{: proportion of exposure variance explained by latent structure (intensity of confounding)}

\item{rho}{: correlation outcome/exposure (direct effect)}

\item{sigma}{:    standard deviation of residual errors}

\item{sd.A}{:    standard deviation for effect sizes (A: M->X)}

\item{mean.A}{:    (vector) mean of effect sizes}

\item{sd.B}{: standard deviation for effect sizes (B: M->Y)}

\item{mean.B}{:    (vector) mean of effect sizes}

\item{sd.U}{: (vector) standard deviations for factors}

\item{sd.V}{: standard deviations for loadings}

\item{strength}{: if you use real data for loading of cell type, strength of the real data}

\item{sd.ct}{: standard deviations for loadings (cell type)}

\item{alpha}{: parameter for the dirichlet distribution (for cell type), default : runif(K.ct)}

\item{prob.bin}{: if you use binairy exposure (X), probability of success on each trial}

\item{prop.causal}{: proportion of causal variables (probes/loci)}
}
\value{
M : matrix of methylation beta values

X : exposure

Y : phenotype/health outcome

A : effect sizes exposure

B : effect sizes phenotype/health outcome

M.bin : matrix of methylation beta values, use if you use the binairy exposure

X.bin : binairy exposure

CT : proportion of cell type

CT.l : loading of cell type

mediators : set of true mediators

causal.x : set of CpGs associated with the exposure

causal.y : set of CpGs associated with the outcome

U : simulated confounders

V : loadings of coufounders

freq : mean methylation values

controls : true control gene (NOT USE for simulation study)
}
\description{
r_mediation_real_cell_type : function to simulate DNA methylation data for mediation analyzis (and cell type)
}
\details{
This function is used to simulate datasets for analysis of mediations.
The simulation model is based on linear relationships.
First, it construct a covariance matrix for X, Y and U using the parameter rho
(direct effect or correlation between X and Y) and propvar
(intensity of the confounders or correlation between Y and U).
Then this matrix is used to simulate via normal laws X, Y and U.
Thereafter, the effect sizes of X (A), Y (B) and U (V) are calculated
using mean parameters of effect sizes (meanA and meanB) and standard deviations (sdA, sdB and sdV).
Note that the effect sizes of X and Y are calculated only for causal mediators with X and/or Y.
For non-causal mediators, the effect sizes is 0.
On the other hand, a residual error matrix is calculated via the sigma (Z) parameter.
For Cell type (CT) and loading of cell type (CT.l) you can use real data.
If Cell type is simulate :
Cell type is simulate with dirichlet distribution, loading of cell type simulate via normal laws.
To finish the methylation matrix is calculated thanks to the formula :
M = VU + CT.l*CT + AX + BY + Z
}
\examples{
# Simulate data :
simu <- r_mediation_cell_type(100, 500, 2, 5)
}
