# Partial correlation networks of Gaussian processes

## Abstract:

In Gaussian graphical models, conditional independence and partial correlations are the natural inferential target for understanding direct relationships in multivariate data. No comparable framework exists for spatial processes, where multivariate analysis defaults to modeling unconditional cross-covariance structure, even though direct relationships remain a scientific quantity of interest in many applied settings.
We address this gap by establishing a characterization of process-level partial correlation for multivariate Gaussian processes that recovers a direct link with Gaussian graphical models. Our analysis proceeds through a class of stationary multivariate processes, termed spectrally inside-out, in which a precision matrix modulates the strength of conditional dependence and yields necessary and sufficient conditions for conditional independence.
Within this class, partial cross-correlation functions factorize into a process-level partial correlation coefficient and an attenuation term independent of cross-process parameters.
This is, to our knowledge, the first partial correlation analysis available at the process level for multivariate Gaussian processes.
The spectrally inside-out class includes the separable coregionalization model, a process convolution construction, and the parsimonious multivariate Mat\'ern, for which such a characterization was previously thought unavailable. We further show that a nonstationary inside-out model satisfies the same factorization and admits the same necessary and sufficient conditions.
Our results clarify the limitations of existing approaches: linear coregionalization models encode conditional independence through the zero pattern of the inverse factor loading matrix and do not result in interpretable partial cross-correlation functions. Low-rank spatial factor models lack a meaningful graphical characterization. Methods that enforce network structure through auxiliary graphical layers only characterize presence or absence of graph edges. We illustrate our results through synthetic and real data.

### Preprint: [https://arxiv.org](https://arxiv.org)


## Requirements:


Code structure: 

-  The `.R` files in the main folders give helper functions for fitting and postprocessing
-  Reproduce Figure 1 via `figure1/Figure_1_unconditional_vs_conditional_v2.r`
-  The illustrations in the paper are split into 3 subsections, each corresponding to a folder. In the order in which they appear: `glasso`, `sim_pmatern`, `jura`.

R packages:

-  meshed: `devtools::install_github("mkln/meshed", ref="ff2c9c0876ca4ef28e08c565982b24f7c5bb64e3")`
-  GpGpm and associated script. Refer to [https://github.com/yf297/GpGp_multi_paper](https://github.com/yf297/GpGp_multi_paper)
-  spiox: `devtools::install_github("mkln/spiox", ref="dcc69264cde394687aff4aebbb726d2ff2ea9085")`
-  `glasso` 1.11, `BDgraph` 2.74.

