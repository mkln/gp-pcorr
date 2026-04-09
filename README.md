# Requirements:


R packages:

-  meshed: `devtools::install_github("mkln/meshed", ref="ff2c9c0876ca4ef28e08c565982b24f7c5bb64e3")`
-  GpGpm and associated script. Refer to [https://github.com/yf297/GpGp_multi_paper](https://github.com/yf297/GpGp_multi_paper)
-  spiox: `devtools::install_github("mkln/spiox", ref="dcc69264cde394687aff4aebbb726d2ff2ea9085")`
-  `glasso` 1.11, `BDgraph` 2.74.


Structure: 

-  The `.R` files in the main folders give helper functions for fitting and postprocessing
-  Reproduce Figure 1 via `figure1/Figure_1_unconditional_vs_conditional_v2.r`
-  The illustrations in the paper are split into 3 subsections, each corresponding to a folder. In the order in which they appear: `glasso`, `sim_pmatern`, `jura`.
