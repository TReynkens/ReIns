<!-- README.md is generated from README.Rmd. Please edit that file -->
ReIns
=====

ReIns contains functions and datasets of the book "Reinsurance: Actuarial and Statistical Aspects" (2015) of Albrecher, Beirlant and Teugels. It contains implementations of

-   Basic extreme value theory (EVT) estimators and graphical methods as described in "Statistics of Extremes: Theory and Applications" (2004) of Beirlant, Goegebeur, Segers and Teugels.

-   EVT estimators and graphical methods adapted for censored and/or truncated data.

-   Splicing of mixed Erlang distributions with EVT distributions (Pareto, GPD).

-   VaR, expected shortfall and excess-loss premium estimates.

All code is written in plain R but optimised by applying vectorisation when possible.

The package is not yet available on CRAN but you can install the latest development version from GitHub using

    install.packages("devtools")

    devtools::install_github("TReynkens/ReIns")
