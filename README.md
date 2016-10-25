<!-- README.md is generated from README.Rmd. Please edit that file -->
ReIns
=====

The *ReIns* package contains functions from the book "Reinsurance: Actuarial and Statistical Aspects" (2017) of Hansjörg Albrecher, Jan Beirlant and Jef Teugels. It contains implementations of

-   Basic extreme value theory (EVT) estimators and graphical methods as described in "Statistics of Extremes: Theory and Applications" (2004) of Jan Beirlant, Yuri Goegebeur, Johan Segers and Jef Teugels.

-   EVT estimators and graphical methods adapted for censored and/or truncated data.

-   Splicing of mixed Erlang distributions with EVT distributions (Pareto, GPD).

-   VaR, expected shortfall and excess-loss premium estimates.

The package is not yet available on CRAN but you can install the latest development version from GitHub. If you work on Windows, make sure first that \[Rtools\]{<https://cran.r-project.org/bin/windows/Rtools/>} is installed. Then, install the latest development version of *ReIns* using

    install.packages("devtools")

    devtools::install_github("TReynkens/ReIns")
