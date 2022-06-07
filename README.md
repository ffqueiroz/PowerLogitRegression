# PowerLogitRegression
###### R codes and datasets for the paper "Power logit regression for modeling bounded data".

This repository contains R codes and datasets used in the applications of the paper "Power logit regression for modeling bounded data" by Queiroz and Ferrari (2022). 

To reproduce the codes, it is necessary to install some packages, including the PLreg package, available at https://github.com/ffqueiroz/PLreg.
To install the PLreg package, you can type in the R console the following commands:

```
install.packages("devtools")
devtools::install_github("ffqueiroz/PLreg")
```
If the devtools package does not work, install the PLreg package using the remotes package, with the following commands:
```
install.packages("remotes")
remotes::install_github("ffqueiroz/PLreg", ref = "main")
```
