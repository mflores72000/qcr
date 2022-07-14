# qcr: Quality Control Review

## Version 1.3

This package implements univariate and multivariate Statistical Quality Control (SQC) tools that completes and increases the SQC techniques available in R.
Apart from integrating different R packages devoted to SQC (`qcc`,`MSQC`), provides nonparametric tools that are highly useful when Gaussian assumption is not met. 
This package computes standard univariate control charts for individual measurements, X-bar, S, R, p, np, c, u, EWMA and CUSUM.


This package also includes functions to perform multivariate control charts such as Hotelling T2, MEWMA and MCUSUM. 
As representative feature, multivariate nonparametric alternatives based on data depth are implemented in this package: r, Q and S control charts. 
In addition, Phase I and II control charts for functional data are included. 

This package also allows the estimation of the most complete set of capability indices from first to fourth generation, covering the nonparametric alternatives, and performing the corresponding capability analysis graphical outputs, including the process capability plots.


## References

Flores, M., Naya, S., Fernández-Casal, R., Zaragoza, S., Raña, P., and Tarrío-Saavedra, J. (2020). Constructing a control chart using functional data. *Mathematics*, **8**, 58,
[DOI](https://doi.org/10.3390/math8010058).

Flores, M., Fernández-Casal, R., Naya, S., & Tarrío-Saavedra, J. (2021). Statistical Quality Control with the qcr Package. *R Journal*, **13**, 194-217. [DOI](http://doi.org/10.32614/rj-2021-034).

Naya, S., Devia-Rivera, A., Saavedra, J. T., & Flores, M. A. (2016). Nueva propuesta de índices de capacidad robustos para el control de la calidad. Dyna, 83(198), 94-101. [DOI](https://doi.org/10.15446/dyna.v83n198.49930).