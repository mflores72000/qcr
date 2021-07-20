# qcr 1.4

* Added a website for the package (with pkgdown).
 
* Added 'README.md', 'NEWS.md' and '_pkgdown.yml'.
 
* Added some vignettes (pkgdown articles): 


# qcr 1.3

* Avoided dependence on `qualityTools` package as it was removed from CRAN 
  (dependencies were included as internal functions in 'qualityTools.R').


# qcr 1.2

* Added utilities for Functional Data Quality Control: `fdqcd()`, `fdqcs.depth()`,     
  `fdqcs.rank()`, and plot methods `plot.fdqcd()`, `plot.fdqcs.depth()`, 
  `plot.fdqcs.rank()`
	

# qcr 1.0

* Added utilities for nonparametric Quality Control: `npqcd()`, `npqcs()`,
  `npqcs.add()`, `npqcs.Q()`, `npqcs.r()`, `npqcs.S()`, `npstate.control()`, and
  plot methods `plot.npqcs()`, `plot.npqcs.Q()`, `plot.npqcs.r()`, `plot.npqcs.S()`.
  
* Added utilities for parametric and nonparametric Capability Analysis: `qcs.ca()`,
  `qcs.cp()`, `qcs.cpn()`, `qcs.hat.cpm()`, `qcs.pcr()`.	


# qcr 0.1-18 

* Initial version in package form.