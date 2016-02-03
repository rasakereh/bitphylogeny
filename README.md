---
## Overview

This README file discusses the code associated with the manuscript:

Ke Yuan, Thomas Sakoparnig, Florian Markowetz and Niko Beerenwinkel. (2015) 
BitPhylogeny: a probabilistic framework for reconstructing intra-tumor phylogenies. Genome Biology. Accepted

The Python code is based on the implementation from 
Ryan Prescott Adams, Zoubin Ghahramani and Michael I. Jordan. (2010)
Tree-Structured Stick Breaking for Hierarchical Data. NIPS.
http://hips.seas.harvard.edu/files/tssb.tgz
---
## License

BitPhylogeny is licensed under the GPL v3, see the LICENSE.txt file for details.

---
## Dependencies

1. Python:
   * NumPy 1.9.2,
   * SciPy 0.15.1,
   * scikit-learn 0.16.1: http://scikit-learn.org/stable/ 
   * rpy2 2.5.6: http://rpy.sourceforge.net/ 
   * pandas 0.16.1: http://pandas.pydata.org/
   * h5py 2.5.0: http://www.h5py.org/
   * python-igraph 0.7: http://igraph.org/python/

2. R:
   rPython, mcclust, e1071, igraph, gplots, riverplot, plyr.
   ```
   sessionInfo()
   R version 3.2.0 (2015-04-16)
   Platform: x86_64-apple-darwin13.4.0 (64-bit)
   Running under: OS X 10.10.5 (Yosemite)

   locale:
   [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

  attached base packages:
   [1] stats     graphics  grDevices utils     datasets  methods   base     

  other attached packages:
   [1] bitphylogenyR_0.99 igraph_0.7.1       rPython_0.0-6      RJSONIO_1.3-0     

  loaded via a namespace (and not attached):
   [1] riverplot_0.5      Rcpp_0.11.6        gtools_3.5.0       class_7.3-12      
   [5] mcclust_1.0        bitops_1.0-6       plyr_1.8.2         e1071_1.6-7       
   [9] KernSmooth_2.23-14 gplots_2.17.0      gdata_2.17.0       tools_3.2.0       
  [13] cluster_2.0.1      caTools_1.17.1 
  ```

---
## Installation

Clone the BitPhylogeny repository:
```
git clone git@bitbucket.org:ke_yuan/bitphylogeny.git 
cd bitphylogeny
```

Install the BitPhylogeny python package:
```
cd python
sudo python setup.py install
```

Install the R package:
```
cd ../R
R CMD INSTALL bitphylogenyR_0.99.tar.gz
```
---
## Useage
We've provided a sample data file 'sample_data.csv' in the python folder. To run the sample file, execute the following code. 
```
BitPhylogeny analyse 'sample_data.csv' 'output' -true_label -n 200 -b 10 -t 5 -mode "methylation" -seed 1234 
```
```
'sample_data.csv': the data file. Please make sure rows represent data points, colums represent dimension 
'output': the output folder where the MCMC traces and tree scripts are stored.
-true_label: data file contains true labels as the last column.
-n: number of mcmc samples in total 
-b: number of burn-in samples
-t: thin number (a current requirement is (number of mcmc samples / thin) > 5)
-mode: "methylation" or "mutation"
-seed: default 1234
-row_names: the first colum of the data file is row names 
```
The output folder contains the mcmc-traces and treescripts. 

For post processing, please use the accompanying R package. Tutorials of the R package are stored in doc/ directory