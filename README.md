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
   NumPy,SciPy
   scikit-learn: http://scikit-learn.org/stable/ 
   rpy2: http://rpy.sourceforge.net/ 
   pandas: http://pandas.pydata.org/
   h5py: http://www.h5py.org/
   igraph

2. R:
   rPython, mcclust, e1071, igraph, gplots, riverplot, plyr.

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
R CMD INSTALL bitphylogenyR_0.1.tar.gz
```
---
## Useage
We've provided a sample data file 'sample_data.csv' in the python folder. To run the sample file, execute the following code. 
```
BitPhylogeny analyse 'sample_data.csv' 'output' -true_label -n 200 -b 10 -t 5 -mode "methylation" -rand_seed 1234 
```
```
'sample_data.csv': the data file
'output': the output folder where the MCMC traces and tree scripts are stored.
-true_label: data file contains true labels as the last column.
-n: number of mcmc samples in total 
-b: number of burn-in samples
-t: thin number (a current requirement is (number of mcmc samples / thin) > 5)
-mode: "methylation" or "mutation"
-rand_seed: default 1234
```
The output folder contains the mcmc-traces and treescripts. 

For post processing, please use the accompanying R package. Tutorials of the R package are stored in doc/ directory