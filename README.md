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
## Dependencies and tested versions

1. Python 2.7.6:
      * NumPy 1.9.2,
      * SciPy 0.15.1,
      * scikit-learn 0.16.1: http://scikit-learn.org/stable/ 
      * rpy2 2.5.6: http://rpy.sourceforge.net/ 
      * pandas 0.16.1: http://pandas.pydata.org/
      * h5py 2.5.0: http://www.h5py.org/
      * python-igraph 0.7.1: http://igraph.org/python/

2. R 3.2.0:
      * rPython 0.0-6
      * mcclust 1.0
      * e1071 1.6-7 
      * igraph 0.7.1 
      * gplots 2.17.0 
      * riverplot 0.5 
      * plyr 1.8.2.

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
-collect_all_trees: option to store all trees (default: false)
```
The output folder contains the mcmc-traces, treescripts, alltreescripts (if -collect_all_trees is flagged). 

For post processing, please use the accompanying R package. Tutorials of the R package are stored in doc/ directory
For instruction on the Bitphylogeny graph, please find doc/bitphylogeny_node.html