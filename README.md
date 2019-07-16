# SPC

## Statistical Parametric Connectome (SPC)

This is the Python implementation of the Statistical Parametric Connectome (SPC) a.k.a Two-step procedure for statistical testing of connectomes introduced by D.-E. Meskaldji.

The Matlab and R versions can be found here: https://miplab.epfl.ch/index.php/software/spc

**Reference:** D.-E. Meskaldji et al, "Improved statistical evaluation of group differences in connectomes by screening-filtering strategy with application to study maturation of brain connections between childhood and adolescence", NeuroImage, 2015 (https://www.sciencedirect.com/science/article/pii/S1053811914009896?via%3Dihub)


The main function is called `main_relaxPvalues.py` and imports the `padjust.py` and `computeK2.py`.

**Dependencies:** `numpy`, `scipy`


## Toy example

```
from padjust import padjust
from computeK2 import computeK2
from main_relaxPvalues import relaxPvalues
    
SubregionsIndices = [1,1,1,2,2,2] # a list 
p_values = [0.001, 0.0005, 0.0001, 0.2, 0.4, 0.2] # a list
assert(len(SubregionsIndices) == len(p_values))
Alpha = 0.05
globalAlpha = 0.05
method = 'fdr'
q = 1 # verbose

subsetPvalues,relaxedpvalues,subsetScores = relaxPvalues(p_values,SubregionsIndices,Alpha,globalAlpha,method,q)
```
