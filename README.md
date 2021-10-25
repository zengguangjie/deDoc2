# deDoc2
deDoc2 is a TAD-like domain(TLD) prediction tool using structural information theory, it treats the Hi-C contact map as a weighted graph, 
and applys dynamic programming algorith to globally optimize the two-dimensional structural entropy of the graph partiton. 
The deDoc2 package consists of two predictors, deDoc2.w and deDoc2.s, to predict higher level, larger and lower level, smaller TLDs, respectively. 
The deDoc2.w minimizes the structural entropy in the whole Hi-C contact map, while the deDoc2.s minimizes the structural entropy in the matrices of sliding windows along the genome. 
We previously developed deDoc for bulk Hi-C TAD predicting at https://github.com/yinxc/structural-information-minimisation.
