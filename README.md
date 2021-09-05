# R-code-for-dual-inverse-regression
The R code is specified for causal inference with the binary exposure, which is the scenario in the numerical studies of the paper. It can be easily adjusted for the general cases. When implemeting the dual inverse regression method using this R code, the readers can choose SIR (Li, 1991), SAVE (Cook and Weisberg, 1991), directional regression (Li and Wang, 2007), or PWSVM (Shin, Wu, Zhang, and Liu, 2017) to construct M_Y and M_T.

The R code includes three main functions:
1. dirm: the function that runs the dual inverse regression methods.
2. ladle.dirm: the function that applies the ladle estimator to determine the dimension of the dual inverse regression subspaces.
3. ladle.glb: the function that determines whether the globally efficient dimension reduction subspace exists (acceptance if the output is 0).

All the other R functions are the prerequisites to run these three functions. The R code for PWSVM is downloaded from https://github.com/liyf1988/PrincipalWeightedSVM/blob/master/PrincipalWeightedSVM.r, which was uploaded by Dr. Liu, one of the authors of Shin, Wu, Zhang, and Liu (2017).

Please see the header of each R function for more details. Please also contact the author at luoluowei@gmail.com or wluostat@yahoo.com for any questions or requests.

References:
Cook, R. D., Weisberg, S., 1991. Discussion of "Sliced inverse regression for dimension reduction". Journal of the American Statistical Association 86, 316-342.
Li, K. C., 1991. Sliced inverse regression for dimension reduction (with discussion). Journal of the American Statistical Association 86, 316-342.
Li, B., Wang, S., 2007. On directional regression for dimension reduction. Journal of the American Statistical Association 102 (479), 997-1008.
Shin, S. J., Wu, Y., Zhang, H., Liu, Y., 2017. Principal weighted support vector machines for sufficient dimension reduction in binary classification. Biometrika 104, 67-81.
