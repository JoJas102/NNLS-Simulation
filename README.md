# NNLS-simulation

## Description

Program to simulate influence of different parameters on NNLS fitting of diffusion data and compare multi-exponential fitting methods. Evaluating the reliability of finding the total number of components contributing to the corresponding multi-exponantial signal and analysing the results by calculating corresponding diffusion parameters, comparing those to the ground truth.\
This fitting routine uses the regularized NNLS algorithm with cross validation from Thorarin Bjarnason for comparison.

## Initial variables mandatory
Essential simulation and acquisition parameters can be edited and adjusted in [InitVar.m](InitVar.m):
* numberOfIter = number of simulations/iteration steps, differentiation for single- or multi-plotting of results
* snrIn        = Input SNR for artificial DWI signal rawSignal()
* b[]          = array of b-values
* Dmin/Dmax    = defining the D-range for diffusion constants/basis values DValues[] and DBasis[]
* m            = number of diffusion coefficients/bins between Dmin and Dmax\
* eventually:
  * f[] and d[]  = arrays with volume fractions [tissue tubules blood]/[slow, inter and fast] and associated  diffusion coefficients for ground truth (combined   in fdInput[]) or
  * dRange[] and fRange[] = reasonable intervals for the random choice of f and d made by DiffParamRandomizer.m

## General information
MATLAB code will be added soon
