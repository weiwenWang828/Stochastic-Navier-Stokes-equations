# An efficient fully explicit scheme for stochastic Navier-Stokes equation driven by multiplicative noise

This repository contains MATLAB implementations of numerical methods for stochastic Navier-Stokes equations as described in our paper.

## Code Organization

### Section 3.3: OAV vs TAV Method Comparison
- **Location**: `Matlab files4SDE/` folder
- **Description**: Code for comparing the OAV (Operator-splitting Auxiliary Variable) and TAV (Time-splitting Auxiliary Variable) methods

### Section 5.1: Accuracy Tests  
- **Location**: `Accuracy test/` folder
- **Dependencies**: Some subroutines are adapted from the spectral methods described in:
  - J. Shen, T. Tang, L.-L. Wang. *Spectral Methods: Algorithms, Analysis and Applications*[M]. Berlin: Springer, 2010.

### Section 5.1: Stochastic Shear Layer Roll-up Problem
- **Location**: `Stochastic shear layer roll-up problem/` folder
- **Description**: Implementation for the stochastic shear layer roll-up test case

## Requirements
- MATLAB (version R2018a or later recommended)
- Basic MATLAB toolboxes

## Usage
Navigate to the respective folder and run the main script files. Please refer to the comments within each file for specific usage instructions.

## Citation
If you use this code in your research, please cite our paper and the referenced spectral methods book.
