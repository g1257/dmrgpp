# Input 


## Introduction

There is a single input file that is passed as the argument to `-f`, like so

    ./dmrg -f input.inp.

Examples of input files can be found under `TestSuite/inputs/`
There are three kinds of parameters in the input file:
(i) model connections (``geometry'') parameters, (ii) model on-site parameters, and (iii) DMRG Solver parameters.
Each type of input parameters is discussed below.

## Geometry Input
@copydoc hide_geometry1
@copydoc hide_geometry2

## Model Input

The Model parameters vary from model to model.

## DMRG Solver parameters

@copydoc hide_ParametersDmrgSolver

## Finite Loops

@copydoc hide_FiniteLoop

