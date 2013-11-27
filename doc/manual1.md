Input File Format
=================================

\section inputfile Input File

There is a single input file that is passed as the first
and only argument to the program.
There are three kinds of parameters in the input file:
(i) model connections (``geometry'') parameters, (ii) model on-site parameters, and (iii) DMRG Solver parameters.
The Model parameters vary from model to model.
The DMRG Solver parameters are discussed below.

@copydoc Dmrg::ParametersDmrgSolver

\section finiteloops Finite Loops

@copydoc Dmrg::FiniteLoop