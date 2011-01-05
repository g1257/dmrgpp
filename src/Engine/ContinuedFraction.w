\documentclass{report}
\usepackage[T1]{fontenc}
\usepackage{bera}

\usepackage[pdftex,usenames,dvipsnames]{color}
\usepackage{listings}
\definecolor{mycode}{rgb}{0.9,0.9,1}
\lstset{language=c++,tabsize=1,basicstyle=\scriptsize,backgroundcolor=\color{mycode}}

\usepackage{hyperref}
\usepackage{fancyhdr}
\pagestyle{fancy}
\usepackage{verbatim}
\begin{document}

%\title{The ContinuedFraction Class}
%\author{G.A.}
%\maketitle

\begin{comment}
@o ContinuedFraction.h -t
@{
/*
@i license.txt
*/
@}
\end{comment}

\section{Continued Fractions for Dynamic Observables}
This class implements the dynamic observable computation described in
Rev. Mod. Phys. 66, 763 (1993), pages 776ff. 
Usually, this algorithm is not at all accurate for the DMRG since
there are so few states due to the truncation.
However, one use of this algorithm is to check the time-step targetting 
(implemented in TimeStepTargetting.w) with a lanczos-only approach.
In this case, one runs the DMRG without truncation (select m large enough).

@o ContinuedFraction.h -t
@{
#ifndef CONTINUED_FRACTION_H 
#define CONTINUED_FRACTION_H
#include <iostream>
#include "ProgressIndicator.h"
#include "BLAS.h"
#include "LanczosSolver.h"
@}

This class is templated on 7 templates, which are:
\begin{enumerate}
\item \verb|LanczosSolverTemplate|, being usually the \verb|LanczosSolver| class.
\item \verb|ModelType| is the model in question. These are classes under the directory Models.
\item \verb|ConcurrenyType| is the type to deal with parallelization or lack thereof.
\item \verb|IoType| is usually the \verb|IoSimple| class, and deals with writing to disk the time-vectors produced by this class.
\item \verb|VectorType| is usually the \verb|std::vector| class that encapsulates
the functionality of a vector.
\end{enumerate}

@o ContinuedFraction.h -t 
@{
namespace Dmrg {
	template<
    	typename ModelType_,
	 	typename ConcurrencyType_>
	class ContinuedFraction  {
	public:
		@<typedefs@>
		@<constructor@>
		@<publicfunctions@>
	
	private:
		@<privatefunctions@>
		@<privatedata@>
	}; // class ContinuedFraction
} // namespace Dmrg

#endif 
@}

A long series of typedefs follow. Need to explain these maybe (FIXME).
@d typedefs 
@{
typedef ModelType_ ModelType;
typedef ConcurrencyType_ ConcurrencyType;
typedef typename ModelType::RealType RealType;
typedef typename ModelType::VectorType VectorType;
typedef typename ModelType::SparseMatrixType SparseMatrixType;
typedef typename VectorType::value_type FieldType;
typedef typename ModelType::BasisType BasisType;

typedef LanczosSolver<RealType,SparseMatrixType,VectorType> LanczosSolverType;
typedef psimag::Matrix<FieldType> MatrixType;
typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;

static const size_t parallelRank_ = 0; // ContF needs to support concurrency FIXME @}

Now let us look at the private data of this class:
@d privatedata 
@{
const ModelType& model_;
ProgressIndicator progress_;
SparseMatrixType hamiltonian_;
RealType gsEnergy_;
VectorType gsVector_; @}
which contains a \verb|model| reference, a progress indicator object, and storage for
the ground state energy and ground state vector.

Now we get to the constructor and stack initialization of this object (note the colon).
The constructor takes as one and only argument the model object. 
We make a reference copy for the model.
We initialize the \verb|progress| object that helps with printing progress to the terminal.

@d constructor 
@{
ContinuedFraction(const ModelType& model)
	: model_(model),progress_("ContinuedFraction",0)
{
	// printHeader();
	// task 1: Compute Hamiltonian and
	// task 2: Compute ground state |phi>
	computeGroundState();
	
} @}

Note the tasks listed above. This will be performed under private functions below.

This is the Green function for site i, and site j, at time t.
It should be complex, FIXME:
@d publicfunctions
@{
@<gsEnergy@>
@<greenFunction@> @}

@d gsEnergy
@{
RealType gsEnergy() const
{
	return gsEnergy_;
} @}

@d greenFunction
@{
void getGreenFunction(TridiagonalMatrixType& ab,RealType& norma,
		size_t i,size_t j) const
{
	// task 3: compute |initVector> =\sum_x c_x|phi>, where
	// c_x are some operator
	VectorType initVector;
	computeInitVector(initVector,i,j);

	// task 4: tridiag H starting with |initVector>
	triDiagonalize(ab,initVector);

	norma = initVector*initVector;
}
@}
	
The private functions are as follows:

@d privatefunctions 
@{
@<computeGroundState@>
@<computeInitVector@>
@<triDiagonalize@>
@}


@d computeGroundState 
@{
void computeGroundState()
{
	model_.setupHamiltonian(hamiltonian_);
	MatrixType fm;
	crsMatrixToFullMatrix(fm,hamiltonian_);
	std::cerr<<fm;
	if (!isHermitian(fm)) throw std::runtime_error("Hamiltonian non Hermitian\n");
	//std::cerr<<hamiltonian_;
	std::cerr<<"Done setting up Hamiltonian\n";

	RealType eps= 0.01*ProgramGlobals::LanczosTolerance;
	size_t iter= ProgramGlobals::LanczosSteps;
	size_t parallelRank = 0;

	LanczosSolverType lanczosSolver(hamiltonian_,iter,eps,parallelRank);
	gsVector_.resize(hamiltonian_.rank());
	lanczosSolver.computeGroundState(gsEnergy_,gsVector_);
} @}

@d computeInitVector
@{
void computeInitVector(VectorType& initVector,size_t i,size_t j) const
{
	initVector.resize(model_.size());
	VectorType tmpVector(initVector.size());
	size_t spin = ModelType::SPIN_UP;
	size_t destruction = ModelType::DESTRUCTOR;
	SparseMatrixType ci;
	model_.getOperator(ci,destruction,i,spin);
	SparseMatrixType cj;
	model_.getOperator(cj,destruction,j,spin);
	ci.matrixVectorProduct(tmpVector,gsVector_);
	cj.matrixVectorProduct(initVector,gsVector_);
	initVector += tmpVector;
} @}

@d triDiagonalize 
@{
void triDiagonalize(TridiagonalMatrixType& ab,const VectorType& initVector) const
{
	// tridiagonalize starting with tmpVector = c^\dagger_i|gsVector>
	MatrixType V;

	RealType eps= 0.01*ProgramGlobals::LanczosTolerance;
	size_t iter= ProgramGlobals::LanczosSteps;
	size_t parallelRank = 0;

	LanczosSolverType lanczosSolver(hamiltonian_,iter,eps,parallelRank);

	lanczosSolver.tridiagonalDecomposition(initVector,ab,V);
	//ab.buildDenseMatrix(T);
	//return lanczosSolver.steps();
} @}



\end{document}
