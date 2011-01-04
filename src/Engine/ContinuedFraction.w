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
#include "ApplyOperatorLocal.h"
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
		template<typename,typename> class LanczosSolverTemplate,
    	typename ModelType_,
	 	typename ConcurrencyType_,
    	typename IoType_,
       	template<typename> class VectorType>
	class ContinuedFraction  {
	public:
		@<typedefs@>
		@<constructor@>
		@<publicfunctions@>
	
	private:
		@<privatefunctions@>
		@<privatedata@>
	};
@}

A long series of typedefs follow. Need to explain these maybe (FIXME).
@d typedefs
@{
typedef ModelType_ ModelType;
typedef ConcurrencyType_ ConcurrencyType;
typedef IoType_ IoType;
typedef typename ModelType::RealType RealType;
typedef std::complex<RealType> ComplexType;
typedef typename ModelType::OperatorsType OperatorsType;
typedef std::vector<ComplexType> ComplexVectorType;
typedef LanczosSolverTemplate<InternalProductType,ComplexVectorType> LanczosSolverType;
typedef std::vector<RealType> VectorType;
typedef psimag::Matrix<ComplexType> ComplexMatrixType;
typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
typedef ComplexVectorType TargetVectorType;

static const size_t parallelRank_ = 0; // ContF needs to support concurrency FIXME
@}

Now let us look at the private data of this class:
@d privatedata
@{
const ModelType& model_;
ProgressIndicator progress_;
SparseMatrixType hamiltonian_;
RealType gsEnergy_;
VectorType gsVector_;
@}
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
	printHeader();
	// task 1: Compute Hamiltonian and
	// task 2: Compute ground state |phi>
	computeGroundState();
	
	// task 3: compute |initVector> =\sum_x c_x|phi>, where 
	// c_x are some operator
	VectorType initVector;
	computeInitVector(initVector);
	MatrixType T;
	
	// task 4: tridiag H starting with |initVector>
	tridiagonalize(T,initVector);
	
	// task 5: diag. T and store the result
	// to be able to produce GreenFuction(i,j)
	// analytically
	diagonalizeAndStore(T,initVector);
}
@}

Note the tasks listed above. This will be performed under private functions below.

This is the Green function for site i, and site j, at time t.
It should be complex, FIXME:
@d publicfunctions
@{
@<gsEnergy@>
@<greenFunction@>
@}

@d gsEnergy
@{
RealType gsEnergy() const
{
	return gsEnergy_;
}
@}

@d greenFunction
@{
// nothing here yet
@}
	
The private functions are as follows:

@d privatefunctions
@{
@<computeHamiltonian@>
@<computeGroundState@>
@<computeInitVector@>
@<triDiagonalize@>
@<diagonalizeAndStore@>
@}

@d computeHamiltonian
@{
void computeHamiltonian(SparseMatrixType& hamiltonian)
{
	BasisType basis1(nsite,npart1,nhilbert1);
	BasisType basis2(nsite,npart2,nhilbert2);
	
	model_.setupHamiltonian(hamiltonian,basis1,basis2);
}
@}

@d computeGroundState
@{
void computeGroundState(const SparseMatrixType& hamiltonian)
{
	computeHamiltonian(hamiltonian_);
				
	RealType eps= 0.01*ProgramGlobals::LanczosTolerance;
	size_t iter= ProgramGlobals::LanczosSteps;
	size_t parallelRank = 0;

	LanczosSolverType lanczosSolver(hamiltonian,iter,eps,parallelRank);

	lanczosSolver.computeGroundState(gsEnergy_,gsVector_);
}
@}

@d computeInitVector
@{
void computeInitVector(VectorType& initVector)
{
	VectorType tmpVector;
	size_t i = 6;
	size_t j = 6;
	size_t spin = ModelType::SPIN_UP;
	size_t destruction = ModelType::DESTRUCTION;
	SparseMatrixType ci = model_.getOperator(destruction,i,spin);
	SparseMatrixType cj = model_.getOperator(destruction,j,spin);
	ci.multiply(tmpVector,gsVector);
	cj.multiply(initVector,gsVector);
	initVector += tmpVector;
}
@}

@d triDiagonalize
@{
void triDiagonalize(MatrixType& T,const VectorType& initVector)
{
	// tridiagonalize starting with tmpVector = c^\dagger_i|gsVector>
	TridiagonalMatrixType ab;
	MatrixType V;
	lanczosSolver.tridiagonalDecomposition(initVector,ab,V);
	ab.buildDenseMatrix(T);
	//return lanczosSolver.steps();
}
@}

@d diagonalizeAndStore
@{
void diagonalizeAndStore(MatrixType& T,const VectorType& initVector)
{
	std::vector<>RealType> eig(T.n_row());
	utils::diag(T,eigs,'V');
	RealType norma = norm2(initVector);
}
@}

@o ContinuedFraction.h -t
@{
} // namespace Dmrg

#endif
@}
\end{document}
