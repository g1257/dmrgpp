\documentclass{report}
\usepackage[T1]{fontenc}
\usepackage{bera}

\usepackage[pdftex,usenames,dvipsnames]{color}
\usepackage{listings}
\definecolor{mycode}{rgb}{0.9,0.9,1}
\lstset{language=c++,tabsize=1,basicstyle=\scriptsize,backgroundcolor=\color{mycode}}

\usepackage{hyperref}
\usepackage{fancyhdr}
\rhead{DensityMatrixLocal.h}
\pagestyle{fancy}
\usepackage{verbatim}
\begin{document}

%\title{The DensityMatrix Class}
%\author{G.A.}
%\maketitle

\section{The DensityMatrixLocal Class}

\begin{comment}
@o DensityMatrixLocal.h -t
@{
/*
@i license.txt
*/
@}
\end{comment}

These are the include guards and included files:

@o DensityMatrixLocal.h -t
@{
#ifndef DENSITY_MATRIX_LOCAL_H
#define DENSITY_MATRIX_LOCAL_H

#include "Utils.h"
#include "BlockMatrix.h"
#include "DensityMatrixBase.h"
@}

This is a templated class, and it is templated on 4 templates.
RealType is the types for reals, like double or float.
DmrgBasisType is a light Hilbert space (no operators, only symmetries), whereas
DmrgBasisWithOperatorsType is a, well, a Hilbert space with operators.
Finally TargettingType has to do with how to target states, {\it i.e.,\/} the functionality
to include this or that state in the density matrix.

@o DensityMatrixLocal.h -t
@{
namespace Dmrg {
	//!
	template<
		typename RealType,
		typename DmrgBasisType,
		typename DmrgBasisWithOperatorsType,
		typename TargettingType
		>
@}

Now come some typedef that make it easy use types without resorting to multi-line names
for these types. In other words, these are aliases to shorten things later.
@o DensityMatrixLocal.h -t
@{
	class DensityMatrixLocal : public DensityMatrixBase<RealType,DmrgBasisType,DmrgBasisWithOperatorsType,TargettingType> {
		typedef typename DmrgBasisWithOperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename TargettingType::VectorWithOffsetType TargetVectorType;
		typedef typename TargettingType::TargetVectorType::value_type DensityMatrixElementType;
		typedef BlockMatrix<DensityMatrixElementType,PsimagLite::Matrix<DensityMatrixElementType> > BlockMatrixType;
		typedef typename DmrgBasisType::FactorsType FactorsType;
		enum {EXPAND_SYSTEM = TargettingType::EXPAND_SYSTEM };
@}

Note the word \verb|public| below. All under public is what this class ``exports'', {\it i.e.\/}, it
is the functionality that this class provides.
@o DensityMatrixLocal.h -t
@{
	public:
		typedef typename BlockMatrixType::BuildingBlockType BuildingBlockType;
@}

OK, we are now ready for the constructor of this class.
The constructor takes 4 arguments. The target object which deals with what states to target, i.e.,
to include in the density matrix. Two \verb|BasisWithOperator|s objects, \verb|pBasis|
and \verb|pBasisSummed|.
Sometimes the system (left block) is a free index and the environment (right block) is summed over.
So, sometimes \verb|pBasis| will be the system and \verb|pBasisSummed| the environment, and other times 
vice-versa.
Then, there's \verb|pSE|, a light Hilbert space object (BasisType), which represents the superblock%'
(system+environment). Remember that superblock objects are always light (i.e. do not
contain operators) due to memory reasons. The argument \verb|direction| indicates if we're expanding%'
the system or expanding the environment instead. Finally  the \verb|verbose| variable tells us
if we want to print informational stuff.

@o DensityMatrixLocal.h -t
@{
		DensityMatrixLocal(
			const TargettingType& target,
			const DmrgBasisWithOperatorsType& pBasis,
			const DmrgBasisWithOperatorsType& pBasisSummed,
			const DmrgBasisType& pSE,
			SizeType direction,bool debug=false,bool verbose=false) 
@}

Note the colon that comes here indicating that we're setting stuff on the stack.
We're constructing one object, \verb|data_|, which is the actual density matrix.
The type of \verb|data_| is a block diagonal matrix.
\verb|data_|'s constructor takes, in turn, two arguments, \verb|pBasis.size()|,%'
and |pBasis.partition()-1|. The first one is the rank (number of rows) of the density matrix.
It coincides with the size of the Hilbert space that is free, {\it i.e.\/} not summed over
in the reduced density matrix formula.
So, if the reduced density formula is
\begin{equation}
\rho_{\alpha,\alpha'} = \sum_{\beta} \psi_{\alpha,\beta}
\psi_{\alpha',\beta},
\label{eq:densitymatrix}
\end{equation}
then the Hilbert space corresponding to $\alpha$ is \verb|pBasis|.
The second argument is the number of blocks for \verb|data_| which is equal to the number
of partitions of basis minus one. (The off-one is due to the partitions having
a terminating partition at the end with the last state of the basis.)

Note that we also initialize the debug and verbose flags, whose meaning are immediate.
@o DensityMatrixLocal.h -t
@{
		: data_(pBasis.size() ,pBasis.partition()-1),
				debug_(debug),verbose_(verbose)
		{
		}
@}

The \verb|operator()| function is a public member function that returns the actual density matrix.
What d'you mean by ``actual density matrix'', you might ask? Well, it's a\ldots 
wait for it\ldots matrix, of
course, you know like $\rho_{\alpha,\alpha'}$.%'
@o DensityMatrixLocal.h -t
@{

		virtual BlockMatrixType& operator()()
		{
			return data_;
		}
@}

The rank function comes below, and returns the rank (as in the number of rows). I don't remember why this
function is needed at all. Need to check where it's used. It seems that having the matrix is
enough
@o DensityMatrixLocal.h -t
@{
		virtual SizeType rank() { return data_.rank(); }
@}

These are just checks to see if everything is OK. They're empty for this local symmetries implementation.%'
@o DensityMatrixLocal.h -t
@{
		virtual void check(int direction)
		{
		}

		virtual void check2(int direction)
		{
		}
@}

The function \verb|diag| diagonalizes the density matrix, which, if you remember, is one
of the key steps of the DMRG algorithm. 
Three arguments are passed here. First \verb|eigs| which will be filled with the eigenvalues
of the density matrix. Next \verb|jobz| which is either 'N' or 'V' indicating if we need
also eigenvectors ('V') or only eigenvalues ('N'). 
Finally a \verb|concurrency| object is also passed that can help with parallelizing this 
diagonalization operation.
We don't do much here, just delegate the work to the \verb|BlockMatrix.h| class
that knows how to diagonalize a block matrix (i.e. block diagonal matrix).%'
@o DensityMatrixLocal.h -t
@{
		template<typename ConcurrencyType>
		void diag(typename PsimagLite::Vector<RealType>::Type& eigs,char jobz,ConcurrencyType& concurrency)
		{
			diagonalise<DensityMatrixElementType,RealType,ConcurrencyType>(data_,eigs,jobz,concurrency);
		}
@}

Below we initialize the density matrix, which basically computes it.
We delay initializing the actual density matrix class, and don't do it in the constructor.%'
This is because we only need one type (either local or SU(2)) per run, 
and initializing both would be a waste of 
resources (both CPU and memory).
The argument this function takes are similar to the constructor of the class and won't be explained again.%'
This function calls \verb|initPartition| to accumulate the density matrix for each target state,
since the density matrix is the sum of the density matrices for each state we want to target.
@o DensityMatrixLocal.h -t
@{		virtual void init(
				const TargettingType& target,
				DmrgBasisWithOperatorsType const &pBasis,
				const DmrgBasisWithOperatorsType& pBasisSummed,
				DmrgBasisType const &pSE,
				int direction)
		{	
			//loop over all partitions:
			for (SizeType m=0;m<pBasis.partition()-1;m++) {
				// size of this partition
				SizeType bs = pBasis.partition(m+1)-pBasis.partition(m);
				
				// density matrix block for this partition:
				BuildingBlockType matrixBlock(bs,bs);
				
				// weight of the ground state:
				RealType w = target.gsWeight();
				
				// if we are to target the ground state do it now:
				if (target.includeGroundStage())
					initPartition(matrixBlock,pBasis,m,target.gs(),pBasisSummed,pSE,direction,w);
				
				// target all other states if any:
				for (SizeType i=0;i<target.size();i++) {
					w = target.weight(i)/target.normSquared(i);
					initPartition(matrixBlock,pBasis,m,target(i),pBasisSummed,pSE,direction,w);
				}
				
				// set this matrix block into data_
				data_.setBlock(m,pBasis.partition(m),matrixBlock);
			}
		}
@}

The following statement (see the semicolon at the end) makes the function \verb|operator<<| a friend
of this class. This function will enable printing this class for debugging purposes.
@o DensityMatrixLocal.h -t
@{
		template<
			typename RealType_,
			typename DmrgBasisType_,
			typename DmrgBasisWithOperatorsType_,
   			typename TargettingType_
			> 
		friend std::ostream& operator<<(std::ostream& os,
				const DensityMatrixLocal<RealType_,
    					DmrgBasisType_,DmrgBasisWithOperatorsType_,TargettingType_>& dm);
@}

This class has 3 data memebers, all of them private. 
We've already seen \verb|data_|, a block diagonal matrix that contains the actual density matrix.%'
@o DensityMatrixLocal.h -t
@{
	private:
		BlockMatrixType data_;
		bool debug_,verbose_;
@}

OK, this function was called above and accumulates the density matrix for a single target vector 
\verb|v| and a single symmetry sector \verb|m|. It also weights the contribution of \verb|v| with \verb|weight|.
The algorithm distingushes if we're expanding the system or the environment.%'
@o DensityMatrixLocal.h -t
@{
		void initPartition(BuildingBlockType& matrixBlock,
				DmrgBasisWithOperatorsType const &pBasis,
				SizeType m,
				const TargetVectorType& v,
				DmrgBasisWithOperatorsType const &pBasisSummed,
				DmrgBasisType const &pSE,
    				SizeType direction,
				RealType weight)
		{
			if (direction!=EXPAND_SYSTEM) 
				initPartitionExpandEnviron(matrixBlock,pBasis,m,v,pBasisSummed,pSE,weight);
			else
				initPartitionExpandSystem(matrixBlock,pBasis,m,v,pBasisSummed,pSE,weight);
		}
@}

So, if we're expanding the environment, we compute $\rho_{i,j}$ but only for matrix block \verb|m|.%'
Note the block matrix offset \verb|pBasis.partition(m)| being sustracted.
We delegate the actual computation for each $i,j$ to the function \verb|densityMatrixExpandEnviron|
described below.
@o DensityMatrixLocal.h -t
@{
		void initPartitionExpandEnviron(BuildingBlockType& matrixBlock,
				DmrgBasisWithOperatorsType const &pBasis,
				SizeType m,
				const TargetVectorType& v,
				DmrgBasisWithOperatorsType const &pBasisSummed,
				DmrgBasisType const &pSE,
				RealType weight)
		{
			
			SizeType ns=pBasisSummed.size();
			SizeType ne=pSE.size()/ns;
			
			for (SizeType i=pBasis.partition(m);i<pBasis.partition(m+1);i++) {
				for (SizeType j=pBasis.partition(m);j<pBasis.partition(m+1);j++) {
						
					matrixBlock(i-pBasis.partition(m),j-pBasis.partition(m)) +=
						densityMatrixExpandEnviron(i,j,v,pBasisSummed,pSE,ns,ne)*weight;
					}
				}
		}
@}

So, if we're expanding the system, we compute $\rho_{i,j}$ but only for matrix block \verb|m|.
Note the block matrix offset \verb|pBasis.partition(m)| being sustracted.
We delegate the actual computation for each $i,j$ to the function \verb|densityMatrixExpandSystem|
described below.
@o DensityMatrixLocal.h -t
@{
		void initPartitionExpandSystem(BuildingBlockType& matrixBlock,
				DmrgBasisWithOperatorsType const &pBasis,
				SizeType m,
				const TargetVectorType& v,
				DmrgBasisWithOperatorsType const &pBasisSummed,
				DmrgBasisType const &pSE,
				RealType weight)
		{
			SizeType ne = pBasisSummed.size();
			SizeType ns = pSE.size()/ne;
			
			for (SizeType i=pBasis.partition(m);i<pBasis.partition(m+1);i++) {
				for (SizeType j=pBasis.partition(m);j<pBasis.partition(m+1);j++) {
						
					matrixBlock(i-pBasis.partition(m),j-pBasis.partition(m)) +=
						densityMatrixExpandSystem(i,j,v,pBasisSummed,pSE,ns,ne)*weight;
				}
			}
		}
@}

So, this function computes the density matrix contribution to $\rho(\alpha_1,\alpha_2)$%'
but only due to symmetry sector \verb|m|, and only due to
target state \verb|v|, and only when we're expanding the environment.%'
This is exactly Eq.~(\ref{eq:densitymatrix}).
@o DensityMatrixLocal.h -t
@{
		DensityMatrixElementType densityMatrixExpandEnviron(
				SizeType alpha1,
				SizeType alpha2,
				const TargetVectorType& v,
				DmrgBasisWithOperatorsType const &pBasisSummed,
				DmrgBasisType const &pSE,
				SizeType ns,
				SizeType ne)
		{
			
			SizeType total=pBasisSummed.size();
			
			DensityMatrixElementType sum=0;
			
			SizeType x2 = alpha2*ns;
			SizeType x1 = alpha1*ns;
			for (SizeType beta=0;beta<total;beta++) {
				SizeType jj = pSE.permutationInverse(beta + x2);
				SizeType ii = pSE.permutationInverse(beta + x1);
				sum += v[ii] * std::conj(v[jj]);
			}
			return sum;
		}
@}
So, this function computes the density matrix contribution to $\rho(\alpha_1,\alpha_2)$
but only due to symmetry sector \verb|m|, and only due to
target state \verb|v|, and only when we're expanding the system.%'
This is exactly Eq.~(\ref{eq:densitymatrix}).
@o DensityMatrixLocal.h -t
@{
		DensityMatrixElementType densityMatrixExpandSystem(
				SizeType alpha1,
				SizeType alpha2,
				const TargetVectorType& v,
				DmrgBasisWithOperatorsType const &pBasisSummed,
				DmrgBasisType const &pSE,
				SizeType ns,
				SizeType ne)
		{
			
			SizeType total=pBasisSummed.size();
			
			DensityMatrixElementType sum=0;

			for (SizeType beta=0;beta<total;beta++) {
				SizeType jj = pSE.permutationInverse(alpha2+beta*ns);
				SizeType ii = pSE.permutationInverse(alpha1+beta*ns);
				sum += v[ii] * std::conj(v[jj]);
			}
			return sum;
		}
	}; // class DensityMatrixLocal
@}

Below is a companion function that prints the density matrix object for debugging purposes.
We loop over all the blocks of the (block diagonal) density matrix contained in \verb|dm.data_|.
Printing of each block is delegated to the type of each block (usually a Psimag matrix).
@o DensityMatrixLocal.h -t
@{
	template<
		typename RealType,
		typename DmrgBasisType,
		typename DmrgBasisWithOperatorsType,
  		typename TargettingType
		> 
	std::ostream& operator<<(std::ostream& os,
				const DensityMatrixLocal<RealType,DmrgBasisType,DmrgBasisWithOperatorsType,TargettingType>& dm)
	{
		for (SizeType m=0;m<dm.data_.blocks();m++) {
			SizeType ne = dm.pBasis_.electrons(dm.pBasis_.partition(m));
			os<<" ne="<<ne<<"\n"; 
			os<<dm.data_(m)<<"\n";
		}
		return os;
	}
} // namespace Dmrg

#endif
@}
\end{document}

 
