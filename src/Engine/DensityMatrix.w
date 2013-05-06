\documentclass{report}
\usepackage[T1]{fontenc}
\usepackage{bera}

\usepackage[pdftex,usenames,dvipsnames]{color}
\usepackage{listings}
\definecolor{mycode}{rgb}{0.9,0.9,1}
\lstset{language=c++,tabsize=1,basicstyle=\scriptsize,backgroundcolor=\color{mycode}}

\usepackage{hyperref}
\usepackage{fancyhdr}
\rhead{DensityMatrix.h}
\pagestyle{fancy}
\usepackage{verbatim}
\begin{document}

%\title{The DensityMatrix Class}
%\author{G.A.}
%\maketitle

\section{The DensityMatrix Class}
\begin{comment}
@o DensityMatrix.h -t
@{
/*
@i license.txt
*/
@}
\end{comment}


 This class contains the implementation of 
the density matrix calculation  and diagonalization for the DMRG algorithm.
This is a lazy class, it doesn't do much but instead delegates the work%'
to either the |DensityMatrixLocal| for when there's only local symmetries,
or to the |DensityMatrixSu2| when there's local symmetries and SU(2) symmetry.

The file starts with the normal define guards:
@o DensityMatrix.h -t
@{
#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H
@}

Next come the include files, we need to include \verb|Utils.h|
since it is need for things like multiplications of vectors, and sorting. 
Most of this stuff is gonna end up in PsimagLite eventually, once I get that in place.
@o DensityMatrix.h -t
@{
#include "Utils.h"
@}

Anyway, we need also to include \verb|BlockMatrix.h| that provides capabilities for a 
block matrix. A block matrix is, well, a matrix composed of blocks, a diagonal matrix
as you might call it. Why in hell do you need a block matrix you might (rightfully) ask?
Well, you see, there are symmetries that block not only the Hamiltonian but also the 
density matrix, which (probably, I think, I need to check) has the same symmetries of the 
Hamiltonian
@o DensityMatrix.h -t
@{
#include "BlockMatrix.h"
@}

Now this class works in conjunction with three more files \verb|DensityMatrixBase.h|,
\verb|DensityMatrixLocal.h|, and \verb|DensityMatrixSu2.h|, as explained in detail below.
@o DensityMatrix.h -t
@{
#include "DensityMatrixLocal.h"
#include "DensityMatrixSu2.h"

namespace Dmrg {
@}

This is a templated class, and it is templated on 4 templates.
RealType is the types for reals, like double or float.
DmrgBasisType is a light Hilbert space (no operators, only symmetries), whereas
DmrgBasisWithOperatorsType is a, well, a Hilbert space with operators.
Finally TargettingType has to do with how to target states, {\it i.e.,\/} the functionality
to include this or that state in the density matrix.
@o DensityMatrix.h -t
@{
	template<
		typename RealType,
		typename DmrgBasisType,
		typename DmrgBasisWithOperatorsType,
		typename TargettingType
		>
	class DensityMatrix {
@}

Now come some typedef that make it easy use types without resorting to multi-line names
for these types. In other words, these are aliases to shorten things later.
@o DensityMatrix.h -t
@{
		enum {EXPAND_SYSTEM = TargettingType::EXPAND_SYSTEM };

		typedef typename DmrgBasisWithOperatorsType::SparseMatrixType
			SparseMatrixType;
		typedef typename TargettingType::TargetVectorType::value_type
			DensityMatrixElementType;
		typedef BlockMatrix<DensityMatrixElementType,
			PsimagLite::Matrix<DensityMatrixElementType> > BlockMatrixType;
		typedef typename DmrgBasisType::FactorsType FactorsType;
		typedef DensityMatrixLocal<RealType,DmrgBasisType,
			DmrgBasisWithOperatorsType, TargettingType>
			DensityMatrixLocalType;
		typedef DensityMatrixSu2<RealType,DmrgBasisType,
			DmrgBasisWithOperatorsType,TargettingType>
			DensityMatrixSu2Type;
		typedef DensityMatrixBase<RealType,DmrgBasisType,
			DmrgBasisWithOperatorsType,TargettingType>
			DensityMatrixBaseType;
@}

Note the word \verb|public| below. All public stuff is the functionality that
this class provides.
@o DensityMatrix.h -t
@{
	public:
		typedef typename BlockMatrixType::BuildingBlockType
			BuildingBlockType;
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

@o DensityMatrix.h -t
@{
		DensityMatrix(
			const TargettingType& target,
			const DmrgBasisWithOperatorsType& pBasis,
			const DmrgBasisWithOperatorsType& pBasisSummed,
			const DmrgBasisType& pSE,
			size_t direction,
			bool debug=false,
			bool verbose=false)
			
@}

Note the colon that comes here indicating that we're setting stuff on the stack.
We're constructing two objects, one to handle local symmetries only, and another one to
handle local and SU(2) symmetries. These objects have very light constructors that are
described in (not sure how to cross reference with literate programming, need to learn more!!).
@o DensityMatrix.h -t
@{
			: densityMatrixLocal_(target,pBasis,pBasisSummed,pSE,
					direction,debug,verbose),
				densityMatrixSu2_(target,pBasis,pBasisSummed,pSE,
					direction,debug,verbose)
		{
@}

OK so here we need to see if we are using only local symmetries or we're also using the %'
SU(2) symmetry. As you can imagine, things are different in each case. 
Both \verb|DensityMatrixLocal| and \verb|DensityMatrixSu2| derive from a common parent class.
The trick is to select which one we need and have the \verb|densityMatrixImpl_| pointer point to the
correct one.
@o DensityMatrix.h -t
@{
			if (DmrgBasisType::useSu2Symmetry()) {
				densityMatrixImpl_ = &densityMatrixSu2_;
			} else {
				densityMatrixImpl_ = &densityMatrixLocal_;
			}
@}
Below we initialize the density matrix, which basically computes it.
We delay initializing the actual density matrix class until here.
This is because we only need one type (either local or SU(2)) per run, 
and initializing both would be a waste of 
resources (both CPU and memory).
Note, however, that created both \verb|densityMatrixLocal_| and \verb|densityMatrixSu2_|; but their creation
is light as explained before.
@o DensityMatrix.h -t
@{
			densityMatrixImpl_->init(target,pBasis,pBasisSummed,pSE,direction);
		}
@}

The \verb|operator()| function is a public member function that returns the actual density matrix.
What d'you mean by ``actual density matrix'', you might ask? Well, it's a\ldots 
wait for it\ldots matrix, of
course, you know like $\rho_{ij}$. So, how do we implement that? Well, we can't, simply%'
because remember that we need to consider two cases, local and local plus SU(2) symmetries.
So, the best we can do here is delegate it to the appropriate class, either 
\verb|DensityMatrixLocal| or \verb|DensityMatrixSu2|. Remember, however, that we already set the pointer
\verb|densityMatrixImpl_| to the right object, so we simply call its \verb|operator()| function.
@o DensityMatrix.h -t
@{
		BlockMatrixType& operator()()
		{
			return densityMatrixImpl_->operator()();
		}
@}
The rank function comes below, and returns the rank (as in the number of rows). I don't remember why this
function is needed at all. Need to check where it's used. It seems that having the matrix is
enough.
Again nothing can be actually done here, we simply delegate to the appropriate class.
@o DensityMatrix.h -t
@{
		size_t rank() { return densityMatrixImpl_->rank(); }
@}

These are just checks to see if everything is OK.
@o DensityMatrix.h -t
@{
		void check(int direction)
		{
			return densityMatrixImpl_->check(direction);
		}
@}

And more checks to see if everything is OK.
@o DensityMatrix.h -t
@{
		void check2(int direction)
		{
			densityMatrixImpl_->check2(direction);
		}
@}
The function \verb|diag| diagonalizes the density matrix, which, if you remember, is one
of the key steps of the DMRG algorithm. 
Three arguments are passed here. First \verb|eigs| which will be filled with the eigenvalues
of the density matrix. Next \verb|jobz| which is either 'N' or 'V' indicating if we need
also eigenvectors ('V') or only eigenvalues ('N'). 
Finally a \verb|concurrency| object is also passed that can help with parallelizing this 
diagonalization operation.
Again, we don't do much here, just delegate, and
it's some other class's problem.%'
@o DensityMatrix.h -t
@{
		template<typename ConcurrencyType>
		void diag(typename PsimagLite::Vector<RealType>::Type& eigs,char jobz,
				ConcurrencyType& concurrency)
		{
			if (!DmrgBasisType::useSu2Symmetry()) {
				densityMatrixLocal_.diag(eigs,jobz,concurrency);
			} else {
				densityMatrixSu2_.diag(eigs,jobz,concurrency);
			}
		}
@}
The following statement (see the semicolon at the end) makes the function \verb|operator<<| a friend
of this class. This function will enable printing this class for debugging purposes.
@o DensityMatrix.h -t
@{
		template<
			typename RealType_,
			typename DmrgBasisType_,
			typename DmrgBasisWithOperatorsType_,
   			typename TargettingType_
			> 
		friend std::ostream& operator<<(std::ostream& os,
			const DensityMatrix<RealType_,DmrgBasisType_,
				DmrgBasisWithOperatorsType_,TargettingType_>&
						dm);
@}

This class has 3 data memeber, all of them private. 
We've already seen \verb|densityMatrixLocal_| that does the real work
for the density matrix when there's local symmetries, and  also
\verb|densityMatrixSu2_| that does the real work when there's local symmetries and the SU(2) symmetry.%'
As we explained above, \verb|densityMatrixImpl_| is a pointer that points to the 
correct object, depending on which symmetries the user chose.
@o DensityMatrix.h -t
@{
	private:
		DensityMatrixLocalType densityMatrixLocal_;
		DensityMatrixSu2Type densityMatrixSu2_;
		DensityMatrixBaseType* densityMatrixImpl_;
	}; // class DensityMatrix
@}

Below is a companion function that prints the density matrix object for debugging purposes.
Note that we just print the pointer \verb|densityMatrixImpl_|, and so, again we're delegating%'
all the work to the appropriate class(es).
@o DensityMatrix.h -t
@{
	template<
		typename RealType,
		typename DmrgBasisType,
		typename DmrgBasisWithOperatorsType,
  		typename TargettingType
		> 
	std::ostream& operator<<(std::ostream& os,
				const DensityMatrix<RealType,DmrgBasisType,
				DmrgBasisWithOperatorsType,TargettingType>& dm)
	{
		os<<(*dm.densityMatrixImpl_);
		return os;
	}
} // namespace Dmrg

#endif

@}
And that's all folks. 
\end{document}


