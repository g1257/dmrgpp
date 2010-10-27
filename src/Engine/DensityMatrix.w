\nocon % omit table of contents
\datethis % print date on listing

@* The |DensityMatrix| Class. This class contains the implemenation of 
the density matrix calculation  and diagonalization for the DMRG algorithm.
This is a lazy class, it doesn't do much but instead delegates the work
to either the |DensityMatrixLocal| for when there's only local symmetries,
or to the |DensityMatrixSu2| when there's local symmetries and SU(2) symmetry.

@ The file starts with the normal define guards:
@c
#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H

@ Next come the include files, we need to include \.{Utils.h}
since it is need for things like multiplications of vectors, and sorting. 
Most of this stuff is gonna end up in PsimagLite eventually, once I get that in place.
@c
#include "Utils.h"

@ Anyway, we need also to include \.{BlockMatrix.h} that provides capabilities for a 
block matrix. A block matrix is, well, a matrix composed of blocks, a diagonal matrix
as you might call it. Why in hell do you need a block matrix you might (rightfully) ask?
Well, you see, there are symmetries that block not only the Hamiltonian but also the 
density matrix, which (probably, I think, I need to check) has the same symmetries of the 
Hamiltonian
@c
#include "BlockMatrix.h"

@ Now this class works in conjunction with three more files \.{DensityMatrixBase.h},
\.{DensityMatrixLocal.h}, and \.{DensityMatrixSu2.h}, as explained in detail below.
@c
#include "DensityMatrixLocal.h"
#include "DensityMatrixSu2.h"

namespace Dmrg {

@ This is a templated class, and it is templated on 4 templates.
RealType is the types for reals, like double or float.
DmrgBasisType is a light Hilbert space (no operators, only symmetries), whereas
DmrgBasisWithOperatorsType is a, well, a Hilbert space with operators.
Finally TargettingType has to do with how to target states, {\it i.e.,\/} the functionality
to include this or that state in the density matrix.
@c
	template<
		typename RealType,
		typename DmrgBasisType,
		typename DmrgBasisWithOperatorsType,
		typename TargettingType
		>
	class DensityMatrix {

@ Now come some typedef that make it easy use types without resorting to multi-line names
for these types. In other words, these are aliases to shorten things later.
@c
		typedef typename DmrgBasisWithOperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename TargettingType::TargetVectorType::value_type DensityMatrixElementType;
		typedef BlockMatrix<DensityMatrixElementType,psimag::Matrix<DensityMatrixElementType> > 
			BlockMatrixType;
		typedef typename DmrgBasisType::FactorsType FactorsType;
		enum {EXPAND_SYSTEM = TargettingType::EXPAND_SYSTEM };
		@/
		typedef DensityMatrixLocal<RealType,DmrgBasisType,DmrgBasisWithOperatorsType,
			TargettingType> DensityMatrixLocalType;
		typedef DensityMatrixSu2<RealType,DmrgBasisType,DmrgBasisWithOperatorsType,
			TargettingType> DensityMatrixSu2Type;
		typedef DensityMatrixBase<RealType,DmrgBasisType,DmrgBasisWithOperatorsType,
			TargettingType> DensityMatrixBaseType;

@ Note the word |public| below. All public stuff is the functionality that
this class provides.
@c
	public:
		typedef typename BlockMatrixType::BuildingBlockType BuildingBlockType;

@ OK, we are now ready for the constructor of this class.
The constructor takes 4 arguments. The target object which deals with what states to target, i.e.,
to include in the density matrix. Two BasisWithOperators objects, pBasis and pBasisSummed.
Sometimes the system (left block) is a free index and the environment (right block) is summed over.
So, sometimes |pBasis| will be the system and |pBasisSummed| the environment, and other times 
vice-versa.
Then, there's pSE, a light Hilbert space object (BasisType), which represents the superblock 
(system+environment). Remember that superblock objects are always light (i.e. do not
contain operators) due to memory reasons. The argument |direction| indicates if we're expanding
the system or expanding the environment instead. Finally  the |verbose| variable tells us
if we want to print informational stuff.

@c
		DensityMatrix(
			const TargettingType& target,
			const DmrgBasisWithOperatorsType& pBasis,
			const DmrgBasisWithOperatorsType& pBasisSummed,
			const DmrgBasisType& pSE,
			size_t direction,
			bool debug=false,
			bool verbose=false)
			
@
Note the semicolon that comes here indicating that we're setting stuff on the stack.
We're constructing two objects, one to handle local symmetries only, and another one to
handle local and SU(2) symmetries. These objects have very light constructors that are
described in (not sure how to cross reference with literate programming, need to learne more!!).
@c
			: densityMatrixLocal_(target,pBasis,pBasisSummed,pSE,direction,debug,verbose),
				densityMatrixSu2_(target,pBasis,pBasisSummed,pSE,direction,debug,verbose)
		{

@ OK so here we need to see if we are using only local symmetries or we're also using the 
SU(2) symmetry. As you can imagine, things are different in each case. 
Both |DensityMatrixLocal| and |DensityMatrixSu2| derive from a common parent class.
The trick is to select which one we need and have the |densityMatrixImpl_| pointer point to the
correct one.
@c
			if (DmrgBasisType::useSu2Symmetry()) {
				densityMatrixImpl_ = &densityMatrixSu2_;
			} else {
				densityMatrixImpl_ = &densityMatrixLocal_;
			}
@ Note that we delay initializing the actual density matrix class we need.
This is because we only need one type per run, and initializing both would be a waste of 
resources (both CPU and memory).
However, we need to create both |densityMatrixLocal_| and |densityMatrixSu2_|; their creation
is light as explained before.
@c
			densityMatrixImpl_->init(target,pBasis,pBasisSummed,pSE,direction);
		}

@ The operator() function is a public member function that returns the actual density matrix.
What you mean by ``actual density matrix'', you might ask? Well, it's a ... 
wait for it ... matrix, of
course, you know like $\rho_{ij}$. So, how do we implement that? Well, we can't, simply
because remember that we need to consider two cases, local and local plus SU(2) symmetries.
So, the best we can do here is delegate it to the appropriate class, either 
|DensityMatrixLocal| or |DensityMatrixSu2|. Remember, however, that we already set the pointer
|densityMatrixImpl_| to the right object, so we simply call its operator() function.
@c
		BlockMatrixType& operator()()
		{
			return densityMatrixImpl_->operator()();
		}

@ The rank function returns the rank (as in the number of rows). I don't remember why this
function is needed at all. Need to check where it's used. It seems that having the matrix is
enough.
Again nothing can be actually done here, we can simply delegate to the appropriate class.
@c
		size_t rank() { return densityMatrixImpl_->rank(); }
		
@ Ah, checks to see if everything is OK.
@c
		void check(int direction)
		{
			return densityMatrixImpl_->check(direction);
		}
		
@ And more checks to see if everything is OK.
@c
		void check2(int direction)
		{
			densityMatrixImpl_->check2(direction);
		}
		
@ The function |diag| diagonalizes the density matrix, which, if you remember, is one
of the key steps of the DMRG algorithm. Again we don't do much here, just delegate, and
it's some other class's problem.
@c
		template<typename ConcurrencyType>
		void diag(std::vector<RealType>& eigs,char jobz,ConcurrencyType& concurrency)
		{
			if (!DmrgBasisType::useSu2Symmetry()) {
				densityMatrixLocal_.diag(eigs,jobz,concurrency);
			} else {
				densityMatrixSu2_.diag(eigs,jobz,concurrency);
			}
			
		}

@ The following statement (see the semicolon at the end) makes the function |operator<<| a friend
of this class. This is just to be able to print this class for debugging purposes
@c
		template<
			typename RealType_,
			typename DmrgBasisType_,
			typename DmrgBasisWithOperatorsType_,
   			typename TargettingType_
			> 
		friend std::ostream& operator<<(std::ostream& os,
				const DensityMatrix<RealType_,
    					DmrgBasisType_,DmrgBasisWithOperatorsType_,TargettingType_>& dm);

@ All data is private. We've already seen the |densityMatrixLocal_| that does the real work
for the density matrix when there's local symmetries, and 
|densityMatrixSu2_| that does the real work when there's local symmetries and the SU(2) symmetry.
|densityMatrixImpl_| is a pointer that points to the correct object depending on which
symmetries the user chose.
@c
	private:
		DensityMatrixLocalType densityMatrixLocal_;
		DensityMatrixSu2Type densityMatrixSu2_;
		DensityMatrixBaseType* densityMatrixImpl_;

		
	}; // class DensityMatrix

@ This is a companion function that prints the density matrix object for debugging purposes.
Note that we just print the pointer |densityMatrixImpl_|, and so, again we're delegating
all the work to the appropriate class(es).
@c
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

@ And that's all folks. An index might follow here. 
 
