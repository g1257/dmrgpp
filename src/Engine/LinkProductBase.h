#ifndef LINKPRODUCTBASE_H
#define LINKPRODUCTBASE_H
#include "Vector.h"
#include "ProgramGlobals.h"
#include "PsimagLite.h"

namespace Dmrg {

template<typename ModelHelperType_, typename GeometryType_>
class LinkProductBase {

public:

	typedef ModelHelperType_ ModelHelperType;
	typedef GeometryType_ GeometryType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename GeometryType::AdditionalDataType AdditionalDataType;
	typedef typename ModelHelperType::RealType RealType;
	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	template<typename SomeInputType>
	LinkProductBase(SomeInputType& io, PsimagLite::String terms)
	{
		PsimagLite::split(termNames_, terms);
		SizeType n = termNames_.size();

		for (SizeType i = 0; i < n; ++i) {
			PsimagLite::String label = "Term" + ttos(i);
			PsimagLite::String str("");

			try {
				io.readline(str, label + "=");
			} catch(std::exception&) {
				std::cerr<<"You should add a " + label + "=" + termNames_[i];
				std::cerr<<" for this term to your input\n";
				continue;
			}

			if (str == termNames_[i]) continue;
			PsimagLite::String msg = label + " expected to be " + termNames_[i];
			err(msg + " but " + str + " found instead\n");
		}
	}

	virtual ~LinkProductBase() {}

	// List of function LinkProduct*.h of each model MUST implement

	// You MUST return the number of Hamiltonian terms your model has
	// term are connections from one site to other site or sites
	virtual SizeType terms() const = 0;

	// For term given in first argument, return how many dofs this term has
	// dofs are sub-terms inside a term
	// You MUST ignore the last argument unless your model has a
	// site dependent Hilbert space (SDHS)
	virtual SizeType dofs(SizeType, const AdditionalDataType&) const = 0;

	// Fill 3 outputs and 3 more if you want to support SU(2)
	// Assumes you are connecting two operators A_i x B_j
	// FERMION or BOSON: Say FERMION if both operators A and B are FERMIONS, BOSON otherwise
	// PairSizeType: give indices of A and B in the Model one-site operator storage
	//  std::pair<char,char>: char must be either N or C to indicate whether
	//                         A or B needs transpose conjugate
	// You don't need to fill AngularMomentum AngularFactor Category
	// but then your model will not support SU(2) symmetry
	virtual void setLinkData(SizeType, // term (INPUT)
	                         SizeType, // dof for term (INPUT)
	                         bool, // isSU2 (INPUT)
	                         ProgramGlobals::FermionOrBosonEnum&, // FERMION or BOSON (OUTPUT)
	                         PairSizeType&, // Pair of operator indices (OUTPUT)
	                         std::pair<char,char>&, // Modifier for each operator (OUTPUT)
	                         SizeType&, // AngularMomentum for SU(2) (OUTPUT)
	                         RealType&, // AngularFactor for SU(2) (OUTPUT)
	                         SizeType&, // Category for SU(2) (OUTPUT)
	                         const AdditionalDataType&) const = 0; // For SDHS (INPUT)

	// List of function LinkProduct*.h of each model MAY override

	// Assumes you are connecting two operators A_{i, alpha} x B_{j, beta}
	// You MUST set edofs[0] = alpha and edofs[1] = beta
	virtual void connectorDofs(VectorSizeType& edofs, // (OUTPUT)
	                           SizeType, // TERM (INPUT)
	                           SizeType, // DOF (INPUT)
	                           const AdditionalDataType&) const // For SDHS (INPUT)
	{
		edofs[0] = edofs[1] = 0;
	}

	// Assumes your are connecting two operators value * A_{i, alpha} x B_{j, beta}
	// You MAY add an extra multiplier value, by setting the first argument
	virtual void valueModifier(ComplexOrRealType&, // value (OUTPUT)
	                           SizeType, // TERM (INPUT)
	                           SizeType, // DOF (INPUT)
	                           bool, // isSU2 (INPUT)
	                           const AdditionalDataType&) const // For SDHS (INPUT)
	{}

	// You MUSTN'T override this function for now
	virtual SizeType dofsAllocationSize() const { return 2; }

protected:

	const VectorStringType& termNames() const { return termNames_; }

private:

	VectorStringType termNames_;
};
}
#endif // LINKPRODUCTBASE_H
