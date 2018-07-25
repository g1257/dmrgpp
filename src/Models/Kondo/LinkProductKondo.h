#ifndef LINKPRODUCTKONDO_H
#define LINKPRODUCTKONDO_H
#include "LinkProductBase.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename ModelHelperType, typename GeometryType>
class LinkProductKondo : public LinkProductBase<ModelHelperType, GeometryType> {

	typedef LinkProductBase<ModelHelperType, GeometryType> BaseType;
	typedef typename GeometryType::AdditionalDataType AdditionalDataType;
	typedef typename BaseType::RealType RealType;

public:

	typedef typename BaseType::PairSizeType PairSizeType;

	enum TermEnum {TERM_HOPPING, TERM_SS, TERM_DENSITY};

	// List of function LinkProduct*.h of each model MUST implement

	// You MUST return the number of Hamiltonian terms your model has
	// term are connections from one site to other site or sites
	SizeType terms() const { return 3; } // see TermEnum above for each of the 3

	// For term given in first argument, return how many dofs this term has
	// dofs are sub-terms inside a term
	// You MUST ignore the last argument unless your model has a
	// site dependent Hilbert space (SDHS)
	SizeType dofs(SizeType term, const AdditionalDataType&) const
	{
		switch (term) {
		case TERM_HOPPING:
			return 2; // up-up and down-down
			break;
		case TERM_SS:
			return 2; // spsm and szsz
			break;
		case TERM_DENSITY:
			return 1; // n n
		default:
			break;
		}

		throw PsimagLite::RuntimeError("dofs(): Invalid term\n");
	}

	// Fill 3 outputs and 3 more if you want to support SU(2)
	// Assumes you are connecting two operators A_i x B_j
	// FERMION or BOSON: Say FERMION if both operators A and B are FERMIONS, BOSON otherwise
	// PairSizeType: give indices of A and B in the Model one-site operator storage
	//  std::pair<char,char>: char must be either N or C to indicate whether
	//                         A or B needs transpose conjugate
	// You don't need to fill AngularMomentum AngularFactor Category
	// but then your model will not support SU(2) symmetry
	void setLinkData(SizeType term, // term (INPUT)
	                 SizeType dof, // dof for term (INPUT)
	                 bool isSu2, // isSU2 (INPUT)
	                 ProgramGlobals::FermionOrBosonEnum& fOrB, // FERMION or BOSON (OUTPUT)
	                 PairSizeType& ops, // Pair of operator indices (OUTPUT)
	                 std::pair<char,char>&, // Modifier for each operator (OUTPUT)
	                 SizeType&, // AngularMomentum for SU(2) (OUTPUT)
	                 RealType&, // AngularFactor for SU(2) (OUTPUT)
	                 SizeType&, // Category for SU(2) (OUTPUT)
	                 const AdditionalDataType&) const // For SDHS (INPUT)
	{
		assert(!isSu2); // no SU(2) support for this model
		switch (term) {
		case TERM_HOPPING:
			fOrB = ProgramGlobals::FERMION;
			assert(dof < 2);
			ops = PairSizeType(dof, dof);
			break;
		case TERM_SS:
			fOrB = ProgramGlobals::BOSON;
			assert(dof < 2);
			ops = PairSizeType(dof + 2, dof + 2);
			break;
		case TERM_DENSITY:
			fOrB = ProgramGlobals::BOSON;
			assert(dof == 0);
			ops = PairSizeType(4, 4);
			break;
		default:
			err("Invalid term\n");
			break;
		}
	}
};
}
#endif // LINKPRODUCTKONDO_H
