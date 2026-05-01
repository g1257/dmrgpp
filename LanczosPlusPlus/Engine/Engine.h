
/*
// BEGIN LICENSE BLOCK
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[Lanczos++, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/

#ifndef ENGINE_H_
#define ENGINE_H_
#include "BLAS.h"
#include "ContinuedFractionCollection.h"
#include "DefaultSymmetry.h"
#include "GetBraOrKet.h"
#include "LabeledOperator.h"
#include "LanczosGlobals.h"
#include "LanczosSolver.h"
#include "OneOperatorSpec.h"
#include "ParametersForSolver.h"
#include "ProgressIndicator.h"
#include "Random48.h"
#include "TypeToString.h"
#include <iostream>

namespace LanczosPlusPlus {
template <typename ModelType_,
          template <typename, typename> class InternalProductTemplate,
          typename SpecialSymmetryType>
class Engine {

public:

	enum
	{
		SPIN_UP   = LanczosGlobals::SPIN_UP,
		SPIN_DOWN = LanczosGlobals::SPIN_DOWN
	};

	typedef ModelType_                                                   ModelType;
	typedef typename ModelType::InputType                                InputType;
	typedef typename ModelType::SparseMatrixType                         SparseMatrixType;
	typedef typename ModelType::BasisBaseType                            BasisType;
	typedef InternalProductTemplate<ModelType, SpecialSymmetryType>      InternalProductType;
	typedef DefaultSymmetry<typename ModelType::GeometryType, BasisType> DefaultSymmetryType;
	typedef InternalProductTemplate<ModelType, DefaultSymmetryType> InternalProductDefaultType;
	typedef typename SpecialSymmetryType::GeometryType              GeometryType;
	typedef typename GeometryType::ComplexOrRealType                ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type      RealType;
	typedef PsimagLite::Random48<RealType>                          RandomType;
	typedef PsimagLite::ParametersForSolver<RealType>               ParametersForSolverType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type    VectorType;
	typedef typename PsimagLite::Vector<VectorType>::Type           VectorVectorType;
	typedef PsimagLite::LanczosSolver<InternalProductType>          LanczosSolverType;
	typedef PsimagLite::LanczosSolver<InternalProductDefaultType>   LanczosSolverDefaultType;
	typedef PsimagLite::Matrix<ComplexOrRealType>                   MatrixType;
	typedef PsimagLite::Matrix<RealType>                            MatrixRealType;
	typedef typename PsimagLite::Vector<SizeType>::Type             VectorSizeType;
	typedef std::pair<SizeType, SizeType>                           PairType;
	typedef typename PsimagLite::Vector<RealType>::Type             VectorRealType;
	typedef typename PsimagLite::Vector<PsimagLite::String>::Type   VectorStringType;
	typedef LabeledOperator                                         LabeledOperatorType;
	typedef PsimagLite::OneOperatorSpec                             OneOperatorSpecType;
	typedef typename ModelType::RahulOperatorType                   RahulOperatorType;
	typedef typename ModelType::VectorRahulOperatorType             VectorRahulOperatorType;
	using ContFractionType           = PsimagLite::ContinuedFraction<RealType>;
	using CollectionContFractionType = PsimagLite::ContinuedFractionCollection<RealType>;

	// ContF needs to support concurrency FIXME
	static const SizeType parallelRank_   = 0;
	static const SizeType CHECK_HERMICITY = 1;

	enum
	{
		PLUS,
		MINUS
	};

	Engine(const ModelType& model, InputType& io)
	    : model_(model)
	    , progress_("Engine")
	    , io_(io)
	    , options_("")
	{
		io_.readline(options_, "SolverOptions=");
		SizeType excited = 0;

		try {
			io_.readline(excited, "Excited=");
		} catch (std::exception&) { }

		computeAllStatesBelow(excited);
	}

	RealType energies(SizeType ind) const
	{
		assert(ind < energies_.size());
		return energies_[ind];
	}

	const VectorType& eigenvector(SizeType ind) const
	{
		assert(ind < vectors_.size());
		return vectors_[ind];
	}

	//! Calc Green function G(isite,jsite)  (still diagonal in spin)
	void spectralFunction(CollectionContFractionType&               cfCollection,
	                      VectorStringType&                         vstr,
	                      const LabeledOperatorType&                lOperator,
	                      int                                       isite,
	                      int                                       jsite,
	                      const PsimagLite::Vector<PairType>::Type& spins,
	                      const PairType&                           orbs) const
	{
		std::cout << "orbitals=" << orbs.first << " " << orbs.second << "\n";
		for (SizeType i = 0; i < spins.size(); i++) {
			std::cout << "spins=" << spins[i].first << " " << spins[i].second << "\n";
			spectralFunction(
			    cfCollection, vstr, lOperator, isite, jsite, spins[i], orbs);
		}
	}

	/* PSIDOC SpectralFunctions
	Here we document the spectral functions and Green function G(isite,jsite)
	(still diagonal in spin)
	*/
	void spectralFunction(CollectionContFractionType& cfCollection,
	                      VectorStringType&           vstr,
	                      const LabeledOperatorType&  lOperator1,
	                      int                         isite,
	                      int                         jsite,
	                      const PairType&             spins,
	                      const PairType&             orbs) const
	{
		if (spins.first != spins.second) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "spectralFunction: no support yet for off-diagonal spin\n";
			throw std::runtime_error(str.c_str());
		}

		assert(0 < vectors_.size());
		const VectorType&         gsVector   = vectors_[0];
		const LabeledOperatorType lOperator2 = lOperator1.transposeConjugate();

		const BasisType* basisNew   = 0;
		bool             isDiagonal = (isite == jsite && orbs.first == orbs.second);
		PairType         oldParts   = model_.basis().parts();
		for (SizeType type = 0; type < lOperator1.numberOfTypes(); ++type) {
			if (isDiagonal && type > 1)
				continue;

			const LabeledOperatorType& lOperator = (type & 1) ? lOperator1 : lOperator2;

			if (lOperator.needsNewBasis()) {
				assert(spins.first == spins.second);
				std::pair<SizeType, SizeType> newParts(0, 0);
				if (!model_.hasNewParts(
				        newParts, oldParts, lOperator, spins.first, orbs.first))
					continue;
				// Create new bases
				basisNew = model_.createBasis(newParts.first, newParts.second);
			} else {
				basisNew = &model_.basis();
			}
			VectorType modifVector;
			getModifiedState(modifVector,
			                 lOperator,
			                 gsVector,
			                 *basisNew,
			                 type,
			                 isite,
			                 jsite,
			                 spins.first,
			                 orbs);

			SpecialSymmetryType symm(*basisNew, model_.geometry(), "");
			InternalProductType matrix(model_, *basisNew, symm);
			ContFractionType    cf(cfCollection.freqType());

			if (PsimagLite::norm(modifVector) < 1e-10) {
				std::cerr << "spectralFunction: modifVector==0, type=" << type
				          << "\n";
			}

			calcSpectral(cf,
			             lOperator.isFermionic(),
			             modifVector,
			             matrix,
			             type,
			             spins.first,
			             isDiagonal);
			PsimagLite::String str = ttos(spins.first) + "," + ttos(type) + ",";
			str += ttos(orbs.first) + "," + ttos(orbs.second);
			vstr.push_back(str);
			cfCollection.push(cf);
		}
	}

	void measure(const VectorStringType& braOpKet) const
	{
		if (braOpKet.size() != 3)
			err("LanczosDriver1: Only dressed brakets allowed (FATAL ERROR)\n");

		const PsimagLite::String meas = braOpKet[1];
		VectorStringType         tokens;
		PsimagLite::split(tokens, meas, ";");
		const SizeType                                  n = tokens.size();
		VectorRahulOperatorType                         vops;
		VectorSizeType                                  vsites(n);
		typedef typename OneOperatorSpecType::SiteSplit SiteSplitType;

		for (SizeType i = 0; i < n; ++i) {
			SiteSplitType siteSplit = OneOperatorSpecType::extractSiteIfAny(tokens[i]);
			if (!siteSplit.hasSiteString)
				err("Operator " + tokens[i] + " needs a site in brackets\n");
			tokens[i] = siteSplit.root;
			assert(i < vsites.size());
			vsites[i] = OneOperatorSpecType::strToNumberOrFail(siteSplit.siteString);

			OneOperatorSpecType opspec(tokens[i]);

			vops.push_back(
			    RahulOperatorType(opspec.label, opspec.dof, opspec.transpose));
		}

		PsimagLite::GetBraOrKet myket("|" + braOpKet[2]);
		SizeType                ketIndex = myket.levelIndex();
		checkBraOrKet(braOpKet[2], ketIndex);
		const VectorType& ketVector = vectors_[ketIndex];
		VectorType        psiNew(ketVector.size());

		model_.rahulMethod(psiNew, vops, vsites, ketVector, model_.basis());

		PsimagLite::GetBraOrKet mybra(braOpKet[0] + "|");
		SizeType                braIndex = mybra.levelIndex();
		checkBraOrKet(braOpKet[0], braIndex);
		const VectorType& braVector = vectors_[braIndex];

		const ComplexOrRealType result = braVector * psiNew;
		std::cout << braOpKet[0] << "|" << meas << "|" << braOpKet[2] << " = " << result
		          << "\n";
	}

	void twoPoint(PsimagLite::Matrix<typename VectorType::value_type>& result,
	              const LabeledOperatorType&                           lOperator,
	              const PsimagLite::Vector<PairType>::Type&            spins,
	              const PairType&                                      orbs,
	              const PairType&                                      braAndKet) const
	{
		for (SizeType i = 0; i < spins.size(); i++) {
			std::cout << "spins=" << spins[i].first << " " << spins[i].second << "\n";
			twoPoint(result, lOperator, spins[i], orbs, braAndKet);
		}
	}

	/* PSIDOC TwoPointCorrelations
	Here we document the two-point correlations
	*/
	void twoPoint(PsimagLite::Matrix<typename VectorType::value_type>& result,
	              const LabeledOperatorType&                           lOperator,
	              const PairType&                                      spins,
	              const PairType&                                      orbs,
	              const PairType&                                      braAndKet) const
	{
		const BasisType* basisNew = 0;
		PairType         oldParts = model_.basis().parts();

		if (lOperator.needsNewBasis()) {
			if (spins.first != spins.second) {
				PsimagLite::String str(__FILE__);
				str += " " + ttos(__LINE__) + "\n";
				str += "twoPoint: no support yet for off-diagonal spin ";
				str += "when needs new basis\n";
				throw std::runtime_error(str.c_str());
			}

			std::pair<SizeType, SizeType> newParts(0, 0);
			if (!model_.hasNewParts(
			        newParts, oldParts, lOperator, spins.first, orbs.first))
				return;

			basisNew = model_.createBasis(newParts.first, newParts.second);

			std::cerr << "basisNew.size=" << basisNew->size() << " ";
			std::cerr << "newparts.first=" << newParts.first << " ";
			std::cerr << "newparts.second=" << newParts.second << "\n";
		} else {
			basisNew = &model_.basis();
		}

		checkBraOrKet("bra", braAndKet.first);
		checkBraOrKet("ket", braAndKet.second);
		const VectorType& braVector = vectors_[braAndKet.first];
		const VectorType& ketVector = vectors_[braAndKet.second];
		SizeType          total     = result.n_row();

		for (SizeType isite = 0; isite < total; isite++)
			for (SizeType jsite = 0; jsite < total; jsite++)
				result(isite, jsite) = -100;

		RealType isign = 1.0;

		typename VectorType::value_type sum = 0;
		std::cout << "orbs=" << orbs.first << " " << orbs.second << "\n";
		for (SizeType isite = 0; isite < total; isite++) {
			VectorType modifVector1(basisNew->size(), 0);
			if (orbs.first >= model_.orbitals(isite))
				continue;
			accModifiedState(modifVector1,
			                 lOperator,
			                 *basisNew,
			                 ketVector,
			                 isite,
			                 spins.first,
			                 orbs.first,
			                 isign);
			for (SizeType jsite = 0; jsite < total; jsite++) {
				VectorType modifVector2(basisNew->size(), 0);
				if (orbs.second >= model_.orbitals(jsite))
					continue;
				accModifiedState(modifVector2,
				                 lOperator,
				                 *basisNew,
				                 braVector,
				                 jsite,
				                 spins.second,
				                 orbs.second,
				                 isign);
				result(isite, jsite) = modifVector2 * modifVector1;
				if (isite == jsite)
					sum += result(isite, isite);
			}
		}
		std::cout << "MatrixDiagonal = " << sum << "\n";
	}

	// many point, fixed sites
	ComplexOrRealType manyPoint(const VectorSizeType&                                sites,
	                            const PsimagLite::Vector<LabeledOperatorType>::Type& what,
	                            const VectorSizeType&                                spins,
	                            const VectorSizeType&                                orbs,
	                            const PairType& braAndKet) const
	{
		checkBraOrKet("ket", braAndKet.second);
		VectorType       tmpVector = vectors_[braAndKet.second];
		const BasisType* basisOld  = &(model_.basis());
		RealType         isign     = 1.0;
		PairType         oldParts  = model_.basis().parts();
		PairType         newParts  = oldParts;

		for (SizeType isite = 0; isite < sites.size(); ++isite) {
			SizeType site = sites[isite];
			if (orbs[isite] >= model_.orbitals(site))
				continue;

			const BasisType* basisNew = getNeededBasis(
			    newParts, oldParts, what[isite], spins[isite], orbs[isite]);

			if (!basisNew)
				return 0.0;

			VectorType modifVector(basisNew->size(), 0);
			accModifiedState_(modifVector,
			                  what[isite],
			                  *basisNew,
			                  tmpVector,
			                  *basisOld,
			                  site,
			                  spins[isite],
			                  orbs[isite],
			                  isign);

			tmpVector = modifVector;
			basisOld  = basisNew;
			oldParts  = newParts;
		}

		oldParts = model_.basis().parts();
		if (oldParts != newParts)
			return 0.0;

		checkBraOrKet("bra", braAndKet.first);
		const VectorType& braVector = vectors_[braAndKet.first];

		return braVector * tmpVector;
	}

	const BasisType* getNeededBasis(PairType&                  newParts,
	                                const PairType&            oldParts,
	                                const LabeledOperatorType& lOperator,
	                                SizeType                   spin,
	                                SizeType                   orb) const
	{
		if (!lOperator.needsNewBasis()) {
			newParts = oldParts;
			return &model_.basis();
		}

		if (!model_.hasNewParts(newParts, oldParts, lOperator, spin, orb))
			return 0;

		BasisType* basisNew = model_.createBasis(newParts.first, newParts.second);

		std::cerr << "basisNew.size=" << basisNew->size() << " ";
		std::cerr << "newparts.first=" << newParts.first << " ";
		std::cerr << "newparts.second=" << newParts.second << "\n";
		return basisNew;
	}

	void accModifiedState_(VectorType&                z,
	                       const LabeledOperatorType& lOperator,
	                       const BasisType&           newBasis,
	                       const VectorType&          srcVector,
	                       const BasisType&           srcBasis,
	                       SizeType                   site,
	                       SizeType                   spin,
	                       SizeType                   orb,
	                       ComplexOrRealType          factor) const
	{
		for (SizeType ispace = 0; ispace < srcBasis.size(); ispace++) {
			LanczosGlobals::WordType    ket1 = srcBasis(ispace, SPIN_UP);
			LanczosGlobals::WordType    ket2 = srcBasis(ispace, SPIN_DOWN);
			LanczosGlobals::PairIntType tempValue
			    = newBasis.getBraIndex(ket1, ket2, lOperator, site, spin, orb);
			int      temp  = tempValue.first;
			RealType value = tempValue.second;
			if (temp >= 0 && SizeType(temp) >= z.size()) {
				PsimagLite::String s = "old basis=" + ttos(srcBasis.size());
				s += " newbasis=" + ttos(newBasis.size());
				s += "\n";
				s += "operatorLabel= " + lOperator.toString()
				    + " spin=" + ttos(spin);
				s += " site=" + ttos(site);
				s += "ket1=" + ttos(ket1) + " and ket2=" + ttos(ket2);
				s += "\n";
				s += "accModifiedState_: z.size=" + ttos(z.size());
				s += " but temp=" + ttos(temp) + "\n";
				throw std::runtime_error(s.c_str());
			}
			if (temp < 0)
				continue;
			RealType mysign = (lOperator.isFermionic())
			    ? srcBasis.doSignGf(ket1, ket2, site, spin, orb)
			    : 1;
			if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SPLUS
			    || lOperator.id() == LabeledOperatorType::Label::OPERATOR_SMINUS)
				mysign *= srcBasis.doSignSpSm(ket1, ket2, site, spin, orb);

			z[temp] += factor * mysign * value * srcVector[ispace];
		}
	}

	void calcSpectral(ContFractionType&          cf,
	                  const bool                 isFermionic,
	                  const VectorType&          modifVector,
	                  const InternalProductType& matrix,
	                  SizeType                   type,
	                  SizeType,
	                  bool isDiagonal) const
	{
		ParametersForSolverType params(io_, "Spectral");

		LanczosSolverType lanczosSolver(matrix, params);

		PsimagLite::TridiagonalMatrix<RealType> ab;

		lanczosSolver.decomposition(modifVector, ab);
		typename VectorType::value_type weight = modifVector * modifVector;

		int      s  = (type & 1) ? -1 : 1;
		RealType s2 = (type > 1) ? -1 : 1;
		if (!isFermionic)
			s2 *= s;
		RealType diagonalFactor = (isDiagonal) ? 1 : 0.5;
		s2 *= diagonalFactor;

		assert(0 < energies_.size());
		const RealType gsEnergy = energies_[0];
		cf.set(ab, gsEnergy, PsimagLite::real(weight * s2), -s);
	}

private:

	void getModifiedState(VectorType&                          modifVector,
	                      const LabeledOperatorType&           lOperator,
	                      const VectorType&                    gsVector,
	                      const BasisType&                     basisNew,
	                      SizeType                             type,
	                      SizeType                             isite,
	                      SizeType                             jsite,
	                      SizeType                             spin,
	                      const std::pair<SizeType, SizeType>& orbs) const
	{
		modifVector.resize(basisNew.size());
		for (SizeType temp = 0; temp < modifVector.size(); temp++)
			modifVector[temp] = 0.0;

		accModifiedState_(modifVector,
		                  lOperator,
		                  basisNew,
		                  gsVector,
		                  model_.basis(),
		                  isite,
		                  spin,
		                  orbs.first,
		                  1.0);
		std::cerr << "isite=" << isite << " type=" << type;
		std::cerr << " modif=" << (modifVector * modifVector) << "\n";
		if (model_.name() == "Tj1Orb.h" && isite == jsite)
			return;

		RealType isign = (type > 1) ? -1.0 : 1.0;
		accModifiedState_(modifVector,
		                  lOperator,
		                  basisNew,
		                  gsVector,
		                  model_.basis(),
		                  jsite,
		                  spin,
		                  orbs.second,
		                  isign);
		std::cerr << "jsite=" << jsite << " type=" << type;
		std::cerr << " modif=" << (modifVector * modifVector) << "\n";
	}

	void accModifiedState(VectorType&                z,
	                      const LabeledOperatorType& lOperator,
	                      const BasisType&           newBasis,
	                      const VectorType&          gsVector,
	                      SizeType                   site,
	                      SizeType                   spin,
	                      SizeType                   orb,
	                      RealType                   isign) const
	{
		LabeledOperatorType opN(LabeledOperatorType::Label::OPERATOR_N);

		if (model_.name() == "Tj1Orb.h")
			accModifiedState_(z,
			                  lOperator,
			                  newBasis,
			                  gsVector,
			                  model_.basis(),
			                  site,
			                  spin,
			                  orb,
			                  isign);

		if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_N) {
			accModifiedState_(z,
			                  lOperator,
			                  newBasis,
			                  gsVector,
			                  model_.basis(),
			                  site,
			                  spin,
			                  orb,
			                  isign);
			return;
		} else if (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SZ) {
			accModifiedState_(z,
			                  opN,
			                  newBasis,
			                  gsVector,
			                  model_.basis(),
			                  site,
			                  SPIN_UP,
			                  orb,
			                  isign * 0.5);
			accModifiedState_(z,
			                  opN,
			                  newBasis,
			                  gsVector,
			                  model_.basis(),
			                  site,
			                  SPIN_DOWN,
			                  orb,
			                  -isign * 0.5);
			return;
		}

		accModifiedState_(
		    z, lOperator, newBasis, gsVector, model_.basis(), site, spin, orb, isign);
	}

	void computeAllStatesBelow(SizeType excited)
	{
		const SizeType excitedPlusOne = excited + 1;
		energies_.resize(excitedPlusOne);
		vectors_.resize(excitedPlusOne);

		SpecialSymmetryType     rs(model_.basis(), model_.geometry(), options_);
		InternalProductType     hamiltonian(model_, rs);
		ParametersForSolverType params(io_, "Lanczos");
		LanczosSolverType       lanczosSolver(hamiltonian, params);

		SizeType offset                 = model_.size();
		SizeType currentOffset          = 0;
		bool     firstNonZeroSectorSeen = false;

		for (SizeType i = 0; i < rs.sectors(); ++i) {
			hamiltonian.specialSymmetrySector(i);
			SizeType n = hamiltonian.rows();
			if (n == 0)
				continue;
			VectorType initial(n);
			PsimagLite::fillRandom(initial);
			VectorVectorType zs(excitedPlusOne, VectorType(n));
			VectorRealType   eigs(excitedPlusOne);

			try {
				lanczosSolver.computeAllStatesBelow(
				    eigs, zs, initial, excitedPlusOne);
			} catch (std::exception&) {

				std::cerr << "Engine: Lanczos Solver failed ";
				std::cerr << " trying exact diagonalization...\n";
				VectorRealType eigs2(n);
				MatrixType     fm;
				hamiltonian.fullDiag(eigs2, fm);
				for (SizeType k = 0; k < excitedPlusOne; ++k) {
					for (SizeType j = 0; j < n; ++j)
						zs[k][j] = fm(j, k);
					eigs[k] = eigs2[k];
				}
			}

			if (eigs[0] < energies_[0] || !firstNonZeroSectorSeen) {
				for (SizeType j = 0; j < excitedPlusOne; ++j) {
					vectors_[j]  = zs[j];
					energies_[j] = eigs[j];
				}

				offset                 = currentOffset;
				firstNonZeroSectorSeen = true;
			}

			currentOffset += zs[0].size();
		}

		rs.transform(vectors_, offset);

		printEnergiesAndNorms();
	}

	void checkBraOrKet(PsimagLite::String name, SizeType ind) const
	{
		if (ind < vectors_.size())
			return;

		err("Wrong " + name + " FATAL ERROR\n");
	}

	void printEnergiesAndNorms() const
	{
		const SizeType excited = energies_.size();
		assert(excited == vectors_.size());
		for (SizeType i = 0; i < excited; ++i) {
			const RealType val = PsimagLite::real(vectors_[i] * vectors_[i]);
			std::cout << "E[" << i << "]=" << energies_[i] << " norm=" << val << "\n";
		}
	}

	const ModelType&              model_;
	PsimagLite::ProgressIndicator progress_;
	InputType&                    io_;
	PsimagLite::String            options_;
	VectorRealType                energies_;
	VectorVectorType              vectors_;
}; // class ContinuedFraction
} // namespace Dmrg

#endif // ENGINE_H_
