#ifndef TERMFORTARGETINGEXPRESSION_H
#define TERMFORTARGETINGEXPRESSION_H

namespace Dmrg {

class TermForTargetingExpression {

	void finalize(SizeType termIndex)
	{
		TermType& term = terms_[termIndex];
		if (term.isFinalized()) return;

		SizeType n = vStr_.size();
		if (n == 0)
			err("AlgebraForTargetingExpression: Cannot finalize an empty object\n");

		PsimagLite::String ket;
		SizeType sitesEqualToCoo = 0;
		VectorSizeType discardedTerms;

		VectorStringType newVstr;
		for (SizeType i = 0; i < n; ++i) {
			PsimagLite::String tmp = vStr_[i];
			if (tmp[0] == '|') { // it's a vector
				if (ket != "")
					err("More than one ket found in " + toString() + "\n");
				ket = tmp;
				if (i + 1 != n)
					err("Vector is not at the end in " + toString() + "\n");

				newVstr.push_back(tmp);
				continue; // == break;
			}

			// it's a matrix or a scalar
			if (tmp[0] == '-' || tmp[0] == '+' || (tmp[0] >= 65 && tmp[0] <= 74)) {
				err("Scalars not supported yet\n");
			}

			// it's a matrix
			if (i != n - 2) {
				newVstr.push_back(tmp);
				continue; // apply in order only IMPORTANT
			              // the n - 1 is the ket
			}

			SiteSplitType siteSplit = OneOperatorSpecType::extractSiteIfAny(tmp);
			if (!siteSplit.hasSiteString)
				err("Each op must have a site\n");

			SizeType site = OneOperatorSpecType::strToNumberOfFail(siteSplit.siteString);
			if (site != coo) { // can only apply at center of orthogonality (coo)
				newVstr.push_back(tmp);
				continue;
			}

			if (sitesEqualToCoo >= 1)
				err("Sites cannot be repeated in term\n");
			++sitesEqualToCoo;
			discardedTerms.push_back(i);

			tmp = siteSplit.root;
			OneOperatorSpecType op = new OneOperatorSpecType(tmp);
			finalizeInternal(ket, op, site);
			delete op;
			op = 0;
		}

		// discard discared terms
		vStr_ = newVstr;
		finalizedInMult_ = true;

		//		const RealType oneReal = 1.0;
		//		if (factor_ != oneReal)
		//			fullVector_ = factor_*fullVector_;
		//		factor_ = 1.0;

		//		if (vwo) {
		//			*vwo = fullVector_;
		//			fullVector_.clear();
		//		}

		//		finalized_ = true;
	}

private:


	void oneOperator(PsimagLite::String ket, const OperatorType& op, SizeType site)
	{
		SizeType currentCoO = getCurrentCoO();

		if (site == currentCoO) {
			const VectorWithOffsetType& srcVwo = getCurrentVector(ket);
			applyInSitu(srcVwo, site, op);
			return;
		}

		err("oneOperator should never reach here\n");
		// Fetch ket at coo site --> into vec
		// Apply op to vec ---> vec2
		// Bring vec2 back to current coo
	}

	// returns A|src1>
	void applyInSitu(const VectorWithOffsetType& src1,
	                 SizeType site,
	                 const OperatorType& A)
	{
		err("applyInSitu unimplemented\n");

		typename PsimagLite::Vector<bool>::Type oddElectrons;
		aux_.model.findOddElectronsOfOneSite(oddElectrons,site);
		FermionSign fs(aux_.lrs.left(), oddElectrons);
		bool b1 = (site == 1 && aux_.direction == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON);
		SizeType n = aux_.model.superGeometry().numberOfSites();
		assert(n > 2);
		bool b2 = (site == n - 2 && aux_.direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM);
		BorderEnumType border = (b1 || b2) ? BorderEnumType::BORDER_YES
		                                   : BorderEnumType::BORDER_NO;
		aux_.aoe.applyOpLocal()(fullVector_, src1, A, fs, aux_.direction, border);
	}

	const VectorWithOffsetType& getCurrentVector(PsimagLite::String braOrKet) const
	{
		PsimagLite::GetBraOrKet getBraOrKet(braOrKet);

		if (getBraOrKet.isPvector()) {
			const SizeType pIndex = getBraOrKet.pIndex();
			if (pIndex >= aux_.pvectors.size())
				err("getVector: out of range for " + braOrKet + "\n");
			return aux_.pvectors[pIndex];
		}

		const SizeType sectorIndex = getBraOrKet.sectorIndex();
		return *(aux_.psi[sectorIndex][getBraOrKet.levelIndex()]);
	}

	SizeType getCurrentCoO() const
	{
		const LeftRightSuperType& lrs = aux_.lrs;
		const SizeType systemBlockSize = lrs.left().block().size();
		assert(systemBlockSize > 0);
		const int maxSystemSite = lrs.left().block()[systemBlockSize - 1];
		return (aux_.direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? maxSystemSite
		                                                                        : maxSystemSite + 1;

	}

	VectorStringType vStr_;
};
}
#endif // TERMFORTARGETINGEXPRESSION_H
