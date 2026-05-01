
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

#ifndef FEBASED_SC_H
#define FEBASED_SC_H

#include "Geometry/GeometryDca.h"
#include "LanczosModelBase.hpp"
#include "Parallelizer.h"
#include "ParametersModelFeAs.h"
#include "SparseRow.h"

namespace LanczosPlusPlus {

template <typename ComplexOrRealType, typename BasisType, typename InputType>
class FeBasedSc
    : public LanczosModelBase<ComplexOrRealType, typename BasisType::GeometryType, InputType> {

	typedef typename BasisType::GeometryType                             GeometryType;
	typedef FeBasedSc<ComplexOrRealType, BasisType, InputType>           ThisType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type           RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType>                        MatrixType;
	typedef LanczosModelBase<ComplexOrRealType, GeometryType, InputType> BaseType;
	typedef PsimagLite::GeometryDca<RealType, GeometryType>              GeometryDcaType;

	enum
	{
		SPIN_UP   = LanczosGlobals::SPIN_UP,
		SPIN_DOWN = LanczosGlobals::SPIN_DOWN
	};

	typedef ParametersModelFeAs<ComplexOrRealType> ParametersModelType;
	typedef typename BasisType::BaseType           BasisBaseType;
	typedef typename BasisType::WordType           WordType;
	typedef typename BaseType::SparseMatrixType    SparseMatrixType;
	typedef typename BaseType::VectorType          VectorType;

	class MatrixVectorHelper {

		typedef PsimagLite::Concurrency ConcurrencyType;

	public:

		MatrixVectorHelper(SizeType             nthreads,
		                   VectorType&          x,
		                   const VectorType&    y,
		                   const BasisBaseType& basis,
		                   const ThisType&      myself)
		    : nthreads_(nthreads)
		    , x_(x)
		    , y_(y)
		    , basis_(basis)
		    , myself_(myself)
		{ }

		void doTask(SizeType taskNumber, SizeType)
		{
			SizeType      nsite = myself_.geometry_.numberOfSites();
			SparseRowType sparseRow;

			WordType ket1 = basis_.operator()(taskNumber, SPIN_UP);
			WordType ket2 = basis_.operator()(taskNumber, SPIN_DOWN);

			x_[taskNumber] += myself_.findS(nsite, ket1, ket2, taskNumber, basis_)
			    * y_[taskNumber];

			for (SizeType i = 0; i < nsite; i++) {
				for (SizeType orb = 0; orb < myself_.mp_.orbitals; orb++) {
					myself_.setHoppingTerm(
					    sparseRow, ket1, ket2, i, orb, basis_);

					if (myself_.mp_.feAsMode
					    == ParametersModelType::IntEnum::INT_PAPER33) {
						myself_.setU2OffDiagonalTerm(
						    sparseRow, ket1, ket2, i, orb, basis_);

						for (SizeType orb2 = orb + 1;
						     orb2 < myself_.mp_.orbitals;
						     orb2++) {
							myself_.setU3Term(sparseRow,
							                  ket1,
							                  ket2,
							                  i,
							                  orb,
							                  orb2,
							                  basis_);
						}

						myself_.setJTermOffDiagonal(
						    sparseRow, ket1, ket2, i, orb, basis_);
					} else if (myself_.mp_.feAsMode
					               == ParametersModelType::IntEnum::INT_V
					           || myself_.mp_.feAsMode
					               == ParametersModelType::IntEnum::INT_CODE2) {
						myself_.setOffDiagonalDecay(
						    sparseRow, ket1, ket2, i, orb, basis_);
					} else if (myself_.mp_.feAsMode
					           == ParametersModelType::IntEnum::INT_IMPURITY) {
						myself_.setOffDiagonalJimpurity(
						    sparseRow, ket1, ket2, i, orb, basis_);
					} else if (myself_.mp_.feAsMode
					           == ParametersModelType::IntEnum::INT_KSPACE) {
						myself_.setOffDiagonalKspace(
						    sparseRow, ket1, ket2, i, orb, basis_);
					}
				}
			}

			x_[taskNumber] += sparseRow.finalize(y_);
		}

		SizeType tasks() const { return x_.size(); }

	private:

		SizeType             nthreads_;
		VectorType&          x_;
		const VectorType&    y_;
		const BasisBaseType& basis_;
		const ThisType&      myself_;
	}; // class MatrixVectorHelper

public:

	typedef PsimagLite::SparseRow<SparseMatrixType> SparseRowType;

	enum TermEnum
	{
		HOPPING,
		J_PM,
		J_ZZ
	};

	FeBasedSc(SizeType nup, SizeType ndown, InputType& io, const GeometryType& geometry)
	    : mp_(io)
	    , geometry_(geometry)
	    , basis_(geometry, nup, ndown, mp_.orbitals)
	    , geometryDca_(geometry, mp_.orbitals)
	{
		checkSpinOrbit();
	}

	~FeBasedSc() { BaseType::deleteGarbage(garbage_); }

	SizeType size() const { return basis_.size(); }

	SizeType orbitals(SizeType) const { return mp_.orbitals; }

	void setupHamiltonian(SparseMatrixType& matrix) const { setupHamiltonian(matrix, basis_); }

	bool hasNewParts(std::pair<SizeType, SizeType>&       newParts,
	                 const std::pair<SizeType, SizeType>& oldParts,
	                 const LabeledOperator&               lOperator,
	                 SizeType                             spin,
	                 SizeType                             orb) const
	{
		return basis_.hasNewParts(newParts, oldParts, lOperator, spin, orb);
	}

	const GeometryType& geometry() const { return geometry_; }

	void setupHamiltonian(SparseMatrixType& matrix, const BasisBaseType& basis) const
	{
		// Calculate diagonal elements AND count non-zero matrix elements
		SizeType                                    hilbert = basis.size();
		typename PsimagLite::Vector<RealType>::Type diag(hilbert);
		calcDiagonalElements(diag, basis);

		SizeType nsite = geometry_.numberOfSites();

		// Setup CRS matrix
		matrix.resize(hilbert, hilbert);

		// Calculate off-diagonal elements AND store matrix
		SizeType nCounter = 0;
		for (SizeType ispace = 0; ispace < hilbert; ispace++) {
			SparseRowType sparseRow;
			matrix.setRow(ispace, nCounter);
			WordType ket1 = basis(ispace, SPIN_UP);
			WordType ket2 = basis(ispace, SPIN_DOWN);
			// Save diagonal
			sparseRow.add(ispace, diag[ispace]);
			for (SizeType i = 0; i < nsite; i++) {
				for (SizeType orb = 0; orb < mp_.orbitals; orb++) {
					setHoppingTerm(sparseRow, ket1, ket2, i, orb, basis);

					if (mp_.feAsMode
					    == ParametersModelType::IntEnum::INT_PAPER33) {
						setU2OffDiagonalTerm(
						    sparseRow, ket1, ket2, i, orb, basis);
						for (SizeType orb2 = 0; orb2 < mp_.orbitals;
						     orb2++) {
							if (orb == orb2)
								continue;

							setU3Term(sparseRow,
							          ket1,
							          ket2,
							          i,
							          orb,
							          orb2,
							          basis);
						}

						setJTermOffDiagonal(
						    sparseRow, ket1, ket2, i, orb, basis);

						setSpinOrbitOffDiagonal(
						    sparseRow, ket1, ket2, i, orb, basis);

					} else if (mp_.feAsMode
					               == ParametersModelType::IntEnum::INT_V
					           || mp_.feAsMode
					               == ParametersModelType::IntEnum::INT_CODE2) {
						setOffDiagonalDecay(
						    sparseRow, ket1, ket2, i, orb, basis);
					} else if (mp_.feAsMode
					           == ParametersModelType::IntEnum::INT_IMPURITY) {
						setOffDiagonalJimpurity(
						    sparseRow, ket1, ket2, i, orb, basis);
					} else if (mp_.feAsMode
					           == ParametersModelType::IntEnum::INT_KSPACE) {
						setOffDiagonalKspace(
						    sparseRow, ket1, ket2, i, orb, basis);
					}
				}
			}

			nCounter += sparseRow.finalize(matrix);
		}

		matrix.setRow(hilbert, nCounter);
	}

	void matrixVectorProduct(VectorType& x, const VectorType& y) const
	{
		matrixVectorProduct(x, y, basis_);
	}

	void
	matrixVectorProduct(VectorType& x, const VectorType& y, const BasisBaseType& basis) const
	{
		// Calculate off-diagonal elements AND store matrix
		typedef MatrixVectorHelper                   HelperType;
		typedef PsimagLite::Parallelizer<HelperType> ParallelizerType;
		ParallelizerType threadObject(PsimagLite::Concurrency::codeSectionParams);
		HelperType       helper(
                    PsimagLite::Concurrency::codeSectionParams.npthreads, x, y, basis, *this);

		std::cout << "Using " << threadObject.name();
		std::cout << " with " << PsimagLite::Concurrency::codeSectionParams.npthreads
		          << " threads.\n";
		threadObject.loopCreate(helper);
	}

	const BasisType& basis() const { return basis_; }

	PsimagLite::String name() const { return __FILE__; }

	BasisType* createBasis(SizeType nup, SizeType ndown) const
	{
		BasisType* ptr = new BasisType(geometry_, nup, ndown, mp_.orbitals);
		garbage_.push_back(ptr);
		return ptr;
	}

	void print(std::ostream& os) const { os << mp_; }

private:

	void setOffDiagonalDecay(SparseRowType&  sparseRow,
	                         const WordType& ket1,
	                         const WordType& ket2,
	                         SizeType        i,
	                         SizeType,
	                         const BasisBaseType& basis) const
	{
		for (SizeType spin1 = 0; spin1 < 2; ++spin1) {
			for (SizeType spin2 = 0; spin2 < 2; ++spin2) {

				if (spin1 != spin2)
					continue;

				if (!basis.isThereAnElectronAt(ket1, ket2, i, spin1, 1))
					continue;
				if (!basis.isThereAnElectronAt(ket1, ket2, i, spin2, 1))
					continue;
				WordType mask = BasisType::bitmask(i * mp_.orbitals + 1);
				WordType bra1 = (spin1 == SPIN_UP) ? (ket1 ^ mask) : ket1;
				WordType bra2 = (spin1 == SPIN_UP) ? ket2 : (ket2 ^ mask);

				if (!basis.isThereAnElectronAt(bra1, bra2, i, spin2, 1))
					continue;
				WordType bra3 = (spin2 == SPIN_UP) ? (bra1 ^ mask) : bra1;
				WordType bra4 = (spin2 == SPIN_UP) ? bra2 : (bra2 ^ mask);

				if (basis.isThereAnElectronAt(bra3, bra4, i, spin2, 0))
					continue;
				mask          = BasisType::bitmask(i * mp_.orbitals + 0);
				WordType bra5 = (spin2 == SPIN_UP) ? (bra3 ^ mask) : bra3;
				WordType bra6 = (spin2 == SPIN_UP) ? bra4 : (bra4 ^ mask);

				if (basis.isThereAnElectronAt(bra5, bra6, i, spin1, 2))
					continue;
				mask          = BasisType::bitmask(i * mp_.orbitals + 2);
				WordType bra7 = (spin1 == SPIN_UP) ? (bra5 ^ mask) : bra5;
				WordType bra8 = (spin1 == SPIN_UP) ? bra6 : (bra6 ^ mask);

				SizeType temp = basis.perfectIndex(bra7, bra8);
				sparseRow.add(temp, mp_.coulombV);
			}
		}
	}

	RealType findSdecay(SizeType,
	                    WordType             ket1,
	                    WordType             ket2,
	                    SizeType             i,
	                    SizeType             orb,
	                    const BasisBaseType& basis) const
	{
		RealType s = mp_.hubbardU[orb + orb * mp_.orbitals]
		    * basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb)
		    * basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb);

		for (SizeType orb2 = orb + 1; orb2 < mp_.orbitals; orb2++) {
			// Hubbard term U1
			s += mp_.hubbardU[orb + orb2 * mp_.orbitals]
			    * nix(ket1, ket2, i, orb, basis) * nix(ket1, ket2, i, orb2, basis);
		}

		return s;
	}

	ComplexOrRealType hoppings(SizeType i, SizeType orb1, SizeType j, SizeType orb2) const
	{
		return -geometry_(i, orb1, j, orb2, TermEnum::HOPPING);
	}

	void setHoppingTerm(SparseRowType&       sparseRow,
	                    const WordType&      ket1,
	                    const WordType&      ket2,
	                    SizeType             i,
	                    SizeType             orb,
	                    const BasisBaseType& basis) const
	{
		SizeType ii  = i * mp_.orbitals + orb;
		WordType s1i = (ket1 & BasisType::bitmask(ii));
		if (s1i > 0)
			s1i = 1;
		WordType s2i = (ket2 & BasisType::bitmask(ii));
		if (s2i > 0)
			s2i = 1;

		SizeType nsite = geometry_.numberOfSites();

		// Hopping term
		for (SizeType j = 0; j < nsite; j++) {
			if (j < i)
				continue;
			for (SizeType orb2 = 0; orb2 < mp_.orbitals; orb2++) {
				SizeType          jj = j * mp_.orbitals + orb2;
				ComplexOrRealType h  = hoppings(i, orb, j, orb2);
				if (PsimagLite::real(h) == 0 && PsimagLite::imag(h) == 0)
					continue;
				WordType s1j = (ket1 & BasisType::bitmask(jj));
				if (s1j > 0)
					s1j = 1;
				WordType s2j = (ket2 & BasisType::bitmask(jj));
				if (s2j > 0)
					s2j = 1;

				if (s1i + s1j == 1) {
					if (s1i == 0)
						h = PsimagLite::conj(h);
					WordType bra1 = ket1
					    ^ (BasisType::bitmask(ii) | BasisType::bitmask(jj));
					SizeType temp = basis.perfectIndex(bra1, ket2);
					RealType extraSign
					    = (s1i == 1) ? LanczosGlobals::FERMION_SIGN : 1;
					RealType tmp2
					    = basis_.doSign(ket1, ket2, i, orb, j, orb2, SPIN_UP);
					ComplexOrRealType cTemp = h * extraSign * tmp2;
					sparseRow.add(temp, cTemp);
				}

				if (s2i + s2j == 1) {
					if (s2i == 0)
						h = PsimagLite::conj(h);
					WordType bra2 = ket2
					    ^ (BasisType::bitmask(ii) | BasisType::bitmask(jj));
					SizeType temp = basis.perfectIndex(ket1, bra2);
					RealType extraSign
					    = (s2i == 1) ? LanczosGlobals::FERMION_SIGN : 1;
					RealType tmp2
					    = basis_.doSign(ket1, ket2, i, orb, j, orb2, SPIN_DOWN);
					ComplexOrRealType cTemp = h * extraSign * tmp2;
					sparseRow.add(temp, cTemp);
				}
			}
		}
	}

	void setU2OffDiagonalTerm(SparseRowType&       sparseRow,
	                          const WordType&      ket1,
	                          const WordType&      ket2,
	                          SizeType             i,
	                          SizeType             orb1,
	                          const BasisBaseType& basis) const
	{
		RealType val = mp_.hubbardU[2] * 0.5;
		for (SizeType orb2 = 0; orb2 < mp_.orbitals; orb2++) {
			if (orb1 == orb2)
				continue;
			RealType sign = jTermSign(ket1, ket2, i, orb1, i, orb2, basis);
			setSplusSminus(sparseRow, ket1, ket2, i, orb1, i, orb2, val * sign, basis);
		}
	}

	// N.B.: orb1!=orb2 here
	void setSplusSminus(SparseRowType&       sparseRow,
	                    const WordType&      ket1,
	                    const WordType&      ket2,
	                    SizeType             i,
	                    SizeType             orb1,
	                    SizeType             j,
	                    SizeType             orb2,
	                    ComplexOrRealType    value,
	                    const BasisBaseType& basis) const
	{
		if (splusSminusNonZero(ket1, ket2, i, orb1, j, orb2, basis) == 0)
			return;

		SizeType ii = i * mp_.orbitals + orb1;
		SizeType jj = j * mp_.orbitals + orb2;
		assert(ii != jj);
		WordType bra1 = ket1 ^ (BasisType::bitmask(ii) | BasisType::bitmask(jj));
		WordType bra2 = ket2 ^ (BasisType::bitmask(ii) | BasisType::bitmask(jj));
		SizeType temp = basis.perfectIndex(bra1, bra2);
		sparseRow.add(temp, value);
	}

	// N.B.: orb1!=orb2 here
	void setU3Term(SparseRowType&       sparseRow,
	               const WordType&      ket1,
	               const WordType&      ket2,
	               SizeType             i,
	               SizeType             orb1,
	               SizeType             orb2,
	               const BasisBaseType& basis) const
	{
		assert(orb1 != orb2);
		if (u3TermNonZero(ket1, ket2, i, orb1, orb2, basis) == 0)
			return;

		SizeType ii   = i * mp_.orbitals + orb1;
		SizeType jj   = i * mp_.orbitals + orb2;
		WordType bra1 = ket1 ^ (BasisType::bitmask(ii) | BasisType::bitmask(jj));
		WordType bra2 = ket2 ^ (BasisType::bitmask(ii) | BasisType::bitmask(jj));
		SizeType temp = basis.perfectIndex(bra1, bra2);
		RealType sign = jTermSign(ket1, ket2, i, orb1, i, orb2, basis);
		sparseRow.add(temp, LanczosGlobals::FERMION_SIGN * mp_.hubbardU[3] * sign);
	}

	void setSpinOrbitOffDiagonal(SparseRowType&       sparseRow,
	                             const WordType&      ketUp,
	                             const WordType&      ketDown,
	                             SizeType             i,
	                             SizeType             orb1,
	                             const BasisBaseType& basis) const
	{
		const RealType zero = 0;
		if (mp_.spinOrbit.n_row() != 4)
			return;
		SizeType orbitals = mp_.orbitals;
		SizeType i1       = i * mp_.orbitals + orb1;
		for (SizeType orb2 = 0; orb2 < orbitals; orb2++) {
			SizeType i2 = i * mp_.orbitals + orb2;
			for (SizeType spin1 = 0; spin1 < 2; spin1++) {

				WordType ket1 = (spin1) ? ketDown : ketUp;

				for (SizeType spin2 = 0; spin2 < 2; spin2++) {

					WordType          ket2  = (spin2) ? ketDown : ketUp;
					WordType          bra1  = ket1;
					WordType          bra2  = ket2;
					ComplexOrRealType value = mp_.spinOrbit(
					    spin1 + spin2 * 2, orb1 + orb2 * orbitals);

					WordType s12 = (ket1 & BasisType::bitmask(i1));
					if (s12 > 0)
						s12 = 1;
					WordType s21 = (ket2 & BasisType::bitmask(i2));
					if (s21 > 0)
						s21 = 1;
					if (s12 == 0 || s21 == 1)
						continue;
					if (spin1 == spin2) {
						bra1 = ket1
						    ^ (BasisType::bitmask(i2)
						       | BasisType::bitmask(i1));
						bra2 = (spin1) ? ketUp : ketDown;
					} else {
						bra1 = ket1 ^ (BasisType::bitmask(i1));
						bra2 = ket2 ^ (BasisType::bitmask(i2));
					}

					SizeType temp = (spin1) ? basis.perfectIndex(bra2, bra1)
					                        : basis.perfectIndex(bra1, bra2);

					RealType s = basis_.doSignSpinOrbit(
					    ket1, ket2, i, spin1, orb1, spin2, orb2);
					ComplexOrRealType cTemp = value * s;
					if (cTemp == zero)
						continue;
					sparseRow.add(temp, cTemp);
				}
			}
		}
	}

	void setJTermOffDiagonal(SparseRowType&       sparseRow,
	                         const WordType&      ket1,
	                         const WordType&      ket2,
	                         SizeType             i,
	                         SizeType             orb,
	                         const BasisBaseType& basis) const
	{
		const RealType zeroPointFive = 0.5;

		for (SizeType j = 0; j < geometry_.numberOfSites(); j++) {
			ComplexOrRealType value = jCoupling(i, j, TermEnum::J_PM) * zeroPointFive;
			if (PsimagLite::real(value) == 0 && PsimagLite::imag(value) == 0)
				continue;
			value *= 0.5; // double counting i,j
			assert(i != j);
			for (SizeType orb2 = 0; orb2 < mp_.orbitals; orb2++) {
				RealType sign = jTermSign(ket1, ket2, i, orb, j, orb2, basis);
				setSplusSminus(
				    sparseRow, ket1, ket2, i, orb, j, orb2, value * sign, basis);
				setSplusSminus(
				    sparseRow, ket1, ket2, j, orb2, i, orb, value * sign, basis);
			}
		}
	}

	int jTermSign(const WordType&      ket1,
	              const WordType&      ket2,
	              SizeType             i,
	              SizeType             orb1,
	              SizeType             j,
	              SizeType             orb2,
	              const BasisBaseType& basis) const
	{
		if (i > j)
			return jTermSign(ket1, ket2, j, orb2, i, orb1, basis);
		int x = basis.doSign(ket1, ket2, i, orb1, j, orb2, SPIN_UP);
		x *= basis.doSign(ket1, ket2, i, orb1, j, orb2, SPIN_DOWN);
		return x;
	}

	void calcDiagonalElements(typename PsimagLite::Vector<RealType>::Type& diag,
	                          const BasisBaseType&                         basis) const
	{
		SizeType hilbert = basis.size();
		SizeType nsite   = geometry_.numberOfSites();

		// Calculate diagonal elements
		for (SizeType ispace = 0; ispace < hilbert; ispace++) {
			WordType ket1 = basis(ispace, SPIN_UP);
			WordType ket2 = basis(ispace, SPIN_DOWN);
			diag[ispace]  = findS(nsite, ket1, ket2, ispace, basis);
		}
	}

	RealType findS(SizeType nsite,
	               WordType ket1,
	               WordType ket2,
	               SizeType,
	               const BasisBaseType& basis) const
	{
		RealType s = 0;
		for (SizeType i = 0; i < nsite; i++) {
			RealType szOrb = 0;

			for (SizeType orb = 0; orb < mp_.orbitals; orb++) {

				if (mp_.feAsMode == ParametersModelType::IntEnum::INT_PAPER33) {
					s += findSnoDecay(nsite, ket1, ket2, i, orb, basis);
				} else if (mp_.feAsMode == ParametersModelType::IntEnum::INT_V
				           || mp_.feAsMode
				               == ParametersModelType::IntEnum::INT_CODE2) {
					s += findSdecay(nsite, ket1, ket2, i, orb, basis);
				} else if (mp_.feAsMode
				           == ParametersModelType::IntEnum::INT_IMPURITY) {
					s += findSImpurity(nsite, ket1, ket2, i, orb, basis);
				} else if (mp_.feAsMode
				           == ParametersModelType::IntEnum::INT_KSPACE) {
					s += findSkspace(nsite, ket1, ket2, i, orb, basis);
				}

				// Potential term
				s += mp_.potentialV[i + (orb + mp_.orbitals * 0) * nsite]
				        * basis.getN(ket1, ket1, i, SPIN_UP, orb)
				    + mp_.potentialV[i + (orb + mp_.orbitals * 1) * nsite]
				        * basis.getN(ket2, ket2, i, SPIN_DOWN, orb);

				szOrb += szTerm(ket1, ket2, i, orb, basis);
			}

			s += mp_.anisotropyD * szOrb * szOrb;
		}

		return s;
	}

	RealType findSnoDecay(SizeType             nsite,
	                      WordType             ket1,
	                      WordType             ket2,
	                      SizeType             i,
	                      SizeType             orb,
	                      const BasisBaseType& basis) const
	{
		const RealType zeroPointFive = 0.5;
		// Hubbard term U0
		ComplexOrRealType s = mp_.hubbardU[0]
		    * basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb)
		    * basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb);

		for (SizeType orb2 = orb + 1; orb2 < mp_.orbitals; orb2++) {
			// Hubbard term U1
			s += mp_.hubbardU[1] * nix(ket1, ket2, i, orb, basis)
			    * nix(ket1, ket2, i, orb2, basis);

			// Diagonal U2 term
			s += mp_.hubbardU[4] * szTerm(ket1, ket2, i, orb, basis)
			    * szTerm(ket1, ket2, i, orb2, basis);

			s += mp_.hubbardU[5]
			    * basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb)
			    * basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb2);

			s += mp_.hubbardU[5]
			    * basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb)
			    * basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb2);
		}

		// JNN and JNNN diagonal part
		for (SizeType j = 0; j < nsite; j++) {
			for (SizeType orb2 = 0; orb2 < mp_.orbitals; orb2++) {
				ComplexOrRealType value = jCoupling(i, j, TermEnum::J_ZZ);
				if (PsimagLite::real(value) == 0 && PsimagLite::imag(value) == 0)
					continue;
				s += value * zeroPointFive * // RealType counting i,j
				    szTerm(ket1, ket2, i, orb, basis)
				    * szTerm(ket1, ket2, j, orb2, basis);
			}
		}

		// spin orbit diagonal in the matrix part
		if (mp_.spinOrbit.n_row() == 4) {
			for (SizeType spin = 0; spin < 2; ++spin)
				s += PsimagLite::real(
				         mp_.spinOrbit(spin + spin * 2, orb + orb * mp_.orbitals))
				    * basis.getN(ket1, ket2, i, spin, orb);
		}

		assert(fabs(PsimagLite::imag(s)) < 1e-12);
		return PsimagLite::real(s);
	}

	RealType findSImpurity(SizeType,
	                       WordType             ket1,
	                       WordType             ket2,
	                       SizeType             i,
	                       SizeType             orb,
	                       const BasisBaseType& basis) const
	{
		if (i > 0)
			return 0.0;

		// Hubbard term U0
		RealType s = mp_.hubbardU[0]
		    * basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb)
		    * basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb);

		for (SizeType orb2 = 0; orb2 < mp_.orbitals; orb2++) {
			// Hubbard term U1
			if (orb == orb2)
				continue;
			for (SizeType spin = 0; spin < 2; ++spin)
				s += 0.5 * mp_.hubbardU[1]
				    * basis.isThereAnElectronAt(ket1, ket2, i, spin, orb)
				    * basis.isThereAnElectronAt(ket1, ket2, i, spin, orb2);
		}

		for (SizeType orb2 = 0; orb2 < mp_.orbitals; orb2++) {
			if (orb == orb2)
				continue;
			// Diagonal U2 term
			s += mp_.hubbardU[4]
			    * basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb)
			    * basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb2);
		}

		return s;
	}

	RealType findSkspace(SizeType,
	                     WordType             ket1,
	                     WordType             ket2,
	                     SizeType             i,
	                     SizeType             orb,
	                     const BasisBaseType& basis) const
	{
		if (i > 0)
			return 0.0;
		RealType s   = 0;
		SizeType ck1 = basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb);
		if (ck1 == 0)
			return 0.0;

		for (SizeType orb2 = 0; orb2 < mp_.orbitals; orb2++) {
			SizeType ck2 = basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb2);
			if (ck2 == 0)
				continue;
			s++;
		}

		return s * mp_.hubbardU[0];
	}

	SizeType splusSminusNonZero(const WordType&      ket1,
	                            const WordType&      ket2,
	                            SizeType             i,
	                            SizeType             orb1,
	                            SizeType             j,
	                            SizeType             orb2,
	                            const BasisBaseType& basis) const
	{
		if (basis.isThereAnElectronAt(ket1, ket2, j, SPIN_UP, orb2) == 0)
			return 0;
		if (basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb1) == 1)
			return 0;
		if (basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb1) == 0)
			return 0;
		if (basis.isThereAnElectronAt(ket1, ket2, j, SPIN_DOWN, orb2) == 1)
			return 0;
		return 1;
	}

	SizeType u3TermNonZero(const WordType&      ket1,
	                       const WordType&      ket2,
	                       SizeType             i,
	                       SizeType             orb1,
	                       SizeType             orb2,
	                       const BasisBaseType& basis) const
	{
		if (basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb2) == 0)
			return 0;
		if (basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb1) == 1)
			return 0;
		if (basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb1) == 1)
			return 0;
		if (basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb2) == 0)
			return 0;
		return 1;
	}

	SizeType nix(const WordType&      ket1,
	             const WordType&      ket2,
	             SizeType             i,
	             SizeType             orb,
	             const BasisBaseType& basis) const
	{
		SizeType sum = 0;
		for (SizeType spin = 0; spin < 2; spin++)
			sum += basis.isThereAnElectronAt(ket1, ket2, i, spin, orb);
		return sum;
	}

	RealType szTerm(const WordType&      ket1,
	                const WordType&      ket2,
	                SizeType             i,
	                SizeType             orb,
	                const BasisBaseType& basis) const
	{
		RealType sz = basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb);
		sz -= basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb);
		return 0.5 * sz;
	}

	ComplexOrRealType jCoupling(SizeType i, SizeType j, SizeType term) const
	{
		if (geometry_.terms() == 1)
			return 0.0;
		return geometry_(i, 0, j, 0, term);
	}

	void setOffDiagonalJimpurity(SparseRowType&       sparseRow,
	                             const WordType&      ket1,
	                             const WordType&      ket2,
	                             SizeType             i,
	                             SizeType             orb1,
	                             const BasisBaseType& basis) const
	{
		if (i > 0)
			return;

		SizeType orbitals = mp_.orbitals;

		for (SizeType type = 0; type < 2; ++type) {
			for (SizeType orb2 = 0; orb2 < orbitals; ++orb2) {

				if (orb1 == orb2)
					continue;

				SizeType orb3 = (type == 0) ? orb2 : orb1;
				SizeType orb4 = (type == 0) ? orb1 : orb2;

				if (!basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb4))
					continue;
				if (basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb3))
					continue;
				if (!basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb2))
					continue;
				if (basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb1))
					continue;

				WordType mask4 = BasisType::bitmask(i * mp_.orbitals + orb4);
				WordType mask3 = BasisType::bitmask(i * mp_.orbitals + orb3);
				WordType bra2  = (ket2 ^ mask4) ^ mask3;

				WordType mask2 = BasisType::bitmask(i * mp_.orbitals + orb2);
				WordType mask1 = BasisType::bitmask(i * mp_.orbitals + orb1);
				WordType bra1  = (ket1 ^ mask2) ^ mask1;

				RealType x = basis.doSign(ket1, ket2, i, orb1, i, orb2, SPIN_UP);
				x *= basis.doSign(ket1, ket2, i, orb3, i, orb4, SPIN_DOWN);

				SizeType temp = basis.perfectIndex(bra1, bra2);
				sparseRow.add(temp, x * mp_.hubbardU[3]);
			}
		}
	}

	void setOffDiagonalKspace(SparseRowType&       sparseRow,
	                          const WordType&      ket1,
	                          const WordType&      ket2,
	                          SizeType             i,
	                          SizeType             orb1,
	                          const BasisBaseType& basis) const
	{
		if (i > 0)
			return;

		SizeType orbitals = mp_.orbitals;

		if (basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb1))
			return;

		for (SizeType orb2 = 0; orb2 < orbitals; ++orb2) {
			if (orb1 == orb2)
				continue;

			if (!basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb2))
				continue;

			for (SizeType orb3 = 0; orb3 < orbitals; ++orb3) {

				if (basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb3))
					continue;

				SizeType orb4 = getMomentum(orb1, orb2, orb3);
				assert(orb3 != orb4);

				if (!basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb4))
					continue;

				WordType mask4 = BasisType::bitmask(i * mp_.orbitals + orb4);
				WordType mask3 = BasisType::bitmask(i * mp_.orbitals + orb3);
				WordType bra2  = (ket2 ^ mask4) ^ mask3;

				WordType mask2 = BasisType::bitmask(i * mp_.orbitals + orb2);
				WordType mask1 = BasisType::bitmask(i * mp_.orbitals + orb1);
				WordType bra1  = (ket1 ^ mask2) ^ mask1;

				RealType x = basis.doSign(ket1, ket2, i, orb1, i, orb2, SPIN_UP);
				x *= basis.doSign(ket1, ket2, i, orb3, i, orb4, SPIN_DOWN);

				SizeType temp = basis.perfectIndex(bra1, bra2);
				sparseRow.add(temp, x * mp_.hubbardU[0]);
			}
		}
	}

	SizeType getMomentum(SizeType orb1, SizeType orb2, SizeType orb3) const
	{
		assert(orb1 < mp_.orbitals);
		assert(orb2 < mp_.orbitals);
		assert(orb3 < mp_.orbitals);

		SizeType tmp  = geometryDca_.kSum(orb3, orb1);
		SizeType orb4 = geometryDca_.kSustract(tmp, orb2);

		assert(orb4 < mp_.orbitals);
		return orb4;
	}

	void checkSpinOrbit() const
	{
		if (mp_.spinOrbit.n_row() == 0)
			return;
		if (mp_.spinOrbit.n_row() != 4)
			throw PsimagLite::RuntimeError("SpinOrbit must have 4 rows\n");

		for (SizeType spin = 0; spin < 2; ++spin) {
			for (SizeType orb = 0; orb < mp_.orbitals; ++orb) {
				RealType s = PsimagLite::imag(
				    mp_.spinOrbit(spin + spin * 2, orb + orb * mp_.orbitals));
				if (s == 0.0)
					continue;
				throw PsimagLite::RuntimeError("SpinOrbit term not Hermitian\n");
			}
		}
	}

	const ParametersModelType                             mp_;
	const GeometryType&                                   geometry_;
	BasisType                                             basis_;
	GeometryDcaType                                       geometryDca_;
	mutable typename PsimagLite::Vector<BasisType*>::Type garbage_;
}; // class FeBasedSc

} // namespace LanczosPlusPlus
#endif
