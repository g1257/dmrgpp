/*
Copyright (c) 2009-2014, UT-Battelle, LLC
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

#ifndef REFLECTION_SYMM_H
#define REFLECTION_SYMM_H
#include "CrsMatrix.h"
#include "LanczosGlobals.h"
#include "ProgressIndicator.h"
#include "Vector.h"
#include <iostream>

namespace LanczosPlusPlus {

class ReflectionItem {

public:

	enum
	{
		DIAGONAL,
		PLUS,
		MINUS
	};

	ReflectionItem(SizeType ii)
	    : i(ii)
	    , j(ii)
	    , type(DIAGONAL)
	{ }

	ReflectionItem(SizeType ii, SizeType jj, SizeType type1)
	    : i(ii)
	    , j(jj)
	    , type(type1)
	{ }

	SizeType i, j, type;

}; // class ReflectionItem

bool operator==(const ReflectionItem& item1, const ReflectionItem& item2);

template <typename GeometryType_, typename BasisType> class ReflectionSymmetry {

	typedef typename GeometryType_::ComplexOrRealType          ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType>              MatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef LanczosGlobals::WordType                           WordType;
	typedef ReflectionItem                                     ItemType;

public:

	typedef GeometryType_                                        GeometryType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType>             SparseMatrixType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename PsimagLite::Vector<RealType>::Type          VectorRealType;
	typedef typename PsimagLite::Vector<VectorType>::Type        VectorVectorType;

	ReflectionSymmetry(const BasisType&    basis,
	                   const GeometryType& geometry,
	                   PsimagLite::String  options)
	    : progress_("ReflectionSymmetry")
	    , transform_(basis.size(), basis.size())
	    , plusSector_(0)
	    , matrixStored_(2)
	    , pointer_(0)
	    , printMatrix_(options.find("printmatrix") != PsimagLite::String::npos)
	{
		SizeType hilbert       = basis.size();
		SizeType numberOfDofs  = basis.dofs();
		SizeType numberOfSites = geometry.numberOfSites();
		SizeType termId        = 0;
		//			SizeType counter=0;
		PsimagLite::Vector<ItemType>::Type buffer;
		for (SizeType ispace = 0; ispace < hilbert; ispace++) {
			typename PsimagLite::Vector<WordType>::Type y(numberOfDofs, 0);
			for (SizeType dof = 0; dof < numberOfDofs; dof++) {
				WordType x = basis(ispace, dof);
				for (SizeType site = 0; site < numberOfSites; site++) {
					SizeType reflectedSite
					    = geometry.term(termId).findReflection(site);
					SizeType thisSiteContent = x & 1;
					x >>= 1; // go to next site
					addTo(y[dof], thisSiteContent, reflectedSite);
					if (!x)
						break;
				}
			}

			SizeType yIndex = basis.perfectIndex(y);
			//				s_.setRow(ispace,counter);
			//				s_.pushCol(yIndex);
			//				s_.pushValue(1.0);
			//				counter++;
			if (yIndex == ispace) { // then S|psi> = |psi>
				ItemType item1(ispace);
				buffer.push_back(item1);
				continue;
			}
			// S|psi> != |psi>
			// Add normalized +
			ItemType item2(ispace, yIndex, ItemType::PLUS);
			buffer.push_back(item2);

			// Add normalized -
			ItemType item3(ispace, yIndex, ItemType::MINUS);
			buffer.push_back(item3);
		}
		//			s_.setRow(s_.rank(),counter);
		setTransform(buffer);
		//			checkTransform();
	}

	template <typename SomeModelType>
	void init(const SomeModelType& model, const BasisType& basis)
	{
		SparseMatrixType matrix2;
		model.setupHamiltonian(matrix2, basis);
		transformMatrix(matrixStored_, matrix2);

		if (matrixStored_.size() == 0)
			return;

		int nrows = matrixStored_[0].rows();
		if (printMatrix_) {
			if (nrows > 40)
				throw PsimagLite::RuntimeError("printMatrix too big\n");
			std::cout << matrixStored_[0].toDense();
		}
	}

	SizeType rows() const { return matrixStored_[pointer_].rows(); }

	void transformMatrix(typename PsimagLite::Vector<SparseMatrixType>::Type& matrix1,
	                     const SparseMatrixType&                              matrix) const
	{
		SparseMatrixType rT;
		transposeConjugate(rT, transform_);

		if (matrix.rows() < 40)
			printFullMatrix(matrix, "originalHam");
		SparseMatrixType tmp;
		multiply(tmp, matrix, rT);

		SparseMatrixType matrix2;
		multiply(matrix2, transform_, tmp);

		assert(matrix1.size() == 2);
		split(matrix1[0], matrix1[1], matrix2);
	}

	void transform(VectorVectorType& zs, SizeType offset)
	{
		VectorType gstmp(transform_.rows());

		const SizeType excitedPlusOne = zs.size();

		for (SizeType i = 0; i < excitedPlusOne; ++i)
			LanczosGlobals::transform(zs[i], offset, gstmp, transform_);
	}

	SizeType sectors() const { return 2; }

	void setPointer(SizeType p) { pointer_ = p; }

	PsimagLite::String name() const { return "reflection"; }

	void fullDiag(VectorRealType& eigs, MatrixType& fm) const
	{
		if (matrixStored_[pointer_].rows() > 1000)
			throw PsimagLite::RuntimeError("fullDiag too big\n");

		fm = matrixStored_[pointer_].toDense();
		diag(fm, eigs, 'V');

		if (!printMatrix_)
			return;

		for (SizeType i = 0; i < eigs.size(); i++)
			std::cout << eigs[i] << "\n";
		std::cout << fm;
	}

	template <typename SomeVectorType>
	void matrixVectorProduct(SomeVectorType& x, SomeVectorType const& y) const
	{
		return matrixStored_[pointer_].matrixVectorProduct(x, y);
	}

private:

	void addTo(WordType& yy, SizeType what, SizeType site) const
	{
		if (what == 0)
			return;
		WordType mask = (1 << site);
		yy |= mask;
	}

	void setTransform(const PsimagLite::Vector<ItemType>::Type& buffer2)
	{
		PsimagLite::Vector<ItemType>::Type buffer;
		makeUnique(buffer, buffer2);
		assert(buffer.size() == transform_.rows());
		SizeType counter      = 0;
		RealType oneOverSqrt2 = 1.0 / sqrt(2.0);
		RealType sign         = 1.0;
		SizeType row          = 0;
		for (SizeType i = 0; i < buffer.size(); i++) {
			if (buffer[i].type == ItemType::MINUS)
				continue;
			transform_.setRow(row++, counter);
			switch (buffer[i].type) {
			case ItemType::DIAGONAL:
				transform_.pushCol(buffer[i].i);
				transform_.pushValue(1);
				counter++;
				break;
			case ItemType::PLUS:
				transform_.pushCol(buffer[i].i);
				transform_.pushValue(oneOverSqrt2);
				counter++;
				transform_.pushCol(buffer[i].j);
				transform_.pushValue(oneOverSqrt2);
				counter++;
				break;
			}
		}

		for (SizeType i = 0; i < buffer.size(); i++) {
			if (buffer[i].type != ItemType::MINUS)
				continue;
			transform_.setRow(row++, counter);
			transform_.pushCol(buffer[i].i);
			transform_.pushValue(oneOverSqrt2 * sign);
			counter++;
			transform_.pushCol(buffer[i].j);
			transform_.pushValue(-oneOverSqrt2 * sign);
			counter++;
		}
		transform_.setRow(transform_.rows(), counter);
		transform_.checkValidity();
	}

	void makeUnique(PsimagLite::Vector<ItemType>::Type&       dest,
	                const PsimagLite::Vector<ItemType>::Type& src)
	{
		SizeType zeros   = 0;
		SizeType pluses  = 0;
		SizeType minuses = 0;
		for (SizeType i = 0; i < src.size(); i++) {
			ItemType item = src[i];
			int      x    = PsimagLite::indexOrMinusOne(dest, item);
			if (x >= 0)
				continue;
			if (item.type == ItemType::DIAGONAL)
				zeros++;
			if (item.type == ItemType::PLUS)
				pluses++;
			if (item.type == ItemType::MINUS)
				minuses++;

			dest.push_back(item);
		}
		PsimagLite::OstringStream msg(std::cout.precision());
		msg() << pluses << " +, " << minuses << " -, " << zeros << " zeros.";
		progress_.printline(msg, std::cout);
		plusSector_ = zeros + pluses;
	}

	void isIdentity(const SparseMatrixType& s, const PsimagLite::String& label) const
	{
		std::cerr << "Checking label=" << label << "\n";
		for (SizeType i = 0; i < s.rank(); i++) {
			for (int k = s.getRowPtr(i); k < s.getRowPtr(i + 1); k++) {
				SizeType col = s.getCol(k);
				RealType val = s.getValue(k);
				if (col == i)
					assert(isAlmostZero(val - 1.0));
				else
					assert(isAlmostZero(val));
			}
		}
	}

	bool isAlmostZero(const RealType& x) const { return (fabs(x) < 1e-6); }

	void split(SparseMatrixType&       matrixA,
	           SparseMatrixType&       matrixB,
	           const SparseMatrixType& matrix) const
	{
		SizeType counter = 0;
		matrixA.resize(plusSector_, plusSector_);
		for (SizeType i = 0; i < plusSector_; i++) {
			matrixA.setRow(i, counter);
			for (int k = matrix.getRowPtr(i); k < matrix.getRowPtr(i + 1); k++) {
				SizeType          col = matrix.getCol(k);
				ComplexOrRealType val = matrix.getValue(k);
				if (col < plusSector_) {
					matrixA.pushCol(col);
					matrixA.pushValue(val);
					counter++;
					continue;
				}
				if (PsimagLite::norm(val) > 1e-12) {
					PsimagLite::String s(__FILE__);
					s += " Hamiltonian has no reflection symmetry.";
					throw std::runtime_error(s.c_str());
				}
			}
		}
		matrixA.setRow(plusSector_, counter);

		SizeType rank        = matrix.rows();
		SizeType minusSector = rank - plusSector_;
		matrixB.resize(minusSector, minusSector);
		counter = 0;
		for (SizeType i = plusSector_; i < rank; i++) {
			matrixB.setRow(i - plusSector_, counter);
			for (int k = matrix.getRowPtr(i); k < matrix.getRowPtr(i + 1); k++) {
				SizeType          col = matrix.getCol(k);
				ComplexOrRealType val = matrix.getValue(k);
				if (col >= plusSector_) {
					matrixB.pushCol(col - plusSector_);
					matrixB.pushValue(val);
					counter++;
					continue;
				}

				if (PsimagLite::norm(val) > 1e-12) {
					PsimagLite::String s(__FILE__);
					s += " Hamiltonian has no reflection symmetry.";
					throw std::runtime_error(s.c_str());
				}
			}
		}
		matrixB.setRow(minusSector, counter);
	}

	PsimagLite::ProgressIndicator                       progress_;
	SparseMatrixType                                    transform_;
	SizeType                                            plusSector_;
	typename PsimagLite::Vector<SparseMatrixType>::Type matrixStored_;
	SizeType                                            pointer_;
	bool                                                printMatrix_;
}; // class ReflectionSymmetry
} // namespace Dmrg

#endif // REFLECTION_SYMM_H
