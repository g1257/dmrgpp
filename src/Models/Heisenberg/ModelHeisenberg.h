/*
Copyright (c) 2009, 2017, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/
/** \ingroup DMRG */
/*@{*/

/*! \file ModelHeisenberg.h
 *
 *  An implementation of the Quantum Heisenberg Model to use with  DmrgSolver
 *
 */

#ifndef DMRG_MODEL_HEISENBERG_HEADER_H
#define DMRG_MODEL_HEISENBERG_HEADER_H

#include <algorithm>
#include "ModelBase.h"
#include "ParametersModelHeisenberg.h"
#include "LinkProductHeisenberg.h"
#include "CrsMatrix.h"
#include "VerySparseMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "ProgramGlobals.h"
#include "ModelCommon.h"
#include "Utils.h"

namespace Dmrg {

template<typename ModelBaseType>
class ModelHeisenberg : public ModelBaseType {

public:

	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::BasisType BasisType;
	typedef typename ModelBaseType::GeometryType GeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkProductStructType LinkProductStructType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelBaseType::TargetQuantumElectronsType TargetQuantumElectronsType;
	typedef	typename ModelBaseType::VectorType VectorType;

private:

	typedef typename ModelBaseType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef unsigned int long WordType;
	typedef LinkProductHeisenberg<ModelHelperType> LinkProductType;
	typedef ModelCommon<ModelBaseType,LinkProductType> ModelCommonType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelBaseType::VectorRealType VectorRealType;

	static const int NUMBER_OF_ORBITALS=1;
	static const int DEGREES_OF_FREEDOM=2; // spin up and down

public:

	typedef typename PsimagLite::Vector<unsigned int long>::Type HilbertBasisType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename MyBasis::SymmetryElectronsSzType SymmetryElectronsSzType;

	ModelHeisenberg(const SolverParamsType& solverParams,
	                InputValidatorType& io,
	                const GeometryType& geometry,
	                PsimagLite::String additional)
	    : ModelBaseType(io,new ModelCommonType(solverParams,geometry)),
	      modelParameters_(io),
	      geometry_(geometry),
	      spinSquared_(spinSquaredHelper_,NUMBER_OF_ORBITALS,DEGREES_OF_FREEDOM)
	{
		if (additional == "Anisotropic") {
			setAnisotropic(solverParams);
		}

		SizeType n = geometry_.numberOfSites();
		SizeType m = modelParameters_.magneticField.size();
		SizeType md = modelParameters_.anisotropy.size();

		if (m > 0 && m != n) {
			PsimagLite::String msg("ModelHeisenberg: If provided, ");
			msg += " MagneticField must be a vector of " + ttos(n) + " entries.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (md > 0 && md != n) {
			PsimagLite::String msg("ModelHeisenberg: If provided, ");
			msg += " Anisotropy must be a vector of " + ttos(n) + " entries.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (BasisType::useSu2Symmetry() && modelParameters_.twiceTheSpin != 1) {
			PsimagLite::String msg("ModelHeisenberg: SU(2) symmetry, ");
			msg += " for spin different than 1/2 is not implemented yet.\n";
			throw PsimagLite::RuntimeError(msg);
		}
	}

	SizeType memResolv(PsimagLite::MemResolv& mres,
	                   SizeType,
	                   PsimagLite::String msg = "") const
	{
		PsimagLite::String str = msg;
		str += "ModelHeisenberg";

		const char* start = reinterpret_cast<const char *>(this);
		const char* end = reinterpret_cast<const char *>(&modelParameters_);
		SizeType total = end - start;
		mres.push(PsimagLite::MemResolv::MEMORY_TEXTPTR,
		          total,
		          start,
		          msg + " ModelHeisenberg vptr");

		start = end;
		end = start + PsimagLite::MemResolv::SIZEOF_HEAPPTR;
		total += mres.memResolv(&modelParameters_, end-start, str + " modelParameters");

		start = end;
		end = reinterpret_cast<const char *>(&spinSquaredHelper_);
		mres.push(PsimagLite::MemResolv::MEMORY_HEAPPTR,
		          PsimagLite::MemResolv::SIZEOF_HEAPREF,
		          start,
		          str + " ref to geometry");
		total += (end - start);
		mres.memResolv(&geometry_, 0, str + " geometry");

		start = end;
		end = reinterpret_cast<const char *>(&spinSquared_);
		total += mres.memResolv(&spinSquaredHelper_,
		                        end-start,
		                        str + " spinSquaredHelper");

		total += mres.memResolv(&spinSquared_,
		                        sizeof(*this) - total,
		                        str + " spinSquared");

		return total;
	}

	void print(std::ostream& os) const { os<<modelParameters_; }

	SizeType hilbertSize(SizeType) const
	{
		return modelParameters_.twiceTheSpin + 1;
	}

	//! set operator matrices for sites in block
	void setOperatorMatrices(VectorOperatorType& operatorMatrices,
	                         const BlockType& block) const
	{
		HilbertBasisType natBasis;
		SparseMatrixType tmpMatrix;

		SymmetryElectronsSzType qq;
		setBasis(natBasis, qq, block);

		operatorMatrices.clear();
		for (SizeType i=0;i<block.size();i++) {
			// Set the operators S^+_i in the natural basis
			tmpMatrix=findSplusMatrices(i,natBasis);

			typename OperatorType::Su2RelatedType su2related;
			su2related.source.push_back(i*DEGREES_OF_FREEDOM);
			su2related.source.push_back(i*DEGREES_OF_FREEDOM+NUMBER_OF_ORBITALS);
			su2related.source.push_back(i*DEGREES_OF_FREEDOM);
			su2related.transpose.push_back(-1);
			su2related.transpose.push_back(-1);
			su2related.transpose.push_back(1);
			su2related.offset = NUMBER_OF_ORBITALS;

			OperatorType myOp(tmpMatrix,1,PairType(2,2),-1,su2related);
			operatorMatrices.push_back(myOp);

			// Set the operators S^z_i in the natural basis
			tmpMatrix = findSzMatrices(i,natBasis);
			typename OperatorType::Su2RelatedType su2related2;
			OperatorType myOp2(tmpMatrix,1,PairType(2,1),1.0/sqrt(2.0),su2related2);
			operatorMatrices.push_back(myOp2);

			if (LinkProductType::terms() == 2) continue;
			// Set the operators S^x_i in the natural basis
			tmpMatrix = findSxMatrices(i,natBasis);
			typename OperatorType::Su2RelatedType su2related3;
			OperatorType myOp3(tmpMatrix,1,PairType(2,1),1.0/sqrt(2.0),su2related3);
			operatorMatrices.push_back(myOp3);
		}
	}

	OperatorType naturalOperator(const PsimagLite::String& what,
	                             SizeType site,
	                             SizeType) const
	{
		BlockType block;
		block.resize(1);
		block[0]=site;
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		setOperatorMatrices(creationMatrix,block);
		assert(creationMatrix.size()>0);

		if (what=="splus") { // S^+
			return creationMatrix[0];
		}

		if (what=="sminus") { // S^-
			creationMatrix[0].conjugate();
			return creationMatrix[0];
		}

		if (what == "z" || what == "sz") { // S^z
			return creationMatrix[1];
		}

		PsimagLite::String str("ModelHeisenberg: naturalOperator: no label ");
		str += what + "\n";
		throw PsimagLite::RuntimeError(str);
	}

	//! Dummy since this model has no fermion sign
	void findElectrons(VectorSizeType& electrons,
	                   const HilbertBasisType& basis,
	                   SizeType) const
	{
		electrons.resize(basis.size());
		for (SizeType i=0;i<electrons.size();i++)
			electrons[i] = 0;
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const VectorOperatorType& cm,
	                                const BlockType& block,
	                                RealType,
	                                RealType factorForDiagonals=1.0)  const
	{
		SizeType linSize = geometry_.numberOfSites();
		SizeType n = block.size();

		if (modelParameters_.magneticField.size() == linSize) {
			for (SizeType i = 0; i < n; ++i) {
				// magnetic field
				RealType tmp = modelParameters_.magneticField[block[i*2]]*factorForDiagonals;
				hmatrix += tmp*cm[1+i*2].data;
			}
		}

		if (modelParameters_.anisotropy.size() != linSize)
			return; // <--- PLEASE NOTE EARLY EXIT HERE

		SparseMatrixType Szsquare;
		for (SizeType i = 0; i < n; ++i) {
			// anisotropy
			RealType tmp = modelParameters_.anisotropy[block[i*2]]*factorForDiagonals;
			multiply(Szsquare,cm[1+i*2].data,cm[1+i*2].data);
			hmatrix += tmp*Szsquare;
		}

	}

	virtual const TargetQuantumElectronsType& targetQuantum() const
	{
		return modelParameters_.targetQuantum;
	}

	//! find all states in the natural basis for a block of n sites
	void setBasis(HilbertBasisType& basis,
	              SymmetryElectronsSzType& qq,
	              const VectorSizeType& block) const
	{
		SizeType total = utils::powUint(modelParameters_.twiceTheSpin + 1,block.size());

		HilbertBasisType basisTmp(total);
		for (SizeType i = 0; i < total; ++i) basisTmp[i] = i;
		setSymmetryRelated(qq, basisTmp, block.size());
		ModelBaseType::orderBasis(basis, basisTmp, qq);
	}

private:

	void setAnisotropic(const SolverParamsType& solverParams)
	{
		LinkProductType::setAnisotropic();
		bool isCanonical = (modelParameters_.targetQuantum.isCanonical);
		bool useTheForce = (solverParams.options.find("useTheForce") !=
		        PsimagLite::String::npos);
		if (!isCanonical) return;

		PsimagLite::String warning("HeisenbergAnisotropic: ");
		warning += "canonical mode in use. ";
		warning += "Results will likely be WRONG.\n";
		warning += "Please delete the TargetSzPlusConst= ";
		warning += "line in the input file.\n";

		if (useTheForce) {
			std::cerr<<"WARNING: "<<warning;
			std::cout<<"WARNING: "<<warning;
			return;
		}

		std::cerr<<"FATAL: "<<warning;
		err("You may useTheForce in SolverOptions to run it anyway\n");
	}

	SizeType logBase2(SizeType x) const
	{
		SizeType counter = 0;
		while (x > 0) {
			x >>= 1;
			counter++;
		}

		return (counter == 0) ? counter : counter - 1;
	}

	//! Find S^+_site in the natural basis natBasis
	SparseMatrixType findSplusMatrices(SizeType site,
	                                   const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total,total);
		RealType j = 0.5*modelParameters_.twiceTheSpin;
		SizeType bitsForOneSite = utils::bitSizeOfInteger(modelParameters_.twiceTheSpin);
		SizeType bits = 1 + logBase2(modelParameters_.twiceTheSpin);
		SizeType mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;
		mask <<= (site*bitsForOneSite);

		for (SizeType ii=0;ii<total;ii++) {
			SizeType ket = natBasis[ii];

			SizeType ketsite = ket & mask;
			ketsite >>= (site*bitsForOneSite);

			assert(ketsite == ket);
			SizeType brasite = ketsite + 1;
			if (brasite >= modelParameters_.twiceTheSpin+1) continue;

			SizeType bra = ket & (~mask);
			assert(bra == 0);
			brasite <<= (site*bitsForOneSite);
			bra |= brasite;
			assert(bra == brasite);

			RealType m = ketsite - j;
			RealType x = j*(j+1)-m*(m+1);
			assert(x>=0);

			cm(ket,bra) = sqrt(x);
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Find S^z_i in the natural basis natBasis
	SparseMatrixType findSzMatrices(SizeType site,
	                                const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total,total);
		RealType j = 0.5*modelParameters_.twiceTheSpin;
		SizeType bitsForOneSite = utils::bitSizeOfInteger(modelParameters_.twiceTheSpin);
		SizeType bits = logBase2(modelParameters_.twiceTheSpin) + 1;
		SizeType mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;
		mask <<= (site*bitsForOneSite);

		for (SizeType ii=0;ii<total;ii++) {
			SizeType ket = natBasis[ii];

			SizeType ketsite = ket & mask;
			ketsite >>= (site*bitsForOneSite);
			assert(ketsite == ket);
			RealType m = ketsite - j;
			cm(ket,ket) = m;
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	SparseMatrixType findSxMatrices(SizeType site,
	                                const HilbertBasisType& natBasis) const
	{
		SparseMatrixType Splus_temp=findSplusMatrices(site,natBasis);
		SparseMatrixType Sminus_temp,Sx;
		transposeConjugate(Sminus_temp,Splus_temp);
		RealType tmp=0.5;

		Sx = tmp*Splus_temp;
		Sx += tmp*Sminus_temp;

		return Sx;
	}

	void setSymmetryRelated(SymmetryElectronsSzType& q,
	                        const HilbertBasisType& basis,
	                        int n) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;
		typename PsimagLite::Vector<PairType>::Type jmvalues;
		VectorSizeType flavors;
		PairType jmSaved;
		jmSaved.first = modelParameters_.twiceTheSpin;
		jmSaved.second = basis[0];
		jmSaved.first++;
		jmSaved.second++;

		SizeType totalElectrons = modelParameters_.targetQuantum.totalElectrons;
		if (!modelParameters_.targetQuantum.isSu2 && totalElectrons > 0) {
			PsimagLite::String msg("Please delete the line ");
			msg += "TargetElectronsTotal= in the input file\n";
			throw PsimagLite::RuntimeError(msg);
		}

		bool isCanonical = (modelParameters_.targetQuantum.isCanonical);
		VectorSizeType szPlusConst(basis.size());
		VectorSizeType electrons(basis.size());
		for (SizeType i=0;i<basis.size();i++) {
			PairType jmpair;
			jmpair.first = modelParameters_.twiceTheSpin;
			jmpair.second = basis[i];
			jmvalues.push_back(jmpair);
			szPlusConst[i] = (isCanonical) ? getSzPlusConst(basis[i],n) : 0;
			electrons[i] = (modelParameters_.targetQuantum.isSu2) ? 1 : 0;
			flavors.push_back(1);
			jmSaved = jmpair;
		}

		q.set(jmvalues,flavors,electrons,szPlusConst);
	}

	SizeType getSzPlusConst(SizeType ket, SizeType n) const
	{
		if (n == 1) return ket;

		SizeType bitsForOneSite = utils::bitSizeOfInteger(modelParameters_.twiceTheSpin);

		SizeType sum = 0;
		for (SizeType site = 0; site < n; ++site) {
			SizeType mask = modelParameters_.twiceTheSpin;
			mask <<= (site*bitsForOneSite);
			SizeType ketsite = ket & mask;
			ketsite >>= (site*bitsForOneSite);
			sum += ketsite;
		}

		return sum;
	}

	//serializr start class ModelHeisenberg
	//serializr vptr
	//serializr normal modelParameters_
	ParametersModelHeisenberg<RealType>  modelParameters_;
	//serializr ref geometry_ start
	GeometryType const &geometry_;
	//serializr normal spinSquaredHelper_
	SpinSquaredHelper<RealType,WordType> spinSquaredHelper_;
	//serializr normal spinSquared_
	SpinSquared<SpinSquaredHelper<RealType,WordType> > spinSquared_;
}; // class ModelHeisenberg

} // namespace Dmrg
/*@}*/
#endif //DMRG_MODEL_HEISENBERG_HEADER_H

