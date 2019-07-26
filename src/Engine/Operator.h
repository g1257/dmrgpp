/*
Copyright (c) 2009-2016-2018-2019, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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
data_, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, data_, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/

/*! \file Operator.h
 *
 *  A class to represent an operator
 *  Contains the actual data_, the (J,M) that indicates
 * how this operator transforms, the fermionSign which
 * indicates if this operator commutes or anticommutes
 * with operators of the same class on different sites, and
 * other properties.
 *
 */
#ifndef OPERATOR_H
#define OPERATOR_H
#include "Su2Related.h"
#include "CrsMatrix.h"
#include "Io/IoSelector.h"
#include "InputNg.h"
#include "InputCheck.h"
#include "CanonicalExpression.h"
#include "Io/IoSerializerStub.h"
#include "ProgramGlobals.h"
#include "OperatorStorage.h"

namespace Dmrg {

// This is a structure, don't add member functions here!
template<typename ComplexOrRealType>
class Operator {

public:

	enum {CAN_BE_ZERO = false, MUST_BE_NONZERO = true};

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef ComplexOrRealType value_type;
	typedef PsimagLite::CrsMatrix<value_type> SparseMatrixType;
	typedef typename PsimagLite::Real<value_type>::Type RealType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef Su2Related Su2RelatedType;
	typedef PsimagLite::Matrix<value_type> DenseMatrixType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef OperatorStorage<ComplexOrRealType> StorageType;

	Operator()
	    : fermionOrBoson_(ProgramGlobals::FermionOrBosonEnum::BOSON), angularFactor_(1)
	{}

	Operator(const SparseMatrixType& data1,
	         ProgramGlobals::FermionOrBosonEnum fermionSign1,
	         const PairType& jm1,
	         RealType angularFactor1,
	         const Su2RelatedType& su2Related1)
	    : data_(data1),
	      fermionOrBoson_(fermionSign1),
	      jm_(jm1),
	      angularFactor_(angularFactor1),
	      su2Related_(su2Related1)
	{}

	Operator(const StorageType& data1,
	         ProgramGlobals::FermionOrBosonEnum fermionSign1,
	         const PairType& jm1,
	         RealType angularFactor1,
	         const Su2RelatedType& su2Related1)
	    : data_(data1),
	      fermionOrBoson_(fermionSign1),
	      jm_(jm1),
	      angularFactor_(angularFactor1),
	      su2Related_(su2Related1)
	{}

	template<typename IoInputType, typename SomeModelType>
	Operator(IoInputType& io,
	         SomeModelType& model,
	         bool checkNonZero,
	         PsimagLite::String prefix)
	{
		/*PSIDOC Operator
		 \item[TSPOperator] [String] One of \{\verb!cooked!, \verb!raw!,
		 \verb!expression!\}, in order to indicate
		 how the operator will be specified.
		 \item[OperatorExpression=] A label containing an operator expression.
		 Specify only if \verb!TSPOperator! was set to \verb!expression!
		 \item[COOKED\_OPERATOR] [String] A label naming the operator. This is model
		 dependent and must be listed in the \verb!naturalOperator! function for
		 the indicated in \verb!Model! in this input file. Do not specify unless
		 \verb!TSPOperator!
		 was set to \verb!cooked!.
		 \item[COOKED\_EXTRA] [VectorInteger] The first number is the number of numbers
		 to follow. The other numbers are degrees of freedom for the cooked operator
		 mentioned, and are passed as arguments (in order)
		 to the \verb!naturalOperator! function for
		 the indicated in \verb!Model! in this input file. Do not specify unless
		 \verb!TSPOperator!
		 was set to \verb!cooked!.
		 \item[RAW\_MATRIX] [MatrixComplexOrRealType] The number of rows and
		 columns of this matrix, followed by the matrix in zig-zag format.
		 Do not specify unless
		 \verb!TSPOperator!
		 was set to \verb!raw!.
		 \item[FERMIONSIGN] [RealType] Either 1 or -1, indicating if this operator
		 commutes or anticommutes at \emph{different} sites. Do not specify if
		 \verb!TSPOperator!
		 was set to \verb!expression!.
		 \item[jm_VALUES] [Integer*2] If not using $SU(2)$ symmetry this is \verb!0 0!.
		 Else it is the $2j$ and $j+m$ for this operator. Do not specify if
		 \verb!TSPOperator!
		 was set to \verb!expression!.
		 \item[angularFactor_] [RealType] If not using $SU(2)$ symmetry this is \verb!1!.
		 Else FIXME. Do not specify if
		 \verb!TSPOperator!
		 was set to \verb!expression!.
		 */
		PsimagLite::String s = "";
		io.readline(s, prefix + "TSPOperator=");

		if (s == "cooked") {
			io.readline(s, prefix + "COOKED_OPERATOR=");
			VectorSizeType v;
			io.read(v,prefix + "COOKED_EXTRA");
			if (v.size() != 2)
				throw PsimagLite::RuntimeError("COOKED_EXTRA must be followed 2 v0 v1\n");
			data_ = model.naturalOperator(s,v[0],v[1]).data_;
		} else if (s == "raw") {
			DenseMatrixType m;
			io.read(m, prefix + "RAW_MATRIX");
			if (checkNonZero) checkNotZeroMatrix(m);
			fullMatrixToCrsMatrix(data_,m);
			PsimagLite::String msg = "WARNING: RAW_MATRIX read, order of basis subject ";
			msg += "to change with DMRG++ version!\n";
			std::cerr<<msg;
			std::cout<<msg;
		} else if (s == "expression") {
			io.readline(s,prefix + "OperatorExpression=");
			int site = 0;
			typedef OperatorSpec<SomeModelType, Operator> OperatorSpecType;
			OperatorSpecType opSpec(model);
			PsimagLite::CanonicalExpression<OperatorSpecType> canonicalExpression(opSpec);
			Operator p;
			const Operator opEmpty;
			canonicalExpression(p, s, opEmpty, site);
			data_ = p.data_;
			fermionOrBoson_ = p.fermionOrBoson_;
			jm_ = p.jm_;
			angularFactor_ = p.angularFactor_;
			// TODO FIXME: deprecate cooked
			return;
		} else {
			PsimagLite::String str(__FILE__);
			str += " : " + ttos(__LINE__) + "\n";
			str += "TSPOperator= must be followed by one of";
			str += "raw, cooked, or expression, not " + s + "\n";
			throw PsimagLite::RuntimeError(str.c_str());
		}

		int fs = 0;
		io.readline(fs,prefix + "FERMIONSIGN=");
		fermionOrBoson_ = (fs < 0) ? ProgramGlobals::FermionOrBosonEnum::FERMION
		                           : ProgramGlobals::FermionOrBosonEnum::BOSON;

		jm_.first = jm_.second = 0;
		angularFactor_ = 1;
		if (!SomeModelType::MyBasis::useSu2Symmetry())
			return;

		VectorSizeType v;
		io.read(v, prefix + "JMVALUES");
		if (v.size() != 2)
			err("FATAL: JMVALUES is not a vector of 2 values\n");
		jm_.first = v[0]; jm_.second = v[1];

		io.readline(angularFactor_,prefix + "angularFactor_=");

		// FIXME: su2Related_ needs to be set properly for when SU(2) is running
	}

	void fromStorage(const PsimagLite::Matrix<ComplexOrRealType>& m)
	{
		data_.fromDense(m);
	}

	void fromStorage(const SparseMatrixType& m)
	{
		fromCRS(data_, m);
	}

	void write(PsimagLite::String label,
	           PsimagLite::IoSerializer& ioSerializer,
	           PsimagLite::IoSerializer::WriteMode mode =
	        PsimagLite::IoNgSerializer::NO_OVERWRITE) const
	{
		if (mode != PsimagLite::IoNgSerializer::ALLOW_OVERWRITE)
			ioSerializer.createGroup(label);

		data_.write(label + "/data_", ioSerializer, mode);
		ioSerializer.write(label + "/fermionOrBoson_", fermionOrBoson_, mode);
		ioSerializer.write(label + "/jm_", jm_, mode);
		ioSerializer.write(label + "/angularFactor_", angularFactor_, mode);
		// su2Related_.write(label + "/su2Related_", ioSerializer);
	}

	void overwrite(PsimagLite::String label,
	               PsimagLite::IoSerializer& ioSerializer) const
	{
		write(label, ioSerializer, PsimagLite::IoNgSerializer::ALLOW_OVERWRITE);
	}

	void read(PsimagLite::String label,
	          PsimagLite::IoSerializer& ioSerializer)
	{
		data_.read(label + "/data_", ioSerializer);
		ioSerializer.read(fermionOrBoson_, label + "/fermionOrBoson_");
		ioSerializer.read(jm_, label + "/jm_");
		ioSerializer.read(angularFactor_, label + "/angularFactor_");
		// su2Related_.read(label + "/su2Related_", ioSerializer);
	}

	void write(std::ostream& os) const
	{
		os<<"TSPOperator=raw\n";
		os<<"RAW_MATRIX\n";
		DenseMatrixType m;
		crsMatrixToFullMatrix(m,data_);
		os<<m;
		int fs = (fermionOrBoson_ == ProgramGlobals::FermionOrBosonEnum::FERMION) ? -1 : 1;
		os<<"FERMIONSIGN="<<fs <<"\n";
		os<<"JMVALUES 2 "<<jm_.first<<" "<<jm_.second<<"\n";
		os<<"angularFactor_="<<angularFactor_<<"\n";
	}

	// operators START
	Operator operator*(const Operator& other) const
	{
		const ProgramGlobals::FermionOrBosonEnum f = ProgramGlobals::multipy(fermionOrBoson_,
		                                                                     other.fermionOrBoson_);
		return Operator(data_*other.data_, f, PairType(0, 0), 1.0, su2Related_);
	}

	Operator& operator*=(value_type x)
	{
		data_ *= x;
		return *this;
	}

	Operator& operator*=(const Operator& other)
	{
		int fSaved = (fermionOrBoson_ == ProgramGlobals::FermionOrBosonEnum::FERMION) ? -1 : 1;
		fermionOrBoson_ = other.fermionOrBoson_;
		if (metaDiff(other) > 0)
			err("operator+= failed for Operator: metas not equal\n");

		StorageType crs = data_*other.data_;
		data_ = crs;

		int fsOther = (other.fermionOrBoson_ == ProgramGlobals::FermionOrBosonEnum::FERMION) ? -1
		                                                                                     : 1;
		int fs = fSaved * fsOther;
		fermionOrBoson_ = (fs < 0) ? ProgramGlobals::FermionOrBosonEnum::FERMION
		                           : ProgramGlobals::FermionOrBosonEnum::BOSON;

		return *this;
	}

	Operator& operator+=(const Operator& other)
	{
		if (metaDiff(other) > 0)
			err("operator+= failed for Operator: metas not equal\n");
		data_ += other.data_;
		return *this;
	}

	// operators END

	void outerProduct(const Operator& A,
	                  SizeType nout,
	                  const VectorRealType& signs,
	                  bool order,
	                  const VectorSizeType& permutationFull)
	{
		externalProduct2(data_,
		                 A.getStorage(),
		                 nout,
		                 signs,
		                 order,
		                 permutationFull);
	}

	SizeType metaDiff(const Operator& op2) const
	{
		const Operator& op1 = *this;

		SizeType code = 0;
		PsimagLite::Vector<bool>::Type b(4, false);

		b[0] = (op1.fermionOrBoson_ != op2.fermionOrBoson_);
		b[1] = (op1.angularFactor_ != op2.angularFactor_);
		b[2] = (op1.jm_ != op2.jm_);
		//b[3] = (op1.su2Related_ != op2.su2Related_);

		SizeType orFactor = 0;
		for (SizeType i = 0; i < b.size(); ++i) {
			if (b[i]) code |= orFactor;
			orFactor <<= 1;
		}

		return code;
	}

	bool isEmpty() const
	{
		return (data_.rows() == 0);
	}

	void clear()
	{
		data_.clear();
	}

	void set(ProgramGlobals::FermionOrBosonEnum fOrB,
	         PairType jm1,
	         RealType af,
	         const Su2RelatedType& su2)
	{
		fermionOrBoson_ = fOrB;
		jm_ = jm1;
		angularFactor_ = af;
		su2Related_ = su2;
	}

	void set(ProgramGlobals::FermionOrBosonEnum fOrB,
	         PairType jm1,
	         RealType af)
	{
		fermionOrBoson_ = fOrB;
		jm_ = jm1;
		angularFactor_ = af;
	}

	// short cut
	const SparseMatrixType& getCRS() const { return data_.getCRS(); }

	const StorageType& getStorage() const { return data_; }

	// FIXME TODO
	StorageType& getStorageNonConst() { return data_; }

	const ProgramGlobals::FermionOrBosonEnum& fermionOrBoson() const
	{
		return fermionOrBoson_;
	}

	const PairType& jm() const { return jm_; }

	const RealType& angularFactor() const { return angularFactor_; }

	const Su2RelatedType& su2Related() const { return su2Related_; }

	// FIXME TODO
	Su2RelatedType& su2RelatedNonConst() { return su2Related_; }

	void conjugate()
	{
		data_.conjugate();
	}

	void transpose()
	{
		data_.transpose();
	}

	void dagger()
	{
		StorageType copy = data_;
		transposeConjugate(data_, copy);
	}

	friend void reorder2(Operator& v, const VectorSizeType& permutation)
	{
		reorder2(v.data_, permutation);
	}

	friend void bcast2(Operator& op)
	{
		bcast(op.data_);
		PsimagLite::MPI::bcast(op.fermionOrBoson_);
		PsimagLite::MPI::bcast(op.jm_);
		PsimagLite::MPI::bcast(op.angularFactor_);
		bcast(op.su2Related_);
	}

private:

	void checkNotZeroMatrix(const DenseMatrixType& m) const
	{
		RealType norma = norm2(m);
		RealType eps = 1e-6;
		if (norma>eps) return;

		PsimagLite::String s(__FILE__);
		s += " : " + ttos(__LINE__) + "\n";
		s += "RAW_MATRIX or COOKED_OPERATOR ";
		s += " is less than " + ttos(eps) + "\n";
		std::cerr<<"WARNING: "<<s;
	}

	StorageType data_;
	// does this operator commute or anticommute with others of the
	// same class on different sites
	ProgramGlobals::FermionOrBosonEnum fermionOrBoson_;
	PairType  jm_; // angular momentum of this operator
	RealType angularFactor_;
	Su2RelatedType su2Related_;
};

template<typename SparseMatrixType,
         template<typename,typename> class SomeVectorTemplate,
         typename SomeAllocator1Type,
         typename SomeAllocator2Type>
void fillOperator(SomeVectorTemplate<SparseMatrixType*,SomeAllocator1Type>& data_,
                  SomeVectorTemplate<Operator<SparseMatrixType>,SomeAllocator2Type>& op)
{
	for (SizeType i = 0; i < data_.size(); ++i)
		data_[i] = &(op[i].data_);
}

template<typename SparseMatrixType>
std::istream& operator>>(std::istream& is,Operator<SparseMatrixType>& op)
{
	is>>op.data_;
	is>>op.fermionSign;
	SizeType theNumber2 = 0;
	is>>theNumber2;
	is>>op.jm_;
	is>>op.angularFactor_;
	is>>op.su2Related_;
	return is;
}

template<typename SparseMatrixType>
std::ostream& operator<<(std::ostream& os,const Operator<SparseMatrixType>& op)
{
	os<<op.data_;
	os<<op.fermionSign<<"\n";
	os<<"2\n"<<op.jm_.first<<" "<<op.jm_.second<<"\n";
	os<<op.angularFactor_<<"\n";
	os<<op.su2Related_;
	return os;
}
} // namespace Dmrg

/*@}*/
#endif

