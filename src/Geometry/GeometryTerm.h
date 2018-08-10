/*
Copyright (c) 2009-2013, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]
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
/** \ingroup PsimagLite */
/*@{*/

/*! \file GeometryTerm.h
 *
 * Each geometry term represents a Hamiltonian connection term of the form
 * X_{ij} A_i B_j, where A and B are operators and X_{ij} are numbers.
 * Note that on-site Hamiltonian terms aren't included in the geometry since
 * they're trivially handled by the DMRG algorithm.
 */
#ifndef GEOMETRY_TERM_H
#define GEOMETRY_TERM_H

#include "GeometryDirection.h"
#include "GeometryBase.h"
#include <cassert>
#include "Ladder.h"
#include "LadderX.h"
#include "LadderBath.h"
#include "KTwoNiFFour.h"
#include "Star.h"
#include "LongChain.h"
#include "LongRange.h"
#include "Honeycomb.h"
#include "ExpressionCalculator.h"
#include "PsimagLite.h"

namespace PsimagLite {

template<typename ComplexOrRealType,typename InputType>
class GeometryTerm {

	typedef GeometryBase<ComplexOrRealType, InputType> GeometryBaseType;
	typedef GeometryDirection<ComplexOrRealType,GeometryBaseType> GeometryDirectionType;

public:

	typedef typename GeometryBaseType::AdditionalDataType AdditionalDataType;
	typedef typename Real<ComplexOrRealType>::Type RealType;

	struct Auxiliary {

		Auxiliary(bool d, SizeType t, SizeType n, SizeType l)
		    : debug(d), termId(t), numberOfTerms(n), linSize(l)
		{}

		void write(PsimagLite::String label, IoSerializer& ioSerializer) const
		{
			ioSerializer.createGroup(label);
			ioSerializer.write(label + "/debug", debug);
			ioSerializer.write(label + "/termId", termId);
			ioSerializer.write(label + "/numberOfTerms", numberOfTerms);
			ioSerializer.write(label + "/linSize", linSize);
		}

		bool debug;
		SizeType termId;
		SizeType numberOfTerms;
		SizeType linSize;
	}; // Auxiliary

	GeometryTerm()
	    : orbitals_(0),geometryBase_(0)
	{}

	/** @class hide_geometry2
	 - DegreesOfFreedom=integer Degrees of freedom on which the connectors depend on.
	 - GeometryKind=string One of chain, chainEx, ladder, ladderx, ladderbath, ktwoniffour,
	   or star.
	 - GeometryOptions=string Either none or ConstantValues needs to explain more FIXME
	*/
	GeometryTerm(InputType& io,
	             const Auxiliary& aux)
	    : aux_(aux),geometryBase_(0),gOptions_("none")
	{
		String savedPrefix = io.prefix();
		io.prefix() += (aux.numberOfTerms > 1) ? "gt" + ttos(aux.termId) + ":" : "";

		int x = -1;
		io.readline(x,   "DegreesOfFreedom=");
		if (x<=0) throw RuntimeError("DegreesOfFreedom<=0 is an error\n");

		SizeType edof = (x > 1);
		String s;
		io.readline(s,  "GeometryKind=");

		io.readline(gOptions_, "GeometryOptions=");
		bool constantValues = (gOptions_.find("ConstantValues") != String::npos);

		if (s == "chain" || s=="longchain") {
			geometryBase_ = new LongChain<ComplexOrRealType, InputType>(aux.linSize,io);
		} else if (s == "chainEx") {
			throw RuntimeError("GeometryTerm::ctor(): ChainEx: no longer supported.\n");
		} else if (s=="ladder") {
			geometryBase_ = new Ladder<ComplexOrRealType, InputType>(aux.linSize,io);
		} else if (s=="ladderx") {
			geometryBase_ = new LadderX<ComplexOrRealType, InputType>(aux.linSize,io);
		} else if (s=="ladderbath") {
			geometryBase_ = new LadderBath<ComplexOrRealType, InputType>(aux.linSize,io);
		} else if (s=="ktwoniffour") {
			geometryBase_ = new KTwoNiFFour<ComplexOrRealType, InputType>(aux.linSize,io);
		} else if (s=="star") {
			geometryBase_ = new Star<ComplexOrRealType, InputType>(aux.linSize,io);
		} else if (s=="LongRange") {
			geometryBase_ = new LongRange<ComplexOrRealType, InputType>(aux.linSize,io);
			edof |= 2;
		} else if (s == "Honeycomb") {
			geometryBase_ = new Honeycomb<ComplexOrRealType, InputType>(aux.linSize,io);
		} else {
			throw RuntimeError("Unknown geometry " + s + "\n");
		}

		for (SizeType i=0;i<geometryBase_->dirs();i++) {
			typename GeometryDirectionType::Auxiliary aux(constantValues, i, edof);

			directions_.push_back(GeometryDirectionType(io, aux, geometryBase_));
		}

		try {
			io.readline(vModifier_,  "GeometryValueModifier=");
		} catch (std::exception&) {}

		orbitals_ = findOrbitals();
		cacheValues();

		io.prefix() = savedPrefix;

		if (aux.debug) {
			std::cerr<<"Cached values:\n";
			std::cerr<<cachedValues_;
			std::cerr<<"-----------\n";
		}
	}

	~GeometryTerm()
	{
		if (geometryBase_) delete geometryBase_;
	}

	void write(PsimagLite::String label, IoSerializer& ioSerializer) const
	{
		ioSerializer.createGroup(label);
		aux_.write(label + "/aux_", ioSerializer);
		ioSerializer.write(label + "/orbitals_", orbitals_);
		// geometryBase_->write(label + "/geometryBase_", ioSerializer);
		ioSerializer.write(label + "/gOptions_", gOptions_);
		ioSerializer.write(label + "/vModifier_", vModifier_);
		ioSerializer.write(label + "/directions_", directions_);
		cachedValues_.write(label + "/cachedValues_", ioSerializer);
	}

	static String import()
	{
		String str("");

		for (SizeType i = 0; i < 9; ++i) {
			String istr = (i < 8) ? "gt" + ttos(i) + ":" : "";
			str += "integer " + istr + "DegreesOfFreedom;\n";
			str += "string " + istr + "GeometryKind;\n";
			str += "string " + istr + "GeometryOptions;\n";
			str += "integer " + istr + "LadderLeg;\n";
			str += "integer " + istr + "LongChainDistance;\n";
			for (SizeType j = 0; j < 9; ++j) {
				String jstr = "dir" + ttos(j) + ":";
				str += "vector " + istr + jstr + "Connectors;\n";
			}
		}

		return str;
	}

	template<class Archive>
	void write(Archive &, const unsigned int)
	{}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&,
	                   SizeType,
	                   String) const
	{
		return 0;
	}

	const ComplexOrRealType& operator()(SizeType i1,
	                                    SizeType edof1,
	                                    SizeType i2,
	                                    SizeType edof2) const
	{
		int k1 = geometryBase_->index(i1,edof1,orbitals_);
		int k2 = geometryBase_->index(i2,edof2,orbitals_);
		assert(k1>=0 && k2>=0);
		return cachedValues_(k1,k2);
	}

	template<typename T>
	typename EnableIf<IsComplexNumber<T>::True || Loki::TypeTraits<T>::isStdFloat,
	T>::Type vModifier(T value, RealType time) const
	{
		if (vModifier_ == "") return value;

		typedef ExpressionCalculator<T> ExpressionCalculatorType;
		typename ExpressionCalculatorType::VectorStringType ve;
		split(ve, vModifier_, ",");

		PrepassData<T> pd;
		typename PrepassData<T>::VectorType vr(2,0);
		vr[0] = time;
		vr[1] = value;
		pd.names = "tv";
		pd.values = vr;

		ExpressionPrepass<PrepassData<T> >::prepass(ve,pd);

		ExpressionCalculatorType ec(ve);
		return ec();
	}

	//assumes 1<smax+1 < emin
	const ComplexOrRealType& operator()(SizeType smax,
	                                    SizeType emin,
	                                    SizeType i1,
	                                    SizeType edof1,
	                                    SizeType i2,
	                                    SizeType edof2) const
	{
		bool bothFringe = (geometryBase_->fringe(i1,smax,emin) &&
		                   geometryBase_->fringe(i2,smax,emin));
		SizeType siteNew1 = i1;
		SizeType siteNew2 = i2;
		SizeType edofNew1 = edof1;
		SizeType edofNew2 = edof2;
		if (bothFringe) {
			if (i2<i1) {
				siteNew1 = i2;
				siteNew2 = i1;
				edofNew1 = edof2;
				edofNew2 = edof1;
			}
			siteNew2 = geometryBase_->getSubstituteSite(smax,emin,siteNew2);
		}
		return operator()(siteNew1,edofNew1,siteNew2,edofNew2);
	}

	bool connected(SizeType smax,SizeType emin,SizeType i1,SizeType i2) const
	{
		if (i1==i2) return false;

		bool bothFringe = (geometryBase_->fringe(i1,smax,emin) &&
		                   geometryBase_->fringe(i2,smax,emin));

		if (!bothFringe) return geometryBase_->connected(i1,i2);
		//std::cerr<<"fringe= "<<i1<<" "<<i2<<"\n";
		return true;
	}

	bool connected(SizeType i1,SizeType i2) const
	{
		return geometryBase_->connected(i1,i2);
	}

	String label() const
	{
		return geometryBase_->label();
	}

	SizeType maxConnections() const
	{
		return geometryBase_->maxConnections();
	}

	void fillAdditionalData(AdditionalDataType& additionalData,
	                        SizeType ind,
	                        SizeType jnd) const
	{
		geometryBase_->fillAdditionalData(additionalData,ind,jnd);
	}

	SizeType findReflection(SizeType site) const
	{
		return geometryBase_->findReflection(site);
	}

	SizeType length(SizeType i) const
	{
		return geometryBase_->length(i);
	}

	SizeType translate(SizeType site,SizeType dir,SizeType amount) const
	{
		return geometryBase_->translate(site,dir,amount);
	}

	void print(std::ostream& os) const
	{
		SizeType linSize = aux_.linSize;

		os<<"#orbital changes first\n";
		for (SizeType i=0;i<linSize;i++) {
			SizeType dofsi = orbitals(i);
			for (SizeType dof1=0;dof1<dofsi;dof1++) {
				for (SizeType j=0;j<linSize;j++) {
					SizeType dofsj = orbitals(j);
					for (SizeType dof2=0;dof2<dofsj;dof2++) {
						if (!connected(i,j)) {
							os<<0<<" ";
							continue;
						}
						os<<operator()(i,dof1,j,dof2)<<" ";
					}
				}

				os<<"\n";
			}
		}
	}

	SizeType handle(SizeType ind, SizeType jnd) const
	{
		return geometryBase_->handle(ind,jnd);
	}

	SizeType directions() const
	{
		return geometryBase_->dirs();
	}

	SizeType calcDir(SizeType i, SizeType j) const
	{
		return geometryBase_->calcDir(i,j);
	}

	String options() const
	{
		return gOptions_;
	}

	template<typename ComplexOrRealType_,typename InputType_>
	friend std::ostream& operator<<(std::ostream& os,
	                                const GeometryTerm<ComplexOrRealType_,
	                                InputType_>& gt);

	SizeType orbitals(SizeType site) const
	{
		return geometryBase_->orbitals(orbitals_, site);
	}

private:

	void cacheValues()
	{
		SizeType linSize = aux_.linSize;
		SizeType matrixRank = geometryBase_->matrixRank(linSize, orbitals_);
		cachedValues_.resize(matrixRank,matrixRank);

		for (SizeType i1 = 0; i1 < linSize; ++i1) {
			for (SizeType i2 = 0; i2 < linSize; ++i2) {
				if (!geometryBase_->connected(i1,i2)) continue;
				for (SizeType edof1=0;edof1<orbitals_;edof1++) {
					int k1 = geometryBase_->index(i1,edof1,orbitals_);
					if (k1<0) continue;
					for (SizeType edof2=0;edof2<orbitals_;edof2++) {
						int k2 = geometryBase_->index(i2,edof2,orbitals_);
						if (k2<0) continue;
						cachedValues_(k1,k2)=calcValue(i1,edof1,i2,edof2);
					}
				}
			}
		}
	}

	SizeType findOrbitals() const
	{
		if (directions_.size() == 0) {
			String str("GeometryTerm: ");
			str += "No directions found for this term.\n";
			throw RuntimeError(str);
		}

		SizeType orbitals = directions_[0].orbitals();
		for (SizeType dir=1;dir<directions_.size();dir++) {
			if (orbitals == directions_[dir].orbitals()) continue;
			String str("GeometryTerm: All directions must have");
			str += " connector matrices of the same rank\n";
			throw RuntimeError(str);
		}

		return orbitals;
	}

	ComplexOrRealType calcValue(SizeType i1,
	                            SizeType edof1,
	                            SizeType i2,
	                            SizeType edof2) const
	{
		if (!geometryBase_->connected(i1,i2)) return 0.0;

		SizeType dir = geometryBase_->calcDir(i1,i2);
		assert(dir<directions_.size());
		return directions_[dir](i1,edof1,i2,edof2);
	}

	GeometryTerm(const GeometryTerm&);

	GeometryTerm& operator=(const GeometryTerm&);

	Auxiliary aux_;
	SizeType orbitals_;
	GeometryBaseType* geometryBase_;
	String gOptions_;
	String vModifier_;
	typename Vector<GeometryDirectionType>::Type directions_;
	Matrix<ComplexOrRealType> cachedValues_;
}; // class GeometryTerm

template<typename ComplexOrRealType,typename InputType>
std::ostream& operator<<(std::ostream& os,
                         const GeometryTerm<ComplexOrRealType,InputType>& gt)
{
	os<<"#GeometryDirections="<<gt.directions_.size()<<"\n";
	for (SizeType i=0;i<gt.directions_.size();i++) os<<gt.directions_[i];
	gt.print(os);
	return os;
}
} // namespace PsimagLite

/*@}*/
#endif // GEOMETRY_TERM_H

