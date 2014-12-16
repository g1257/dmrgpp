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
#include "String.h"
#include "Ladder.h"
#include "LadderX.h"
#include "LadderBath.h"
#include "KTwoNiFFour.h"
#include "Star.h"
#include "LongChain.h"

namespace PsimagLite {

template<typename ComplexOrRealType,typename InputType>
class GeometryTerm {

	typedef GeometryBase<InputType> GeometryBaseType;
	typedef GeometryDirection<ComplexOrRealType,GeometryBaseType> GeometryDirectionType;

	enum {NUMBERS = GeometryDirectionType::NUMBERS,
	      MATRICES = GeometryDirectionType::MATRICES};

public:

	typedef typename GeometryBaseType::AdditionalDataType AdditionalDataType;
	typedef typename Real<ComplexOrRealType>::Type RealType;

	GeometryTerm()
	    : linSize_(0),maxEdof_(0),geometryBase_(0)
	{}

	/** @class hide_geometry2
	 - DegreesOfFreedom=integer Degrees of freedom on which the connectors depend on.
	 - GeometryKind=string One of chain, chainEx, ladder, ladderx, ladderbath, ktwoniffour,
	   or star.
	 - GeometryOptions=string Either none or ConstantValues needs to explain more FIXME
	*/
	GeometryTerm(InputType& io,SizeType,SizeType linSize,bool debug=false) :
	    linSize_(linSize),maxEdof_(0),geometryBase_(0),gOptions_("none")
	{
		int x;
		io.readline(x,"DegreesOfFreedom=");
		if (x<=0) throw RuntimeError("DegreesOfFreedom<=0 is an error\n");

		SizeType edof = (x==1) ? NUMBERS : MATRICES;
		String s;
		io.readline(s,"GeometryKind=");

		io.readline(gOptions_,"GeometryOptions=");

		if (s == "chain" || s=="longchain") {
			geometryBase_ = new LongChain<InputType>(linSize,io);
		} else if (s == "chainEx") {
			throw RuntimeError("GeometryTerm::ctor(): ChainEx: no longer supported.\n");
		} else if (s=="ladder") {
			geometryBase_ = new Ladder<InputType>(linSize,io);
		} else if (s=="ladderx") {
			geometryBase_ = new LadderX<InputType>(linSize,io);
		} else if (s=="ladderbath") {
			geometryBase_ = new LadderBath<InputType>(linSize,io);
		} else if (s=="ktwoniffour") {
			geometryBase_ = new KTwoNiFFour<InputType>(linSize,io);
		} else if (s=="star") {
			geometryBase_ = new Star<InputType>(linSize,io);
		} else {
			throw RuntimeError("Unknown geometry " + s + "\n");
		}

		for (SizeType i=0;i<geometryBase_->dirs();i++) {
			directions_.push_back(GeometryDirectionType(io,
			                                            i,
			                                            edof,
			                                            gOptions_,
			                                            geometryBase_));
		}

		findMaxEdof();
		cacheValues();

		if (debug) {
			std::cerr<<"Cached values:\n";
			std::cerr<<cachedValues_;
			std::cerr<<"-----------\n";
		}
	}

	~GeometryTerm()
	{
		if (geometryBase_) delete geometryBase_;
	}

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar.template register_type<LongChain<InputType> >();
		ar.template register_type<Ladder<InputType> >();
		ar.template register_type<LadderX<InputType> >();
		ar.template register_type<LadderBath<InputType> >();
		ar.template register_type<KTwoNiFFour<InputType> >();
		ar.template register_type<Star<InputType> >();
		ar & linSize_;
		ar & maxEdof_;
		ar & geometryBase_;
		ar & directions_;
		ar & cachedValues_;
	}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType x,
	                   PsimagLite::String msg) const
	{
		PsimagLite::String str = msg;
		str += "GeometryTerm";
		const char* start = (const char *)&linSize_;
		const char* end = (const char*)&maxEdof_;
		SizeType total = mres.memResolv(&linSize_,end-start,str + " linSize");

		start = end;
		end = (const char*)&geometryBase_;
		total += mres.memResolv(&maxEdof_,end-start,str + " maxEdof");

		start = end;
		end = (const char*)&directions_;
		total += (end - start);
		mres.push(SomeMemResolvType::MEMORY_HEAPPTR,
		                               sizeof(geometryBase_),
			                           &geometryBase_,
			                           str + " geometryBase");

		start = end;
		end = (const char*)&cachedValues_;
		total += mres.memResolv(&directions_,end-start,str + " directions");

		total += mres.memResolv(&cachedValues_,sizeof(*this)-total, str+" cachedValues");

		geometryBase_->memResolv(mres,x,str);

		return total;
	}

	const ComplexOrRealType& operator()(SizeType i1,
	                                    SizeType edof1,
	                                    SizeType i2,
	                                    SizeType edof2) const
	{
		int k1 = geometryBase_->index(i1,edof1,maxEdof_);
		int k2 = geometryBase_->index(i2,edof2,maxEdof_);
		assert(k1>=0 && k2>=0);
		return cachedValues_(k1,k2);
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
		SizeType linSize = linSize_;
		SizeType dofs = 1;
		for (SizeType dof1=0;dof1<dofs;dof1++) {
			for (SizeType dof2=0;dof2<dofs;dof2++) {
				os<<"dof1="<<dof1<<" dof2="<<dof2<<"\n";
				for (SizeType i=0;i<linSize;i++) {
					for (SizeType j=0;j<linSize;j++) {
						if (!connected(i,j)) {
							os<<0<<" ";
							continue;
						}
						os<<operator()(i,dof1,j,dof2)<<" ";
					}
					os<<"\n";
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

	PsimagLite::String options() const
	{
		return gOptions_;
	}

	template<typename ComplexOrRealType_,typename InputType_>
	friend std::ostream& operator<<(std::ostream& os,
	                                const GeometryTerm<ComplexOrRealType_,
	                                                   InputType_>& gt);

private:

	void cacheValues()
	{
		SizeType matrixRank = geometryBase_->matrixRank(linSize_,maxEdof_);
		cachedValues_.resize(matrixRank,matrixRank);

		for (SizeType i1=0;i1<linSize_;i1++) {
			for (SizeType i2=0;i2<linSize_;i2++) {
				if (!geometryBase_->connected(i1,i2)) continue;
				for (SizeType edof1=0;edof1<maxEdof_;edof1++) {
					int k1 = geometryBase_->index(i1,edof1,maxEdof_);
					if (k1<0) continue;
					for (SizeType edof2=0;edof2<maxEdof_;edof2++) {
						int k2 = geometryBase_->index(i2,edof2,maxEdof_);
						if (k2<0) continue;
						cachedValues_(k1,k2)=calcValue(i1,edof1,i2,edof2);
					}
				}
			}
		}
	}

	void findMaxEdof()
	{
		maxEdof_ = 0;
		for (SizeType dir=0;dir<directions_.size();dir++) {
			maxEdof_ = directions_[dir].nRow();
			if (maxEdof_< directions_[dir].nCol())
				maxEdof_ = directions_[dir].nCol();
		}
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

	SizeType linSize_;
	SizeType maxEdof_;
	GeometryBaseType* geometryBase_;
	PsimagLite::String gOptions_;
	typename Vector<GeometryDirectionType>::Type directions_;
	PsimagLite::Matrix<ComplexOrRealType> cachedValues_;
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

